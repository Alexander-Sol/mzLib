#nullable enable
using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using Chemistry;
using TopDownSimulator.Model;

namespace TopDownSimulator.Simulation;

/// <summary>
/// Turns the forward model into centroid spectra the way an instrument does: one centroid per local
/// maximum of the summed profile, at the height of that maximum.
/// </summary>
/// <remarks>
/// <para>
/// The alternative this replaces evaluated the summed profile at the union of every proteoform's
/// theoretical isotopologue m/z and wrote every sample. That is exact for an isolated envelope and
/// badly wrong for a crowded one: a single real peak is sampled at every nearby isotopologue position
/// of every other species, its shoulders survive as peaks between the isotopologues, and anything
/// downstream that merges nearby centroids adds the same profile height up several times. Measured
/// on rep2 fract7 (H2B, m/z 860-864) it produced 256 peaks where the instrument reported 72.
/// </para>
/// <para>
/// Each scan is handled independently. The (proteoform, charge, isotopologue) Gaussians active at
/// that retention time are sorted by centre and split wherever two neighbours are more than
/// <see cref="InteractionWindowInSigmas"/> σ apart, since beyond that neither can move the other's
/// maximum. Each cluster is rendered on a grid of σ/<see cref="GridPointsPerSigma"/>, and each grid
/// maximum is polished by Newton's method on the analytic profile, so the reported position is the
/// true maximum and the reported height is the forward model evaluated there.
/// </para>
/// <para>
/// Thread-safe: everything taken from <see cref="IsotopeEnvelopeKernel"/> is copied into immutable
/// arrays in the constructor, so scans can be processed in parallel.
/// </para>
/// </remarks>
public sealed class ProfileCentroider
{
    /// <summary>Gaussians further apart than this cannot shift each other's maxima measurably.</summary>
    private const double InteractionWindowInSigmas = 6.0;

    /// <summary>Grid density. Two resolvable maxima are at least ~2σ apart, so this cannot skip one.</summary>
    private const int GridPointsPerSigma = 4;

    /// <summary>Same cut-off <see cref="ForwardModel"/> applies to a species' weight in a scan.</summary>
    private const double MinimumContributionWeight = 1e-18;

    private const int MaxNewtonIterations = 20;

    /// <summary>One (proteoform, charge) envelope, as arrays that are safe to share across threads.</summary>
    private sealed record Envelope(
        int Proteoform, double ChargeWeight, double[] Centers, double[] UnitHeights, double[] Inv2Sigma2, double[] Sigmas);

    private readonly ProteoformModel[] _proteoforms;
    private readonly Envelope[] _envelopes;
    private readonly IPeakWidthModel _widthModel;

    public ProfileCentroider(
        IReadOnlyList<ProteoformModel> proteoforms, int minCharge, int maxCharge, IPeakWidthModel widthModel)
    {
        if (proteoforms is null) throw new ArgumentNullException(nameof(proteoforms));
        if (minCharge < 1 || maxCharge < minCharge)
            throw new ArgumentException("Charge range must satisfy 1 ≤ minCharge ≤ maxCharge.");
        _widthModel = widthModel ?? throw new ArgumentNullException(nameof(widthModel));
        _proteoforms = new ProteoformModel[proteoforms.Count];

        var envelopes = new List<Envelope>();
        for (int p = 0; p < proteoforms.Count; p++)
        {
            var model = proteoforms[p];
            _proteoforms[p] = model;
            var kernel = new IsotopeEnvelopeKernel(model.MonoisotopicMass);
            int n = kernel.IsotopologueCount;

            for (int z = minCharge; z <= maxCharge; z++)
            {
                double fz = model.ChargeDistribution.Evaluate(z);
                if (!(fz > 0) || n == 0) continue;

                var centers = new double[n];
                var heights = new double[n];
                var inv2Sigma2 = new double[n];
                var sigmas = new double[n];
                for (int i = 0; i < n; i++)
                {
                    centers[i] = kernel.NeutralMass(i).ToMz(z);
                    double sigma = widthModel.SigmaAt(centers[i]);
                    if (!(sigma > 0) || !double.IsFinite(sigma))
                        throw new ArgumentOutOfRangeException(nameof(widthModel),
                            $"σ_m must be finite and positive; got {sigma} at m/z {centers[i]}.");

                    sigmas[i] = sigma;
                    // Identical to IsotopeEnvelopeKernel.Evaluate: unit-area Gaussians weighted by
                    // the normalized isotopologue abundance.
                    heights[i] = kernel.Intensity(i) / (sigma * Math.Sqrt(2 * Math.PI));
                    inv2Sigma2[i] = 1.0 / (2 * sigma * sigma);
                }

                envelopes.Add(new Envelope(p, fz, centers, heights, inv2Sigma2, sigmas));
            }
        }

        _envelopes = envelopes.ToArray();
    }

    /// <summary>Centroids every scan, in parallel. Output is parallel to <paramref name="scanTimes"/>.</summary>
    public (double[] Mz, double[] Intensity)[] Centroid(double[] scanTimes)
    {
        if (scanTimes is null) throw new ArgumentNullException(nameof(scanTimes));
        var result = new (double[] Mz, double[] Intensity)[scanTimes.Length];
        Parallel.For(0, scanTimes.Length, s => result[s] = Centroid(scanTimes[s]));
        return result;
    }

    /// <summary>The centroids of the summed profile at one retention time, in ascending m/z.</summary>
    public (double[] Mz, double[] Intensity) Centroid(double scanTime)
    {
        var components = ActiveComponents(scanTime);
        var mz = new List<double>();
        var intensity = new List<double>();
        if (components.Length == 0)
            return (Array.Empty<double>(), Array.Empty<double>());

        int start = 0;
        for (int k = 1; k <= components.Length; k++)
        {
            bool split = k == components.Length
                || components[k].Center - components[k - 1].Center
                   > InteractionWindowInSigmas * Math.Max(components[k].Sigma, components[k - 1].Sigma);
            if (!split) continue;

            PickCluster(components, start, k, mz, intensity);
            start = k;
        }

        return (mz.ToArray(), intensity.ToArray());
    }

    private readonly record struct Component(double Center, double Height, double Inv2Sigma2, double Sigma);

    private Component[] ActiveComponents(double scanTime)
    {
        var weights = new double[_proteoforms.Length];
        for (int p = 0; p < _proteoforms.Length; p++)
        {
            double rt = _proteoforms[p].RtProfile.Evaluate(scanTime);
            weights[p] = rt > 0 ? _proteoforms[p].Abundance * rt : 0;
        }

        var components = new List<Component>();
        foreach (var e in _envelopes)
        {
            double weight = weights[e.Proteoform];
            if (weight <= MinimumContributionWeight) continue;
            weight *= e.ChargeWeight;
            if (weight <= MinimumContributionWeight) continue;

            for (int i = 0; i < e.Centers.Length; i++)
            {
                double height = weight * e.UnitHeights[i];
                if (height > 0)
                    components.Add(new Component(e.Centers[i], height, e.Inv2Sigma2[i], e.Sigmas[i]));
            }
        }

        var array = components.ToArray();
        Array.Sort(array, static (a, b) => a.Center.CompareTo(b.Center));
        return array;
    }

    /// <summary>Finds the maxima of the profile made by components [<paramref name="from"/>, <paramref name="to"/>).</summary>
    private void PickCluster(Component[] all, int from, int to, List<double> mzOut, List<double> intensityOut)
    {
        // A lone Gaussian peaks exactly at its centre; no grid needed.
        if (to - from == 1)
        {
            mzOut.Add(all[from].Center);
            intensityOut.Add(all[from].Height);
            return;
        }

        double maxSigma = 0;
        for (int k = from; k < to; k++)
            maxSigma = Math.Max(maxSigma, all[k].Sigma);
        double window = InteractionWindowInSigmas * maxSigma;

        double lo = all[from].Center - all[from].Sigma;
        double hi = all[to - 1].Center + all[to - 1].Sigma;

        double prev2 = double.NaN, prev1 = double.NaN, prevX = double.NaN;
        double lastAcceptedMz = double.NegativeInfinity;
        for (double x = lo; ; )
        {
            double y = Evaluate(all, from, to, x, window).S;

            // prev1 is a maximum when it rises from prev2 and does not fall short of y.
            if (!double.IsNaN(prev2) && prev1 > prev2 && prev1 >= y)
            {
                var (peakMz, peakHeight) = Polish(all, from, to, prevX, window, _widthModel.SigmaAt(prevX) / GridPointsPerSigma);
                if (peakMz - lastAcceptedMz > 1e-9 * peakMz)
                {
                    mzOut.Add(peakMz);
                    intensityOut.Add(peakHeight);
                    lastAcceptedMz = peakMz;
                }
                else if (peakHeight > intensityOut[^1])
                {
                    intensityOut[^1] = peakHeight;
                }
            }

            if (x >= hi) break;
            prev2 = prev1;
            prev1 = y;
            prevX = x;
            x = Math.Min(hi, x + _widthModel.SigmaAt(x) / GridPointsPerSigma);
        }
    }

    /// <summary>
    /// Newton's method on S'(x) = 0 from a grid maximum, with each step held inside one grid cell so
    /// it cannot jump to a neighbouring maximum.
    /// </summary>
    private static (double Mz, double Height) Polish(Component[] all, int from, int to, double x0, double window, double step)
    {
        double x = x0;
        for (int iteration = 0; iteration < MaxNewtonIterations; iteration++)
        {
            var (_, d1, d2) = Evaluate(all, from, to, x, window);
            if (!(d2 < 0)) break;

            double delta = Math.Clamp(-d1 / d2, -step, step);
            double next = Math.Clamp(x + delta, x0 - step, x0 + step);
            if (Math.Abs(next - x) < 1e-12 * x)
            {
                x = next;
                break;
            }

            x = next;
        }

        return (x, Evaluate(all, from, to, x, window).S);
    }

    /// <summary>The profile and its first two derivatives at <paramref name="x"/>.</summary>
    private static (double S, double D1, double D2) Evaluate(Component[] all, int from, int to, double x, double window)
    {
        int k = LowerBound(all, from, to, x - window);
        double s = 0, d1 = 0, d2 = 0;
        for (; k < to && all[k].Center <= x + window; k++)
        {
            var c = all[k];
            double d = x - c.Center;
            double g = c.Height * Math.Exp(-d * d * c.Inv2Sigma2);
            double twoQ = 2 * c.Inv2Sigma2;
            s += g;
            d1 -= g * twoQ * d;
            d2 += g * (twoQ * twoQ * d * d - twoQ);
        }

        return (s, d1, d2);
    }

    private static int LowerBound(Component[] all, int from, int to, double value)
    {
        int lo = from, hi = to;
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            if (all[mid].Center < value) lo = mid + 1;
            else hi = mid;
        }

        return lo;
    }
}
