#nullable enable
using System;
using System.Collections.Generic;

namespace TopDownSimulator.Comparison;

/// <summary>
/// How closely a simulated centroid spectrum resembles a real one over one m/z range.
/// </summary>
/// <param name="RealPeaks">Real centroids in the range.</param>
/// <param name="SimulatedPeaks">Simulated centroids in the range.</param>
/// <param name="RealMatchedFraction">
/// Fraction of real peaks with a simulated partner within tolerance. Low means the simulation is
/// missing peaks the instrument reported.
/// </param>
/// <param name="SimulatedMatchedFraction">
/// Fraction of simulated peaks with a real partner within tolerance. Low means the simulation is
/// reporting peaks the instrument did not, which is how spurious peaks between isotopologues show up.
/// </param>
/// <param name="Cosine">Cosine over matched pairs, unmatched peaks paired with zero. Dominated by the tallest peaks.</param>
/// <param name="SqrtCosine">
/// The same cosine on square-rooted intensities, which gives the many faint peaks a say. This is the
/// one that moves when the floor between isotopologues is wrong.
/// </param>
/// <param name="TicRatio">Simulated over real summed intensity. 1 means the absolute scale agrees.</param>
/// <param name="RelativeIntensityKs">
/// Kolmogorov-Smirnov distance between the two distributions of log10(intensity / that spectrum's
/// own maximum in the range). Scale-free, so it isolates the dynamic-range shape: 0 is identical,
/// 1 is disjoint.
/// </param>
/// <param name="RealFractionBelowOnePercent">Fraction of real peaks under 1 % of the range's real maximum.</param>
/// <param name="SimulatedFractionBelowOnePercent">The same for the simulated spectrum.</param>
public sealed record CentroidSpectrumComparison(
    int RealPeaks,
    int SimulatedPeaks,
    double RealMatchedFraction,
    double SimulatedMatchedFraction,
    double Cosine,
    double SqrtCosine,
    double TicRatio,
    double RelativeIntensityKs,
    double RealFractionBelowOnePercent,
    double SimulatedFractionBelowOnePercent)
{
    /// <summary>Simulated over real peak count.</summary>
    public double PeakCountRatio => RealPeaks > 0 ? SimulatedPeaks / (double)RealPeaks : double.NaN;

    /// <summary>
    /// Compares two centroid spectra over [<paramref name="minMz"/>, <paramref name="maxMz"/>].
    /// Both m/z arrays must be ascending.
    /// </summary>
    /// <remarks>
    /// Peaks are paired one-to-one by a two-pointer pass that takes the closest available partner
    /// within <paramref name="ppmTolerance"/>. The default of 10 ppm covers the measured per-scan
    /// common-mode offset (σ 1.9 ppm) plus per-peak scatter at modest S/N; much wider and unrelated
    /// peaks in a dense spectrum start pairing by chance.
    /// </remarks>
    public static CentroidSpectrumComparison Compare(
        double[] realMz, double[] realIntensity,
        double[] simulatedMz, double[] simulatedIntensity,
        double minMz = double.NegativeInfinity, double maxMz = double.PositiveInfinity,
        double ppmTolerance = 10.0)
    {
        ArgumentNullException.ThrowIfNull(realMz);
        ArgumentNullException.ThrowIfNull(realIntensity);
        ArgumentNullException.ThrowIfNull(simulatedMz);
        ArgumentNullException.ThrowIfNull(simulatedIntensity);
        if (realMz.Length != realIntensity.Length || simulatedMz.Length != simulatedIntensity.Length)
            throw new ArgumentException("Each m/z array must be the same length as its intensity array.");
        if (!(ppmTolerance > 0))
            throw new ArgumentOutOfRangeException(nameof(ppmTolerance), ppmTolerance, "Tolerance must be positive.");

        var (rMz, rI) = Slice(realMz, realIntensity, minMz, maxMz);
        var (sMz, sI) = Slice(simulatedMz, simulatedIntensity, minMz, maxMz);

        var (realPartner, simMatched) = Match(rMz, sMz, ppmTolerance);

        double dot = 0, rNorm = 0, sNorm = 0;
        double sqrtDot = 0, rSqrtNorm = 0, sSqrtNorm = 0;
        for (int i = 0; i < rMz.Length; i++)
        {
            double r = rI[i];
            double s = realPartner[i] >= 0 ? sI[realPartner[i]] : 0;
            dot += r * s;
            sqrtDot += Math.Sqrt(r) * Math.Sqrt(s);
            rNorm += r * r;
            rSqrtNorm += r;
        }

        for (int j = 0; j < sMz.Length; j++)
        {
            sNorm += sI[j] * sI[j];
            sSqrtNorm += sI[j];
        }

        int realMatched = 0;
        foreach (int p in realPartner)
            if (p >= 0) realMatched++;

        double rTic = Sum(rI), sTic = Sum(sI);

        return new CentroidSpectrumComparison(
            RealPeaks: rMz.Length,
            SimulatedPeaks: sMz.Length,
            RealMatchedFraction: rMz.Length > 0 ? realMatched / (double)rMz.Length : double.NaN,
            SimulatedMatchedFraction: sMz.Length > 0 ? simMatched / (double)sMz.Length : double.NaN,
            Cosine: Ratio(dot, rNorm, sNorm),
            SqrtCosine: Ratio(sqrtDot, rSqrtNorm, sSqrtNorm),
            TicRatio: rTic > 0 ? sTic / rTic : double.NaN,
            RelativeIntensityKs: KolmogorovSmirnov(RelativeLog(rI), RelativeLog(sI)),
            RealFractionBelowOnePercent: FractionBelow(rI, 0.01),
            SimulatedFractionBelowOnePercent: FractionBelow(sI, 0.01));
    }

    private static double Ratio(double dot, double aNorm, double bNorm)
    {
        if (aNorm <= 0 && bNorm <= 0) return 1;
        if (aNorm <= 0 || bNorm <= 0) return 0;
        return dot / Math.Sqrt(aNorm * bNorm);
    }

    private static (double[] Mz, double[] Intensity) Slice(double[] mz, double[] intensity, double min, double max)
    {
        int lo = LowerBound(mz, min);
        int hi = lo;
        while (hi < mz.Length && mz[hi] <= max) hi++;
        return (mz[lo..hi], intensity[lo..hi]);
    }

    private static int LowerBound(double[] sorted, double value)
    {
        int lo = 0, hi = sorted.Length;
        while (lo < hi)
        {
            int mid = (lo + hi) >>> 1;
            if (sorted[mid] < value) lo = mid + 1;
            else hi = mid;
        }

        return lo;
    }

    /// <summary>
    /// One-to-one pairing. Each real peak takes the nearest still-unclaimed simulated peak within
    /// tolerance; because both lists are ascending, a claimed partner never needs revisiting.
    /// </summary>
    private static (int[] RealPartner, int SimMatched) Match(double[] realMz, double[] simMz, double ppm)
    {
        var partner = new int[realMz.Length];
        Array.Fill(partner, -1);
        int next = 0, matched = 0;

        for (int i = 0; i < realMz.Length; i++)
        {
            double tol = realMz[i] * ppm * 1e-6;
            while (next < simMz.Length && simMz[next] < realMz[i] - tol) next++;

            int best = -1;
            double bestDelta = double.PositiveInfinity;
            for (int j = next; j < simMz.Length && simMz[j] <= realMz[i] + tol; j++)
            {
                double delta = Math.Abs(simMz[j] - realMz[i]);
                if (delta < bestDelta)
                {
                    bestDelta = delta;
                    best = j;
                }
            }

            if (best < 0) continue;

            // A later real peak may be closer to this candidate, but taking it here keeps the pass
            // linear; in a centroid list two real peaks within 10 ppm of one simulated peak are rare.
            partner[i] = best;
            matched++;
            next = best + 1;
        }

        return (partner, matched);
    }

    private static double[] RelativeLog(double[] intensity)
    {
        double max = 0;
        foreach (double v in intensity) max = Math.Max(max, v);

        var values = new List<double>(intensity.Length);
        if (max > 0)
            foreach (double v in intensity)
                if (v > 0) values.Add(Math.Log10(v / max));

        var array = values.ToArray();
        Array.Sort(array);
        return array;
    }

    private static double KolmogorovSmirnov(double[] a, double[] b)
    {
        if (a.Length == 0 || b.Length == 0) return double.NaN;

        int i = 0, j = 0;
        double d = 0;
        while (i < a.Length && j < b.Length)
        {
            double x = Math.Min(a[i], b[j]);
            while (i < a.Length && a[i] <= x) i++;
            while (j < b.Length && b[j] <= x) j++;
            d = Math.Max(d, Math.Abs(i / (double)a.Length - j / (double)b.Length));
        }

        return d;
    }

    private static double FractionBelow(double[] intensity, double fractionOfMax)
    {
        if (intensity.Length == 0) return double.NaN;
        double max = 0;
        foreach (double v in intensity) max = Math.Max(max, v);

        int below = 0;
        foreach (double v in intensity)
            if (v < fractionOfMax * max) below++;

        return below / (double)intensity.Length;
    }

    private static double Sum(double[] values)
    {
        double total = 0;
        foreach (double v in values) total += v;
        return total;
    }
}
