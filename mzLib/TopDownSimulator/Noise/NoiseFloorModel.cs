#nullable enable
using System;
using System.Collections.Generic;

namespace TopDownSimulator.Noise;

/// <summary>
/// The unstructured FT noise floor: how many spurious centroids an Orbitrap reports per scan, where
/// in m/z they land, and how bright they are.
/// </summary>
/// <remarks>
/// <para>
/// Both curves below are <b>measured</b>, not derived. The analytic prediction — noise peaks
/// uniformly distributed in frequency, so density ∝ (m/z)^−1.5 because f ∝ (m/z)^−1/2 — is wrong by
/// two orders of magnitude against real data: the observed density rises from m/z 600 to a maximum
/// near 750 and then falls by a factor of ~500 out to m/z 2000, far steeper than any power law,
/// because what gets *reported* is local maxima above a detection threshold rather than every
/// occupied frequency bin. Measuring the curve reproduces reality by construction; deriving it does
/// not. See <c>Test/TopDownSimulator/NoiseCharacterization.cs</c> for the measurement.
/// </para>
/// <para>
/// Calibrated against the 50–55 min window (post-elution, where the floor is unobscured) of four
/// Jurkat top-down runs: rep1 fract6, rep1 fract7, rep2 fract5, rep2 fract7. The density curve and
/// the <i>shape</i> of the level curve agreed across all four to within a few percent. The
/// <i>amplitude</i> of the level curve did not — it spans 663 to 2050 across those files — so it is
/// the one free parameter, <see cref="NoiseLevelAtReferenceMz"/>.
/// </para>
/// <para>
/// The noise is genuinely unstructured, which the same measurement establishes rather than assumes:
/// in the quiet window the instrument failed to assign a charge to 98.4 % of peaks, the fraction
/// with a +1.00335/z partner was ~45 % <i>identically for every z from 1 to 7</i> (i.e. chance
/// collision at this peak density, not isotopic structure), and the brightest 10 % of peaks
/// persisted to the next scan no more often than random ones. In the busy window that last number
/// jumps to 77 %, which is what real analyte looks like. So this model emits independent singlets
/// per scan, and there is no persistent chemical-background component.
/// </para>
/// </remarks>
public sealed class NoiseFloorModel
{
    /// <summary>Width of every bin in the calibration table.</summary>
    private const double BinWidth = 100.0;

    /// <summary>Lower edge of the first bin in the calibration table.</summary>
    private const double FirstBinStart = 600.0;

    /// <summary>
    /// Reported noise peaks per scan in each 100-Th bin from m/z 600, averaged over the four
    /// calibration files. Sums to ~16 900 peaks per scan, against a measured median of 16 400–17 600.
    /// </summary>
    private static readonly double[] JurkatPeaksPerScanPerBin =
    {
        3527, 4303, 3504, 2519, 1337, 745, 383, 235, 141, 87.5, 50.8, 29.0, 16.4, 8.2,
    };

    /// <summary>Number of 100-Th bins in the density table, from m/z 600.</summary>
    public static int BinCount => JurkatPeaksPerScanPerBin.Length;

    /// <summary>Lower m/z edge of density bin <paramref name="bin"/>.</summary>
    public static double BinStart(int bin) => FirstBinStart + bin * BinWidth;

    /// <summary>
    /// Noise amplitude at <see cref="ReferenceMz"/> times injection time, in intensity·ms, for the
    /// Jurkat runs.
    /// </summary>
    /// <remarks>
    /// Reported intensities are normalised by injection time, so a fixed noise charge shows up as an
    /// amplitude proportional to 1/IT. Over rep2 fract7 this product stays between 1.5e4 and 2.5e4
    /// wherever AGC limits the fill, while the amplitude itself spans 180 to 49 000. It drops to
    /// about 9e3 only in empty scans at the maximum injection time.
    /// </remarks>
    public const double JurkatNoiseTimesInjectionTime = 2.0e4;

    private readonly double[] _peaksPerBin;

    /// <summary>
    /// Noise amplitude in each bin relative to the first, averaged over the four calibration files.
    /// Rises ~4.6× to a plateau near m/z 1450–1650 and then falls away.
    /// </summary>
    private static readonly double[] RelativeNoiseLevelPerBin =
    {
        1.00, 1.47, 2.04, 2.68, 3.26, 3.77, 4.04, 4.47, 4.63, 4.57, 4.55, 4.07, 3.82, 3.32,
    };

    /// <summary>Geometric mean of the measured noise amplitude at m/z 650 across the four files.</summary>
    public const double JurkatNoiseLevelAtReferenceMz = 1120.0;

    /// <summary>Centre of the first calibration bin, where <see cref="NoiseLevelAtReferenceMz"/> applies.</summary>
    public const double ReferenceMz = 650.0;

    private readonly double[] _binCentres;

    /// <summary>
    /// Lognormal σ of the noise amplitude at a <i>fixed</i> m/z, which the calibration table's
    /// per-bin medians do not capture.
    /// </summary>
    /// <remarks>
    /// Measured directly: pooling peaks in 5-Th slices and dividing each by its own slice median
    /// gives a right-skewed multiplicative spread with p25 ≈ 0.72, p75 ≈ 1.58, p95 ≈ 2.8, i.e.
    /// σ_log ≈ 0.6. Feeding that number in verbatim over-widens the simulated intensity
    /// distribution, because within a slice the local noise estimate and the S/N a peak achieves are
    /// negatively correlated and this model treats them as independent. The default below is
    /// therefore fitted to reproduce the measured <i>intensity</i> quantile ratios, which is what a
    /// feature finder actually sees, rather than to the amplitude spread in isolation.
    /// </remarks>
    public const double JurkatLevelDispersion = 0.42;

    /// <param name="peaksPerScanPerBin">
    /// Expected noise peaks per scan in each of the <see cref="BinCount"/> 100-Th bins from m/z 600,
    /// or null for the Jurkat quiet-window calibration. A per-scan table is how
    /// <see cref="ScanNoiseConditions"/> makes the density follow the run.
    /// </param>
    public NoiseFloorModel(
        double noiseLevelAtReferenceMz = JurkatNoiseLevelAtReferenceMz,
        double densityScale = 1.0,
        SignalToNoiseDistribution? signalToNoise = null,
        double levelDispersion = JurkatLevelDispersion,
        double[]? peaksPerScanPerBin = null)
    {
        if (peaksPerScanPerBin is not null)
        {
            if (peaksPerScanPerBin.Length != BinCount)
                throw new ArgumentException($"A density table needs {BinCount} bins.", nameof(peaksPerScanPerBin));
            foreach (double n in peaksPerScanPerBin)
                if (!(n >= 0) || !double.IsFinite(n))
                    throw new ArgumentOutOfRangeException(nameof(peaksPerScanPerBin), n,
                        "Bin densities must be finite and non-negative.");
        }

        _peaksPerBin = (double[])(peaksPerScanPerBin ?? JurkatPeaksPerScanPerBin).Clone();

        if (!(noiseLevelAtReferenceMz > 0) || !double.IsFinite(noiseLevelAtReferenceMz))
            throw new ArgumentOutOfRangeException(nameof(noiseLevelAtReferenceMz), noiseLevelAtReferenceMz,
                "Noise level must be finite and positive.");
        if (!(densityScale >= 0) || !double.IsFinite(densityScale))
            throw new ArgumentOutOfRangeException(nameof(densityScale), densityScale,
                "Density scale must be finite and non-negative.");
        if (!(levelDispersion >= 0) || !double.IsFinite(levelDispersion))
            throw new ArgumentOutOfRangeException(nameof(levelDispersion), levelDispersion,
                "Level dispersion must be finite and non-negative.");

        NoiseLevelAtReferenceMz = noiseLevelAtReferenceMz;
        DensityScale = densityScale;
        LevelDispersion = levelDispersion;
        SignalToNoise = signalToNoise ?? SignalToNoiseDistribution.Orbitrap;

        _binCentres = new double[_peaksPerBin.Length];
        for (int i = 0; i < _binCentres.Length; i++)
            _binCentres[i] = FirstBinStart + (i + 0.5) * BinWidth;
    }

    /// <summary>
    /// Noise amplitude at <see cref="ReferenceMz"/>, in the same intensity units as the simulated
    /// signal. The only parameter that has to be recalibrated per instrument or per run; the
    /// measured spread across the four Jurkat files was a factor of 3.
    /// </summary>
    public double NoiseLevelAtReferenceMz { get; }

    /// <summary>
    /// Multiplier on the number of noise peaks emitted. 1.0 is the measured density of ~16 900 peaks
    /// per scan, which is realistic but expensive: it roughly doubles the peak count of a
    /// thousand-proteoform simulation and takes a full-run mzML into the hundreds of megabytes, in
    /// line with the ~548 MB a real converted Jurkat run occupies. Lower it to iterate quickly.
    /// </summary>
    public double DensityScale { get; }

    /// <summary>The S/N a reported noise peak carries. Its intensity is this times the local level.</summary>
    public SignalToNoiseDistribution SignalToNoise { get; }

    /// <summary>
    /// Lognormal σ of the peak-to-peak variation in noise amplitude at a fixed m/z. Zero makes every
    /// peak at a given m/z share one amplitude, which measurably under-disperses the intensities.
    /// See <see cref="JurkatLevelDispersion"/>.
    /// </summary>
    public double LevelDispersion { get; }

    /// <summary>Lowest m/z the calibration covers.</summary>
    public double MinMz => FirstBinStart;

    /// <summary>Highest m/z the calibration covers.</summary>
    public double MaxMz => FirstBinStart + _peaksPerBin.Length * BinWidth;

    /// <summary>Expected noise peaks per scan in density bin <paramref name="bin"/>, before <see cref="DensityScale"/>.</summary>
    public double PeaksPerScanInBin(int bin) => _peaksPerBin[bin];

    /// <summary>
    /// A copy with a different amplitude and, optionally, density table. Dispersion, the S/N
    /// distribution and <see cref="DensityScale"/> are kept.
    /// </summary>
    public NoiseFloorModel With(double noiseLevelAtReferenceMz, double[]? peaksPerScanPerBin = null) =>
        new(noiseLevelAtReferenceMz, DensityScale, SignalToNoise, LevelDispersion, peaksPerScanPerBin ?? _peaksPerBin);

    /// <summary>Expected number of noise peaks in one scan, over the whole calibrated range.</summary>
    public double ExpectedPeaksPerScan
    {
        get
        {
            double total = 0;
            foreach (double n in _peaksPerBin)
                total += n;
            return total * DensityScale;
        }
    }

    /// <summary>
    /// Noise amplitude at <paramref name="mz"/>, linearly interpolated between bin centres and held
    /// flat outside the calibrated range. This is the σ against which a peak's S/N is measured, so
    /// it is also the natural detection threshold for signal peaks — see <see cref="NoiseRelativeFloor"/>
    /// and <see cref="Simulation.ScanReductionOptions.Floor"/>.
    /// </summary>
    public double NoiseLevelAt(double mz)
    {
        double relative = Interpolate(_binCentres, RelativeNoiseLevelPerBin, mz);
        return NoiseLevelAtReferenceMz * relative;
    }

    /// <summary>
    /// Draws one scan's worth of noise peaks into <paramref name="into"/>, in ascending m/z.
    /// </summary>
    /// <remarks>
    /// Sampled bin by bin: the count in each bin is Poisson about its calibrated mean and the
    /// positions within a bin are uniform, which is exact for a piecewise-constant density and
    /// avoids inverting a CDF. Peaks are independent between scans, as the persistence measurement
    /// says they should be.
    /// </remarks>
    public void SampleScan(Random rng, List<SimulatedPeak> into)
    {
        if (rng is null) throw new ArgumentNullException(nameof(rng));
        if (into is null) throw new ArgumentNullException(nameof(into));

        for (int bin = 0; bin < _peaksPerBin.Length; bin++)
        {
            double expected = _peaksPerBin[bin] * DensityScale;
            if (expected <= 0) continue;

            int count = SamplePoisson(rng, expected);
            if (count == 0) continue;

            double binStart = FirstBinStart + bin * BinWidth;
            int firstIndex = into.Count;
            for (int i = 0; i < count; i++)
            {
                double mz = binStart + rng.NextDouble() * BinWidth;
                double level = NoiseLevelAt(mz);
                if (LevelDispersion > 0)
                    level *= Math.Exp(LevelDispersion * SampleStandardNormal(rng));

                double intensity = level * SignalToNoise.Sample(rng);
                into.Add(new SimulatedPeak(mz, intensity, PeakOrigin.FtNoise));
            }

            // Sorting each bin's slice keeps the whole list ascending, since bins are visited in
            // ascending m/z and never overlap. Cheaper than one sort over all ~17 000 peaks.
            into.Sort(firstIndex, count, MzComparer.Instance);
        }
    }

    private sealed class MzComparer : IComparer<SimulatedPeak>
    {
        public static readonly MzComparer Instance = new();
        public int Compare(SimulatedPeak a, SimulatedPeak b) => a.Mz.CompareTo(b.Mz);
    }

    /// <summary>
    /// Poisson draw. Knuth's product method below λ=30 and a rounded normal approximation above it,
    /// where the product method would both underflow and cost O(λ) uniforms per draw.
    /// </summary>
    private static int SamplePoisson(Random rng, double lambda)
    {
        if (lambda < 30)
        {
            double limit = Math.Exp(-lambda);
            double product = rng.NextDouble();
            int count = 0;
            while (product > limit)
            {
                count++;
                product *= rng.NextDouble();
            }

            return count;
        }

        double normal = SampleStandardNormal(rng);
        return Math.Max(0, (int)Math.Round(lambda + Math.Sqrt(lambda) * normal));
    }

    internal static double SampleStandardNormal(Random rng)
    {
        // Box-Muller. One of the two variates is discarded; at these call rates that is cheaper than
        // carrying the cached-spare state through a method that is also called from parallel scans.
        double u1 = 1.0 - rng.NextDouble();
        double u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }

    /// <summary>Piecewise-linear interpolation over ascending <paramref name="xs"/>, clamped at both ends.</summary>
    private static double Interpolate(double[] xs, double[] ys, double x)
    {
        if (x <= xs[0]) return ys[0];
        if (x >= xs[^1]) return ys[^1];

        int hi = Array.BinarySearch(xs, x);
        if (hi >= 0) return ys[hi];
        hi = ~hi;

        int lo = hi - 1;
        double weight = (x - xs[lo]) / (xs[hi] - xs[lo]);
        return ys[lo] + weight * (ys[hi] - ys[lo]);
    }
}
