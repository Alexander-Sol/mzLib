#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using MassSpectrometry;
using MzLibUtil;
using TopDownSimulator.Model;
using TopDownSimulator.Simulation;

namespace TopDownSimulator.Noise;

/// <summary>Peak counts by provenance, summed over a whole simulation.</summary>
public sealed record NoiseInjectionSummary(
    long SignalPeaks, long NoisePeaks, long MergedPeaks, long DroppedByJitter = 0)
{
    /// <summary>Signal peaks that survived jitter, plus noise peaks, before merging.</summary>
    public long TotalPeaks => SignalPeaks + NoisePeaks;
}

/// <summary>
/// Applies the measurement error an instrument imposes on an ideal spectrum: jitters the true peaks,
/// adds a sampled noise floor, and collapses peaks that could not have been resolved apart.
/// </summary>
/// <remarks>
/// <para>
/// Noise must be injected <b>after</b> <see cref="Simulation.SimulatedScanReducer"/>, not before.
/// The reducer's default floor is a fraction of the brightest peak in the entire run, which for a
/// realistic simulation lands far above the noise amplitude, so noise added first would be
/// thresholded away in its entirety.
/// </para>
/// <para>
/// Jitter is applied only to peaks tagged <see cref="PeakOrigin.Signal"/>. Noise peaks are drawn
/// from a distribution measured on centroids a real instrument already reported, so their scatter is
/// baked in and jittering them again would double-count it.
/// </para>
/// <para>
/// Scans are generated independently and in parallel, each from a stream derived by mixing the seed
/// with the scan index. Output therefore does not depend on scheduling, and re-running with the same
/// seed reproduces the file exactly.
/// </para>
/// </remarks>
public sealed class NoiseInjector
{
    private readonly NoiseFloorModel? _model;
    private readonly IReadOnlyList<NoiseFloorModel>? _scanModels;
    private readonly IPeakWidthModel _widthModel;
    private readonly PeakJitterModel? _jitter;
    private readonly IIntensityFloor? _floor;
    private readonly IReadOnlyList<IIntensityFloor>? _scanFloors;
    private readonly double _mergeWithinSigmas;
    private readonly int _seed;

    /// <param name="jitter">
    /// Scan-to-scan measurement error on the true peaks, or null to leave them exact. Note that
    /// noise density and jitter are independent: <see cref="NoiseFloorModel.DensityScale"/> of 0
    /// gives jitter with no added noise floor.
    /// </param>
    /// <param name="mergeWithinSigmas">
    /// Two centroids closer than this many σ are reported by a real instrument as one peak. The
    /// default of 1.0 is a little under half the FWHM.
    /// </param>
    /// <param name="floor">
    /// The detection limit a jittered peak must still clear to be kept. Defaults to the noise
    /// model's own <see cref="NoiseRelativeFloor"/>, which is the floor reduction used.
    /// </param>
    public NoiseInjector(
        NoiseFloorModel model,
        IPeakWidthModel widthModel,
        int seed = 0,
        double mergeWithinSigmas = 1.0,
        PeakJitterModel? jitter = null,
        IIntensityFloor? floor = null)
    {
        _model = model ?? throw new ArgumentNullException(nameof(model));
        _widthModel = widthModel ?? throw new ArgumentNullException(nameof(widthModel));
        if (!(mergeWithinSigmas >= 0) || !double.IsFinite(mergeWithinSigmas))
            throw new ArgumentOutOfRangeException(nameof(mergeWithinSigmas), mergeWithinSigmas,
                "Merge width must be finite and non-negative.");

        _jitter = jitter is { IsActive: true } ? jitter : null;
        _floor = floor ?? new NoiseRelativeFloor(model);
        _mergeWithinSigmas = mergeWithinSigmas;
        _seed = seed;
    }

    /// <summary>
    /// Noise conditions that change from scan to scan: <paramref name="scanModels"/>[s] and
    /// <paramref name="scanFloors"/>[s] apply to scan s, and both must have one entry per scan
    /// passed to <see cref="Apply"/>. See <see cref="ScanNoiseConditions"/>.
    /// </summary>
    /// <param name="scanFloors">Per-scan detection limits, or null for each model's own <see cref="NoiseRelativeFloor"/>.</param>
    public NoiseInjector(
        IReadOnlyList<NoiseFloorModel> scanModels,
        IPeakWidthModel widthModel,
        int seed = 0,
        double mergeWithinSigmas = 1.0,
        PeakJitterModel? jitter = null,
        IReadOnlyList<IIntensityFloor>? scanFloors = null)
    {
        _scanModels = scanModels ?? throw new ArgumentNullException(nameof(scanModels));
        _widthModel = widthModel ?? throw new ArgumentNullException(nameof(widthModel));
        if (!(mergeWithinSigmas >= 0) || !double.IsFinite(mergeWithinSigmas))
            throw new ArgumentOutOfRangeException(nameof(mergeWithinSigmas), mergeWithinSigmas,
                "Merge width must be finite and non-negative.");
        if (scanFloors is not null && scanFloors.Count != scanModels.Count)
            throw new ArgumentException("There must be one floor per scan model.", nameof(scanFloors));

        _jitter = jitter is { IsActive: true } ? jitter : null;
        _scanFloors = scanFloors ?? scanModels.Select(m => (IIntensityFloor)new NoiseRelativeFloor(m)).ToArray();
        _mergeWithinSigmas = mergeWithinSigmas;
        _seed = seed;
    }

    private NoiseFloorModel ModelFor(int scan) => _scanModels is null ? _model! : _scanModels[scan];

    private IIntensityFloor FloorFor(int scan) => _scanFloors is null ? _floor! : _scanFloors[scan];

    public (MsDataScan[] Scans, NoiseInjectionSummary Summary) Apply(MsDataScan[] scans)
    {
        if (scans is null) throw new ArgumentNullException(nameof(scans));
        if (_scanModels is not null && _scanModels.Count != scans.Length)
            throw new ArgumentException(
                $"The injector was built for {_scanModels.Count} scans but was given {scans.Length}.", nameof(scans));

        var result = new MsDataScan[scans.Length];
        var signalCounts = new long[scans.Length];
        var noiseCounts = new long[scans.Length];
        var mergedCounts = new long[scans.Length];
        var droppedCounts = new long[scans.Length];

        Parallel.For(0, scans.Length, s =>
        {
            var peaks = new List<SimulatedPeak>();

            var spectrum = scans[s].MassSpectrum;
            if (_jitter is null)
            {
                for (int i = 0; i < spectrum.XArray.Length; i++)
                    peaks.Add(new SimulatedPeak(spectrum.XArray[i], spectrum.YArray[i], PeakOrigin.Signal));
            }
            else
            {
                droppedCounts[s] = AddJitteredSignal(
                    spectrum.XArray, spectrum.YArray, StreamForScan(_seed, s, JitterStream), peaks,
                    ModelFor(s), FloorFor(s));
            }

            signalCounts[s] = peaks.Count;

            int beforeNoise = peaks.Count;
            ModelFor(s).SampleScan(StreamForScan(_seed, s, NoiseStream), peaks);
            noiseCounts[s] = peaks.Count - beforeNoise;

            // Jitter can reorder signal peaks past one another, and SampleScan appends its own
            // ascending run, so the list is not sorted overall. A full sort is the simplest correct
            // fix and is not the bottleneck next to generating the peaks.
            peaks.Sort(static (a, b) => a.Mz.CompareTo(b.Mz));

            var (mz, intensities, merged) = Collapse(peaks);
            mergedCounts[s] = merged;
            result[s] = CloneWithSpectrum(scans[s], mz, intensities);
        });

        long signal = 0, noise = 0, mergedTotal = 0, dropped = 0;
        for (int s = 0; s < scans.Length; s++)
        {
            signal += signalCounts[s];
            noise += noiseCounts[s];
            mergedTotal += mergedCounts[s];
            dropped += droppedCounts[s];
        }

        return (result, new NoiseInjectionSummary(signal, noise, mergedTotal, dropped));
    }

    /// <summary>
    /// Adds one scan's signal peaks with measurement error applied, returning how many fell below
    /// the detection floor and were dropped.
    /// </summary>
    /// <remarks>
    /// The common-mode m/z offset and intensity factor are drawn once for the scan and applied to
    /// every peak in it, which is what makes them common-mode; the per-peak terms are drawn
    /// individually and scale with each peak's own S/N.
    /// </remarks>
    private long AddJitteredSignal(
        double[] mz, double[] intensities, Random rng, List<SimulatedPeak> into,
        NoiseFloorModel model, IIntensityFloor floor)
    {
        var jitter = _jitter!;
        double commonModePpm = jitter.CommonModeMzPpm > 0
            ? jitter.CommonModeMzPpm * NoiseFloorModel.SampleStandardNormal(rng)
            : 0;
        double commonModeIntensity = PeakJitterModel.IntensityFactor(rng, jitter.CommonModeLogIntensitySigma);

        long dropped = 0;
        for (int i = 0; i < mz.Length; i++)
        {
            double noiseLevel = model.NoiseLevelAt(mz[i]);
            double signalToNoise = noiseLevel > 0 ? intensities[i] / noiseLevel : double.PositiveInfinity;

            double ppm = commonModePpm + jitter.MzPpmSigma(signalToNoise) * NoiseFloorModel.SampleStandardNormal(rng);
            double jitteredMz = mz[i] * (1.0 + ppm * 1e-6);

            double jitteredIntensity = intensities[i]
                * commonModeIntensity
                * PeakJitterModel.IntensityFactor(rng, jitter.LogIntensitySigma(signalToNoise));

            if (jitter.DropPeaksFallingBelowFloor && jitteredIntensity < floor.At(jitteredMz))
            {
                dropped++;
                continue;
            }

            into.Add(new SimulatedPeak(jitteredMz, jitteredIntensity, PeakOrigin.Signal));
        }

        return dropped;
    }

    /// <summary>
    /// Collapses groups of peaks the instrument could not have reported separately into one centroid
    /// each, at their intensity-weighted m/z with their summed intensity.
    /// </summary>
    /// <remarks>
    /// <para>
    /// A group is bounded by its own first member rather than grown transitively from the last one
    /// accepted. Single-linkage chaining would let a dense run of peaks each just inside the
    /// tolerance collapse into a single centroid spanning many σ, which at the measured noise
    /// density — mean spacing of about 2σ around m/z 700 — swallows most of a spectrum.
    /// </para>
    /// <para>
    /// A group made up purely of noise is left alone. The density curve those peaks were drawn from
    /// was measured from the centroids a real instrument <i>reported</i>, so it already accounts for
    /// whatever the instrument's own peak detection merged; merging again here would apply that
    /// reduction twice. Groups containing signal are collapsed, because there the coincidence is
    /// genuine — a noise peak landing on an isotopologue does hide it.
    /// </para>
    /// </remarks>
    private (double[] Mz, double[] Intensities, long Merged) Collapse(List<SimulatedPeak> sorted)
    {
        var mz = new List<double>(sorted.Count);
        var intensities = new List<double>(sorted.Count);
        long merged = 0;

        int i = 0;
        while (i < sorted.Count)
        {
            double groupStart = sorted[i].Mz;
            double tolerance = _mergeWithinSigmas * _widthModel.SigmaAt(groupStart);

            int j = i + 1;
            while (j < sorted.Count && sorted[j].Mz - groupStart <= tolerance)
                j++;

            int count = j - i;
            if (count == 1)
            {
                mz.Add(sorted[i].Mz);
                intensities.Add(sorted[i].Intensity);
                i = j;
                continue;
            }

            bool hasSignal = false;
            for (int k = i; k < j && !hasSignal; k++)
                hasSignal = sorted[k].Origin == PeakOrigin.Signal;

            if (!hasSignal)
            {
                for (int k = i; k < j; k++)
                {
                    mz.Add(sorted[k].Mz);
                    intensities.Add(sorted[k].Intensity);
                }

                i = j;
                continue;
            }

            double weightedMz = 0, summed = 0;
            for (int k = i; k < j; k++)
            {
                weightedMz += sorted[k].Mz * sorted[k].Intensity;
                summed += sorted[k].Intensity;
            }

            // Falls back to the first peak's position if the whole group is zero-intensity, which
            // keeps the m/z finite rather than 0/0.
            mz.Add(summed > 0 ? weightedMz / summed : sorted[i].Mz);
            intensities.Add(summed);
            merged += count - 1;
            i = j;
        }

        return (mz.ToArray(), intensities.ToArray(), merged);
    }

    /// <summary>Stream identifiers, so jitter and noise draw from independent sequences.</summary>
    /// <remarks>
    /// Sharing one stream would couple the two: enabling jitter consumes draws before the noise
    /// sampler runs, so the same seed would produce a different noise floor. That makes an
    /// otherwise-clean A/B comparison between jittered and unjittered output impossible.
    /// </remarks>
    internal const int NoiseStream = 0;

    internal const int JitterStream = 1;

    /// <summary>
    /// A stream for one scan, derived from the seed by SplitMix64 so that adjacent scan indices do
    /// not produce correlated sequences the way adjacent <see cref="Random"/> seeds can.
    /// </summary>
    internal static Random StreamForScan(int seed, int scanIndex, int stream = NoiseStream)
    {
        unchecked
        {
            ulong z = (ulong)(uint)seed * 0x9E3779B97F4A7C15UL + (ulong)(uint)scanIndex;
            z += (ulong)(uint)stream * 0xD1B54A32D192ED03UL;
            z += 0x9E3779B97F4A7C15UL;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
            z ^= z >> 31;
            return new Random((int)(z & 0x7FFFFFFF));
        }
    }

    private static MsDataScan CloneWithSpectrum(MsDataScan source, double[] mz, double[] intensities)
    {
        double tic = 0;
        for (int i = 0; i < intensities.Length; i++)
            tic += intensities[i];

        var scanWindowRange = mz.Length > 0
            ? new MzRange(mz[0], mz[^1])
            : source.ScanWindowRange;

        return new MsDataScan(
            massSpectrum: new MzSpectrum(mz, intensities, false),
            oneBasedScanNumber: source.OneBasedScanNumber,
            msnOrder: source.MsnOrder,
            isCentroid: true,
            polarity: source.Polarity,
            retentionTime: source.RetentionTime,
            scanWindowRange: scanWindowRange,
            scanFilter: source.ScanFilter,
            mzAnalyzer: source.MzAnalyzer,
            totalIonCurrent: tic,
            injectionTime: source.InjectionTime,
            noiseData: null,
            nativeId: source.NativeId);
    }
}
