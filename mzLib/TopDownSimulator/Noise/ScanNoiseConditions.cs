#nullable enable
using System;
using System.Collections.Generic;
using MassSpectrometry;

namespace TopDownSimulator.Noise;

/// <summary>
/// Builds one <see cref="NoiseFloorModel"/> per scan from the real scans a simulation reuses the
/// retention times of, so the noise follows the run instead of being fixed at one calibration.
/// </summary>
/// <remarks>
/// <para>
/// Measured on rep2 fract7, a single model is wrong in both of the quantities it fixes:
/// </para>
/// <list type="bullet">
/// <item><b>Amplitude</b> goes as 1/injection time, so across the run it spans 180 to 49 000 at
/// m/z 650 (<see cref="NoiseFloorModel.JurkatNoiseTimesInjectionTime"/>). In a busy scan with a
/// 0.3 ms fill the floor sits about 100x above a quiet-window calibration, which is what hides faint
/// signal in the real data.</item>
/// <item><b>Density</b> runs from ~350 noise-like peaks per scan before elution to ~17 000 in the
/// wash, and its distribution over m/z changes too. It follows neither injection time nor TIC: the
/// first 25 min and the wash both run at the 50 ms maximum. Most of it is presumably faint,
/// unidentified analyte and chemical background, which no model of the identified proteoforms can
/// predict.</item>
/// </list>
/// <para>
/// So both are taken from the source scan. The amplitude comes from its injection time. The density
/// in each 100-Th bin is its count of peaks under <see cref="DefaultMaxSignalToNoise"/> times that
/// amplitude: the peaks a feature finder would have to reject as noise. Only m/z, intensity and
/// injection time are used, all of which mzLib's readers keep.
/// </para>
/// </remarks>
public static class ScanNoiseConditions
{
    /// <summary>
    /// Real peaks below this S/N count towards the density. Above it a peak is far more likely to be
    /// analyte the simulation models explicitly; the sampled noise itself almost never exceeds it
    /// (P ≈ 0.7 % under <see cref="SignalToNoiseDistribution.Orbitrap"/>).
    /// </summary>
    public const double DefaultMaxSignalToNoise = 10.0;

    /// <summary>
    /// One model per source scan, parallel to <paramref name="sourceScans"/>.
    /// </summary>
    /// <param name="template">Supplies the S/N distribution, dispersion and density scale.</param>
    /// <param name="noiseTimesInjectionTime">Amplitude at m/z 650 times injection time, intensity·ms.</param>
    /// <param name="conditionDensity">
    /// Take each scan's density from the source. False keeps the template's density table and
    /// conditions the amplitude only.
    /// </param>
    /// <remarks>
    /// A scan without a usable injection time keeps the template's amplitude, and its density is
    /// measured against that.
    /// </remarks>
    /// <param name="conditionAmplitude">
    /// Also take the amplitude in each 100-Th bin from the source scan: the median of its peaks below
    /// S/N 10 there, over <see cref="NoiseFloorModel.MedianIntensityOverLevel"/>. Injection time alone
    /// gets it wrong by a factor that depends on the regime and the m/z (rep2 fract7: 0.4–1.2× in the
    /// quiet scans before elution, 1.2–3.5× in the busy region, 1.8–3.3× in the wash, all at the same
    /// or similar injection time), and simulated noise came out up to 9× too bright at high m/z.
    /// Bins with fewer than <see cref="MinimumPeaksForLevel"/> such peaks keep the injection-time
    /// level, scaled by the median correction of the bins that had enough.
    /// </param>
    public static NoiseFloorModel[] FromSourceScans(
        IReadOnlyList<MsDataScan> sourceScans,
        NoiseFloorModel template,
        double noiseTimesInjectionTime = NoiseFloorModel.JurkatNoiseTimesInjectionTime,
        bool conditionDensity = true,
        double maxSignalToNoise = DefaultMaxSignalToNoise,
        bool conditionAmplitude = true)
    {
        if (sourceScans is null) throw new ArgumentNullException(nameof(sourceScans));
        if (template is null) throw new ArgumentNullException(nameof(template));
        if (!(noiseTimesInjectionTime > 0) || !double.IsFinite(noiseTimesInjectionTime))
            throw new ArgumentOutOfRangeException(nameof(noiseTimesInjectionTime), noiseTimesInjectionTime,
                "The amplitude constant must be finite and positive.");
        if (!(maxSignalToNoise > 0))
            throw new ArgumentOutOfRangeException(nameof(maxSignalToNoise), maxSignalToNoise, "Must be positive.");

        double medianOverLevel = template.MedianIntensityOverLevel;
        var models = new NoiseFloorModel[sourceScans.Count];
        for (int s = 0; s < models.Length; s++)
        {
            var scan = sourceScans[s];
            double? it = scan.InjectionTime;
            double level = it is > 0 && double.IsFinite(it.Value)
                ? noiseTimesInjectionTime / it.Value
                : template.NoiseLevelAtReferenceMz;

            var itModel = template.With(level);
            var model = itModel;
            double[]? levels = null;
            var mz = scan.MassSpectrum.XArray;
            var y = scan.MassSpectrum.YArray;
            if (conditionAmplitude)
            {
                // Twice: the second pass classifies peaks against the level the first one measured.
                for (int pass = 0; pass < 2; pass++)
                {
                    var measured = MeasureLevels(mz, y, model, maxSignalToNoise, medianOverLevel, itModel);
                    if (measured is null) break;
                    levels = measured;
                    model = template.WithLevels(levels);
                }
            }

            if (conditionDensity)
            {
                var density = MeasureDensity(mz, y, model, maxSignalToNoise);
                model = levels is null ? model.With(level, density) : model.WithLevels(levels, density);
            }

            models[s] = model;
        }

        return models;
    }

    /// <summary>Bins need this many peaks below the S/N cut for their median to set the level.</summary>
    public const int MinimumPeaksForLevel = 20;

    /// <summary>
    /// The noise level in each 100-Th bin implied by the median of the peaks below
    /// <paramref name="maxSignalToNoise"/> times <paramref name="classifier"/>'s level, or null when no
    /// bin has <see cref="MinimumPeaksForLevel"/> of them.
    /// </summary>
    /// <param name="injectionTimeModel">The amplitude bins without enough peaks fall back on, scaled.</param>
    public static double[]? MeasureLevels(
        double[] mz, double[] intensity, NoiseFloorModel classifier, double maxSignalToNoise,
        double medianOverLevel, NoiseFloorModel injectionTimeModel)
    {
        int bins = NoiseFloorModel.BinCount;
        double first = NoiseFloorModel.BinStart(0);
        double width = NoiseFloorModel.BinStart(1) - first;
        var perBin = new List<double>[bins];
        for (int b = 0; b < bins; b++) perBin[b] = new List<double>();
        for (int i = 0; i < mz.Length; i++)
        {
            int bin = (int)Math.Floor((mz[i] - first) / width);
            if (bin < 0 || bin >= bins) continue;
            if (intensity[i] < maxSignalToNoise * classifier.NoiseLevelAt(mz[i]))
                perBin[bin].Add(intensity[i]);
        }

        var levels = new double[bins];
        var corrections = new List<double>();
        var measured = new bool[bins];
        for (int b = 0; b < bins; b++)
        {
            if (perBin[b].Count < MinimumPeaksForLevel) continue;
            perBin[b].Sort();
            levels[b] = perBin[b][perBin[b].Count / 2] / medianOverLevel;
            measured[b] = levels[b] > 0;
            if (measured[b]) corrections.Add(levels[b] / injectionTimeModel.NoiseLevelAt(first + (b + 0.5) * width));
        }

        if (corrections.Count == 0) return null;
        corrections.Sort();
        double correction = corrections[corrections.Count / 2];
        for (int b = 0; b < bins; b++)
            if (!measured[b])
                levels[b] = correction * injectionTimeModel.NoiseLevelAt(first + (b + 0.5) * width);
        return levels;
    }

    /// <summary>
    /// Peaks per 100-Th bin with intensity under <paramref name="maxSignalToNoise"/> times the
    /// model's level at that m/z.
    /// </summary>
    public static double[] MeasureDensity(double[] mz, double[] intensity, NoiseFloorModel model, double maxSignalToNoise)
    {
        var counts = new double[NoiseFloorModel.BinCount];
        double first = NoiseFloorModel.BinStart(0);
        double width = NoiseFloorModel.BinStart(1) - first;

        for (int i = 0; i < mz.Length; i++)
        {
            int bin = (int)Math.Floor((mz[i] - first) / width);
            if (bin < 0 || bin >= counts.Length) continue;
            if (intensity[i] < maxSignalToNoise * model.NoiseLevelAt(mz[i]))
                counts[bin]++;
        }

        return counts;
    }
}
