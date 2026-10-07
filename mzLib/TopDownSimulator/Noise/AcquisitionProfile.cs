#nullable enable
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;
using MassSpectrometry;

namespace TopDownSimulator.Noise;

/// <summary>
/// Injection time as automatic gain control sets it: long enough to collect the target charge at
/// the current ion flux, capped at the method's maximum.
/// </summary>
/// <remarks>
/// Reported intensities are per unit injection time, so a scan's TIC is its ion flux and
/// IT = min(max IT, <see cref="ChargeTarget"/> / (TIC + <see cref="BackgroundTic"/>)).
/// <see cref="BackgroundTic"/> is ion load the reported spectrum does not show: ions outside the
/// scan range, and what falls below the detection threshold.
/// </remarks>
public sealed record AutomaticGainControl(double MaxInjectionTime, double ChargeTarget, double BackgroundTic)
{
    public double InjectionTime(double tic) =>
        Math.Min(MaxInjectionTime, ChargeTarget / (Math.Max(0, tic) + BackgroundTic));

    /// <summary>
    /// Fits the model to a run's MS1 scans: the cap is the largest injection time seen, and for each
    /// candidate background the target is the median IT·(TIC + B) over the scans below the cap; the
    /// background with the smallest RMS error in log10 IT over all scans wins.
    /// </summary>
    public static AutomaticGainControl Fit(IReadOnlyList<MsDataScan> scans, IEnumerable<double>? backgroundCandidates = null)
    {
        var points = scans
            .Where(s => s.InjectionTime is > 0 && double.IsFinite(s.InjectionTime.Value))
            .Select(s => (It: s.InjectionTime!.Value, Tic: s.MassSpectrum.SumOfAllY))
            .ToArray();
        if (points.Length < 10) throw new ArgumentException("At least ten scans with an injection time are needed.", nameof(scans));

        double maxIt = points.Max(p => p.It);
        var free = points.Where(p => p.It < 0.98 * maxIt).ToArray();
        if (free.Length == 0) free = points;

        AutomaticGainControl? best = null;
        double bestError = double.PositiveInfinity;
        foreach (double background in backgroundCandidates ?? DefaultBackgrounds())
        {
            double target = Median(free.Select(p => p.It * (p.Tic + background)));
            var candidate = new AutomaticGainControl(maxIt, target, background);
            double error = points.Average(p => Math.Pow(Math.Log10(candidate.InjectionTime(p.Tic) / p.It), 2));
            if (error < bestError) (best, bestError) = (candidate, error);
        }

        return best!;
    }

    private static IEnumerable<double> DefaultBackgrounds()
    {
        yield return 0;
        for (double b = 1e4; b <= 1e10; b *= Math.Sqrt(10)) yield return b;
    }

    internal static double Median(IEnumerable<double> values)
    {
        var v = values.Where(double.IsFinite).OrderBy(x => x).ToArray();
        if (v.Length == 0) throw new InvalidOperationException("No finite values.");
        return v.Length % 2 == 1 ? v[v.Length / 2] : 0.5 * (v[v.Length / 2 - 1] + v[v.Length / 2]);
    }
}

/// <summary>
/// What a run's acquisition looks like as a function of retention time, learned from one run of a
/// method so that another run of it can be simulated without a template: when MS1 scans come, how
/// many noise peaks they report, and how injection time follows the ion load.
/// </summary>
/// <remarks>
/// Noise density follows neither injection time nor the bright ion load (R² 0.02–0.10 against the
/// TIC above S/N 10 on rep2 fract7, rep1 fract7 and rep2 fract6), and its curve over retention time is
/// nearly the same in all three, from ~400 peaks per scan before elution to ~16 000 in the wash. So it
/// is kept as a curve over retention time: a property of the gradient, not of the sample. The noise
/// amplitude does follow injection time, and injection time follows the simulated ion load through
/// <see cref="Agc"/>, whose target agreed across the three runs (1.17–1.29·10⁹).
/// </remarks>
/// <param name="BinWidth">Width of the retention-time bins, minutes.</param>
/// <param name="FirstBinStart">Retention time at which the first bin starts.</param>
/// <param name="Ms1Interval">Median time between MS1 scans in each bin, minutes.</param>
/// <param name="PeaksPerScanPerBin">Per RT bin, median noise peaks per scan in each 100-Th m/z bin.</param>
/// <param name="FirstScanTime">Retention time of the first MS1 scan.</param>
/// <param name="LastScanTime">Retention time of the last MS1 scan.</param>
public sealed record AcquisitionProfile(
    double BinWidth,
    double FirstBinStart,
    double[] Ms1Interval,
    double[][] PeaksPerScanPerBin,
    double FirstScanTime,
    double LastScanTime,
    AutomaticGainControl Agc,
    double NoiseTimesInjectionTime,
    string Source,
    double DimLoadRatio = 0)
{
    /// <summary>Learns the profile from one run's MS1 scans.</summary>
    public static AcquisitionProfile Learn(
        IReadOnlyList<MsDataScan> ms1Scans, string source, double binWidth = 1.0,
        double noiseTimesInjectionTime = NoiseFloorModel.JurkatNoiseTimesInjectionTime)
    {
        var scans = ms1Scans.Where(s => s.MsnOrder == 1).OrderBy(s => s.RetentionTime).ToArray();
        if (scans.Length < 10) throw new ArgumentException("At least ten MS1 scans are needed.", nameof(ms1Scans));

        double first = scans[0].RetentionTime, last = scans[^1].RetentionTime;
        double start = Math.Floor(first / binWidth) * binWidth;
        int bins = (int)Math.Ceiling((last - start) / binWidth + 1e-9);
        var noise = ScanNoiseConditions.FromSourceScans(scans, new NoiseFloorModel(), noiseTimesInjectionTime);

        var interval = new double[bins];
        var density = new double[bins][];
        for (int b = 0; b < bins; b++)
        {
            double lo = start + b * binWidth, hi = lo + binWidth;
            var idx = Enumerable.Range(0, scans.Length).Where(i => scans[i].RetentionTime >= lo && scans[i].RetentionTime < hi).ToArray();
            var gaps = idx.Where(i => i + 1 < scans.Length).Select(i => scans[i + 1].RetentionTime - scans[i].RetentionTime).ToArray();
            interval[b] = gaps.Length > 0 ? AutomaticGainControl.Median(gaps) : double.NaN;
            density[b] = idx.Length > 0
                ? Enumerable.Range(0, NoiseFloorModel.BinCount).Select(m => AutomaticGainControl.Median(idx.Select(i => noise[i].PeaksPerScanInBin(m)))).ToArray()
                : null!;
        }

        // Empty bins borrow their nearest neighbour's values.
        for (int b = 0; b < bins; b++)
        {
            if (double.IsNaN(interval[b]) || density[b] is null)
            {
                int near = Enumerable.Range(0, bins).Where(o => !double.IsNaN(interval[o]) && density[o] is not null)
                    .OrderBy(o => Math.Abs(o - b)).First();
                if (double.IsNaN(interval[b])) interval[b] = interval[near];
                density[b] ??= density[near];
            }
        }

        var agc = AutomaticGainControl.Fit(scans);

        // Ion load below S/N 10 per unit of load above it, over the scans where AGC is active. Rep2
        // fract7 2–4.6 in the busy region; a simulation renders the bright part and not the rest.
        var ratios = new List<double>();
        for (int i = 0; i < scans.Length; i++)
        {
            if (!(scans[i].InjectionTime < 0.98 * agc.MaxInjectionTime)) continue;
            double bright = 0, dim = 0;
            var mz = scans[i].MassSpectrum.XArray;
            var y = scans[i].MassSpectrum.YArray;
            for (int p = 0; p < mz.Length; p++)
            {
                if (y[p] >= 10 * noise[i].NoiseLevelAt(mz[p])) bright += y[p];
                else dim += y[p];
            }

            if (bright > 0) ratios.Add(dim / bright);
        }

        double dimLoadRatio = ratios.Count > 0 ? AutomaticGainControl.Median(ratios) : 0;

        return new AcquisitionProfile(binWidth, start, interval, density, first, last,
            agc, noiseTimesInjectionTime, source, dimLoadRatio);
    }

    private int BinOf(double rt) => Math.Clamp((int)Math.Floor((rt - FirstBinStart) / BinWidth), 0, Ms1Interval.Length - 1);

    /// <summary>MS1 scan times from the first to the last, stepping by each bin's interval.</summary>
    public double[] ScanTimes()
    {
        var times = new List<double>();
        for (double t = FirstScanTime; t <= LastScanTime; t += Ms1Interval[BinOf(t)])
            times.Add(t);
        return times.ToArray();
    }

    /// <summary>Noise peaks per scan in each 100-Th bin at <paramref name="rt"/>.</summary>
    public double[] DensityAt(double rt) => PeaksPerScanPerBin[BinOf(rt)];

    /// <summary>
    /// One noise model per scan of <paramref name="signalScans"/>: amplitude from the injection time
    /// AGC would set for that scan's ion load, density from the curve at its retention time. The ion
    /// load is the rendered signal's TIC times (1 + <see cref="DimLoadRatio"/>), for the load below
    /// S/N 10 that the signal model does not render.
    /// </summary>
    /// <param name="signalScans">The rendered signal, before noise.</param>
    /// <param name="template">Supplies dispersion, the S/N distribution and the density scale.</param>
    public NoiseFloorModel[] NoiseModels(IReadOnlyList<MsDataScan> signalScans, NoiseFloorModel? template = null)
    {
        template ??= new NoiseFloorModel();
        var models = new NoiseFloorModel[signalScans.Count];
        for (int s = 0; s < models.Length; s++)
        {
            double it = Agc.InjectionTime((1 + DimLoadRatio) * signalScans[s].MassSpectrum.SumOfAllY);
            models[s] = template.With(NoiseTimesInjectionTime / it, DensityAt(signalScans[s].RetentionTime));
        }

        return models;
    }

    public void Save(string path) =>
        File.WriteAllText(path, JsonSerializer.Serialize(this, new JsonSerializerOptions { WriteIndented = true }));

    public static AcquisitionProfile Load(string path) =>
        JsonSerializer.Deserialize<AcquisitionProfile>(File.ReadAllText(path))
        ?? throw new InvalidDataException($"{path} does not hold an acquisition profile.");
}
