using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using TopDownSimulator.Noise;

namespace Test.TopDownSimulator;

/// <summary>
/// What a simulation without a template run has to reproduce of the acquisition itself: how the
/// injection time follows the ion load (automatic gain control), how often MS1 scans come, and how
/// many noise peaks each scan reports.
/// </summary>
[TestFixture]
public class AgcCharacterization
{
    private static readonly string[] Raws =
    {
        @"D:\JurkatTopdown\Rep2_Raw\02-18-20_jurkat_td_rep2_fract7.raw",
        @"D:\JurkatTopdown\02-18-20_jurkat_td_rep1_fract7.raw",
        @"D:\JurkatTopdown\Rep2_Raw\02-18-20_jurkat_td_rep2_fract6.raw",
    };

    [Test]
    [Explicit("Measures injection time against TIC, MS1 cadence and noise density on rep2 and rep1 fract7")]
    public static void CharacterizeAcquisition()
    {
        foreach (string raw in Raws)
        {
            var file = MsDataFileReader.GetDataFile(raw);
            file.LoadAllStaticData();
            var ms1 = file.GetAllScansList().Where(s => s.MsnOrder == 1).OrderBy(s => s.RetentionTime).ToArray();
            Console.WriteLine($"=== {Path.GetFileName(raw)}: {ms1.Length} MS1 scans");

            var it = ms1.Select(s => s.InjectionTime ?? double.NaN).ToArray();
            var tic = ms1.Select(s => s.MassSpectrum.SumOfAllY).ToArray();
            double maxIt = it.Where(double.IsFinite).Max();
            int capped = it.Count(v => v >= 0.98 * maxIt);
            Console.WriteLine($"max IT {maxIt:F1} ms, scans at the cap {capped} ({capped / (double)ms1.Length:P0})");

            // Thermo reports intensities per unit injection time, so the ions collected are TIC x IT.
            var charge = ms1.Select((s, i) => tic[i] * it[i]).ToArray();
            var free = Enumerable.Range(0, ms1.Length).Where(i => it[i] < 0.98 * maxIt && tic[i] > 0).ToArray();
            Report("TIC x IT, uncapped scans", free.Select(i => charge[i]));
            Report("TIC x IT, capped scans", Enumerable.Range(0, ms1.Length).Where(i => it[i] >= 0.98 * maxIt).Select(i => charge[i]));
            var (a, b, r2) = Line(free.Select(i => Math.Log10(tic[i])).ToArray(), free.Select(i => Math.Log10(it[i])).ToArray());
            Console.WriteLine($"log10 IT = {a:F3} + {b:F3} log10 TIC over uncapped scans, R2 {r2:F3}");

            // IT = min(maxIT, c / (TIC + B)): B is ion load the reported spectrum does not show.
            double bestB = 0, bestC = 0, bestErr = double.PositiveInfinity;
            foreach (double bg in new[] { 0.0, 1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8, 3e8, 1e9 })
            {
                double c = Median(free.Select(i => it[i] * (tic[i] + bg)));
                double err = Enumerable.Range(0, ms1.Length).Where(i => tic[i] > 0 && double.IsFinite(it[i]))
                    .Average(i => Math.Pow(Math.Log10(Math.Min(maxIt, c / (tic[i] + bg))) - Math.Log10(it[i]), 2));
                Console.WriteLine($"  B {bg,8:G2}: c {c:G3}, rms log10 IT error {Math.Sqrt(err):F3}");
                if (err < bestErr) (bestErr, bestB, bestC) = (err, bg, c);
            }

            Console.WriteLine($"best: IT = min({maxIt:F0}, {bestC:G3} / (TIC + {bestB:G2})), rms log10 error {Math.Sqrt(bestErr):F3}");

            var template = new NoiseFloorModel(663.0);
            var noise = ScanNoiseConditions.FromSourceScans(ms1, template);

            // How the load below S/N 10 relates to the bright load above it, scan by scan.
            var brightTic = ms1.Select((s, i) => SignalPeaksTic(s, noise[i], true)).ToArray();
            var dimTic = ms1.Select((s, i) => SignalPeaksTic(s, noise[i], false)).ToArray();
            var density = noise.Select(m => Enumerable.Range(0, NoiseFloorModel.BinCount).Sum(b => m.PeaksPerScanInBin(b))).ToArray();
            var loaded = Enumerable.Range(0, ms1.Length).Where(i => brightTic[i] > 0 && dimTic[i] > 0).ToArray();
            var (da, db, dr2) = Line(loaded.Select(i => Math.Log10(brightTic[i])).ToArray(), loaded.Select(i => Math.Log10(dimTic[i])).ToArray());
            Console.WriteLine($"log10 TIC below S/N 10 = {da:F3} + {db:F3} log10 TIC above, R2 {dr2:F3} ({loaded.Length} scans)");
            var (na, nb, nr2) = Line(loaded.Select(i => Math.Log10(brightTic[i])).ToArray(), loaded.Select(i => Math.Log10(density[i])).ToArray());
            Console.WriteLine($"log10 noise peaks = {na:F3} + {nb:F3} log10 TIC above S/N 10, R2 {nr2:F3}");
            var (ra, rb, rr2) = Line(loaded.Select(i => Math.Log10(tic[i])).ToArray(), loaded.Select(i => Math.Log10(density[i])).ToArray());
            Console.WriteLine($"log10 noise peaks = {ra:F3} + {rb:F3} log10 TIC, R2 {rr2:F3}");

            // Cadence and density by retention time, in 5-minute bins.
            Console.WriteLine($"  {"RT",7} {"scans",6} {"MS1 dt (s)",10} {"IT ms",7} {"log10 TIC",9} {"peaks",7} {"noise peaks",11} {"log TIC>10",10} {"log TIC<10",10}");
            for (double t0 = 0; t0 < ms1[^1].RetentionTime; t0 += 5)
            {
                var idx = Enumerable.Range(0, ms1.Length).Where(i => ms1[i].RetentionTime >= t0 && ms1[i].RetentionTime < t0 + 5).ToArray();
                if (idx.Length < 2) continue;
                double dt = (ms1[idx[^1]].RetentionTime - ms1[idx[0]].RetentionTime) / (idx.Length - 1) * 60;
                Console.WriteLine($"  {t0,3:F0}-{t0 + 5,-3:F0} {idx.Length,6} {dt,10:F2} {Median(idx.Select(i => it[i])),7:F1} " +
                                  $"{Median(idx.Select(i => Math.Log10(1 + tic[i]))),9:F2} {Median(idx.Select(i => (double)ms1[i].MassSpectrum.XArray.Length)),7:F0} " +
                                  $"{Median(idx.Select(i => density[i])),11:F0} {Median(idx.Select(i => Math.Log10(1 + brightTic[i]))),10:F2} {Median(idx.Select(i => Math.Log10(1 + dimTic[i]))),10:F2}");
            }
        }
    }

    /// <summary>
    /// The shape of the noise intensity distribution on rep2 fract7, free of its amplitude: every
    /// peak below S/N 10 divided by the median of its scan's peaks in the same 100-Th band, pooled by
    /// regime. Also the measured median against the injection-time model's level, per band.
    /// </summary>
    [Test]
    [Explicit("Measures the scale-free noise intensity shape and the amplitude against 2e4/IT on rep2 fract7")]
    public static void CharacterizeNoiseIntensity()
    {
        var file = MsDataFileReader.GetDataFile(Raws[0]);
        file.LoadAllStaticData();
        var ms1 = file.GetAllScansList().Where(s => s.MsnOrder == 1).OrderBy(s => s.RetentionTime).ToArray();
        var regimes = new (string Name, double From, double To)[] { ("quiet 0-28", 0, 28), ("busy 35-45", 35, 45), ("wash 55-70", 55, 70) };
        double[] qs = { 0.1, 0.25, 0.75, 0.9, 0.99 };

        foreach (var (name, from, to) in regimes)
        {
            var pooled = new List<double>();
            var ratioByBand = Enumerable.Range(0, NoiseFloorModel.BinCount).Select(_ => new List<double>()).ToArray();
            foreach (var s in ms1.Where(s => s.RetentionTime >= from && s.RetentionTime < to))
            {
                var model = new NoiseFloorModel(NoiseFloorModel.JurkatNoiseTimesInjectionTime / (s.InjectionTime ?? 50));
                var mz = s.MassSpectrum.XArray;
                var y = s.MassSpectrum.YArray;
                for (int b = 0; b < NoiseFloorModel.BinCount; b++)
                {
                    double lo = NoiseFloorModel.BinStart(b), hi = lo + 100;
                    var v = new List<double>();
                    for (int i = 0; i < mz.Length; i++)
                        if (mz[i] >= lo && mz[i] < hi && y[i] < 10 * model.NoiseLevelAt(mz[i])) v.Add(y[i]);
                    if (v.Count < 20) continue;
                    v.Sort();
                    double median = v[v.Count / 2];
                    pooled.AddRange(v.Select(x => x / median));
                    ratioByBand[b].Add(median / model.NoiseLevelAt(lo + 50));
                }
            }

            pooled.Sort();
            Console.WriteLine($"{name}: {pooled.Count} peaks; intensity / scan-band median at " +
                              string.Join(", ", qs.Select(q => $"p{q * 100:F0} {pooled[(int)(q * pooled.Count)]:F2}")));
            Console.WriteLine("  median intensity / (2e4/IT model level), by band from m/z 600: " +
                              string.Join(" ", ratioByBand.Select(r => r.Count == 0 ? "-" : Median(r).ToString("F2"))));
        }

        // The model's own shape for comparison, at the current defaults.
        var rng = new Random(1);
        var current = new NoiseFloorModel(1000.0);
        var sample = new List<SimulatedPeak>();
        for (int i = 0; i < 20; i++) current.SampleScan(rng, sample);
        var band = sample.Where(p => p.Mz is >= 700 and < 800).Select(p => p.Intensity).OrderBy(x => x).ToList();
        double m = band[band.Count / 2];
        Console.WriteLine($"model at defaults (m/z 700-800): " + string.Join(", ", qs.Select(q => $"p{q * 100:F0} {band[(int)(q * band.Count)] / m:F2}")) +
                          $"; median / level {m / current.NoiseLevelAt(750):F2}");

        // Candidate shapes: S/N = threshold + exponential excess of mean ratio x threshold, times a
        // lognormal level spread of sigma d. Target: the wash and busy shapes above.
        double[] target = { 0.63, 0.775, 1.37, 1.85, 2.94 };
        var best = (Ratio: 0.0, Dispersion: 0.0, Error: double.PositiveInfinity);
        foreach (double ratio in new[] { 0.5, 0.6, 0.7, 0.8, 0.9 })
        foreach (double d in new[] { 0.0, 0.1, 0.15, 0.2, 0.25, 0.3 })
        {
            var distribution = new SignalToNoiseDistribution(1.565, ratio * 1.565);
            var r = new Random(2);
            var v = Enumerable.Range(0, 200000).Select(_ => distribution.Sample(r) * Math.Exp(d * NormalSample(r))).OrderBy(x => x).ToArray();
            double med = v[v.Length / 2];
            var shape = qs.Select(q => v[(int)(q * v.Length)] / med).ToArray();
            double error = shape.Zip(target, (a, b) => Math.Pow(Math.Log(a / b), 2)).Sum();
            Console.WriteLine($"  ratio {ratio:F1} d {d:F2}: " + string.Join(" ", shape.Select(x => x.ToString("F2"))) + $"  median {med:F2}  err {error:F4}");
            if (error < best.Error) best = (ratio, d, error);
        }

        Console.WriteLine($"best: excess ratio {best.Ratio:F1}, dispersion {best.Dispersion:F2}");
    }

    private static double NormalSample(Random r) =>
        Math.Sqrt(-2 * Math.Log(1 - r.NextDouble())) * Math.Cos(2 * Math.PI * r.NextDouble());

    private static double SignalPeaksTic(MsDataScan scan, NoiseFloorModel noise, bool above)
    {
        double total = 0;
        var mz = scan.MassSpectrum.XArray;
        var y = scan.MassSpectrum.YArray;
        for (int i = 0; i < mz.Length; i++)
            if ((y[i] >= 10 * noise.NoiseLevelAt(mz[i])) == above)
                total += y[i];
        return total;
    }

    private static (double Intercept, double Slope, double R2) Line(double[] x, double[] y)
    {
        double mx = x.Average(), my = y.Average();
        double sxx = x.Sum(v => (v - mx) * (v - mx)), sxy = x.Zip(y, (u, v) => (u - mx) * (v - my)).Sum();
        double syy = y.Sum(v => (v - my) * (v - my));
        double b = sxy / sxx;
        return (my - b * mx, b, sxy * sxy / (sxx * syy));
    }

    private static double Median(IEnumerable<double> values)
    {
        var v = values.Where(double.IsFinite).OrderBy(x => x).ToArray();
        return v.Length == 0 ? double.NaN : v[v.Length / 2];
    }

    private static void Report(string label, IEnumerable<double> values)
    {
        var v = values.Where(double.IsFinite).OrderBy(x => x).ToArray();
        if (v.Length == 0) { Console.WriteLine($"  {label}: none"); return; }
        double Q(double q) => v[(int)Math.Min(v.Length - 1, Math.Floor(q * v.Length))];
        Console.WriteLine($"  {label,-28} n={v.Length,5}  p5 {Q(0.05):G3}  p25 {Q(0.25):G3}  p50 {Q(0.5):G3}  p75 {Q(0.75):G3}  p95 {Q(0.95):G3}");
    }
}
