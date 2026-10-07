using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using MassSpectrometry;
using NUnit.Framework;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.FilterEnums;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoFisher.CommonCore.RawFileReader;
using TopDownSimulator.Comparison;
using TopDownSimulator.Model;
using TopDownSimulator.Noise;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

/// <summary>
/// The iteration loop for making simulated scans look like real ones: re-simulates a handful of
/// chosen scans from saved model parameters and scores each against the real scan at the same
/// retention time.
/// </summary>
/// <remarks>
/// <para>
/// Re-simulating a few scans takes seconds rather than the minutes a full export needs, because
/// the models are read back from a ground-truth sidecar instead of being refitted. That isolates
/// everything downstream of fitting — the noise floor, jitter, the detection floor and merging —
/// which is where the noise model lives. Changes to fitting itself need a fresh export and a new
/// sidecar.
/// </para>
/// <para>
/// Each <see cref="Variant"/> is one way of producing the simulated scans. To try an idea, add a
/// variant and compare its rows against <c>current</c> and against <c>file</c>, the scans from the
/// mzML already on disk. <c>file</c> and <c>current</c> should agree closely, since they differ only
/// in the noise draw; if they do not, the sidecar no longer describes the file.
/// </para>
/// <para>
/// Real scans are read from the Thermo centroid stream directly so that the per-peak noise and
/// the injection time are available; mzLib's reader discards both.
/// </para>
/// </remarks>
[TestFixture]
public class ScanRealismComparison
{
    private const string RealPath = @"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw";
    /// <summary>
    /// Which export to compare against: <c>full</c> plus MZLIB_TOPDOWN_SIM_OUTPUT_TAG, the same tag
    /// <c>AnalysisExample.ExportRep2FullNoisySimulation</c> appends to the files it writes.
    /// </summary>
    private static string ExportLabel => "full" + (Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_OUTPUT_TAG")?.Trim() ?? "");

    private static string SimulatedPath => $@"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.{ExportLabel}.noisy.simulated.mzML";
    private static string ModelsPath => $@"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.{ExportLabel}.noisy.simulated.groundtruth.tsv";
    private const string OutputDirectory = @"D:\JurkatTopdown\scan-realism";

    /// <summary>The noise amplitude at m/z 650 the full noisy export was written with.</summary>
    private const double ExportNoiseLevel = 663.0;

    /// <summary>
    /// Scans covering the regimes the noise model has to handle, each with the regions worth
    /// scoring on their own. The whole calibrated range is always scored as well.
    /// </summary>
    private static readonly (int Scan, string Regime, (double Lo, double Hi)[] Regions)[] Targets =
    {
        (2534, "busy, IT 0.3 ms", new[] { (860.4, 864.4) }),
        (1840, "busy, IT 2.8 ms", Array.Empty<(double, double)>()),
        (3680, "post-elution, IT 48 ms", Array.Empty<(double, double)>()),
        (1200, "pre-elution, IT 50 ms", Array.Empty<(double, double)>()),
    };

    private static readonly (double Lo, double Hi) WholeScan = (600.0, 2000.0);

    public sealed record RealScan(
        int ScanNumber, int Ms1Ordinal, double RetentionTime, double InjectionTimeMs,
        double[] Mz, double[] Intensity, double[] Noise);

    public sealed record SimulationContext(
        ProteoformModel[] Models, int MinCharge, int MaxCharge, IPeakWidthModel WidthModel, RealScan[] Targets)
    {
        public double[] ScanTimes => Targets.Select(t => t.RetentionTime).ToArray();
    }

    /// <summary>One way of producing simulated scans, returned in the same order as the targets.</summary>
    public sealed record Variant(string Label, Func<SimulationContext, MsDataScan[]> Simulate);

    private static IEnumerable<Variant> Variants()
    {
        yield return new Variant("file", ReadFromSimulatedFile);

        yield return new Variant("current", ctx => new Simulator().PrepareMs1(
            ctx.Models, ctx.MinCharge, ctx.MaxCharge, ctx.WidthModel, ctx.ScanTimes,
            noise: new NoiseFloorModel(ExportNoiseLevel)).Ms1Scans);

        // Same detection floor and jitter, no noise peaks: what the signal model alone contributes.
        yield return new Variant("signal-only", ctx => new Simulator().PrepareMs1(
            ctx.Models, ctx.MinCharge, ctx.MaxCharge, ctx.WidthModel, ctx.ScanTimes,
            noise: new NoiseFloorModel(ExportNoiseLevel, densityScale: 0)).Ms1Scans);

        // Noise conditioned on each source scan (ScanNoiseConditions): amplitude from injection
        // time alone, then density as well.
        yield return new Variant("it-amplitude", ctx => Conditioned(ctx, conditionDensity: false));
        yield return new Variant("it+density", ctx => Conditioned(ctx, conditionDensity: true));
        yield return new Variant("signal-it", ctx => Conditioned(ctx, conditionDensity: false, densityScale: 0));

        // A stand-in for species-level deduplication until the sidecar comes from an export that
        // did it during fitting; models from such a sidecar are already deduplicated.
        yield return new Variant("it+density+dedup", ctx => Conditioned(
            ctx with { Models = DeduplicateByMass(ctx.Models) }, conditionDensity: true));
    }

    private static MsDataScan[] Conditioned(SimulationContext ctx, bool conditionDensity, double densityScale = 1.0)
    {
        var sources = ctx.Targets.Select(t => new MsDataScan(
            new MzSpectrum(t.Mz, t.Intensity, false), t.ScanNumber, 1, true, Polarity.Positive, t.RetentionTime,
            new MzLibUtil.MzRange(600, 2000), "real", MZAnalyzerType.Orbitrap, t.Intensity.Sum(),
            t.InjectionTimeMs, null, $"scan={t.ScanNumber}")).ToArray();
        var scanNoise = ScanNoiseConditions.FromSourceScans(
            sources, new NoiseFloorModel(ExportNoiseLevel, densityScale), conditionDensity: conditionDensity);

        return new Simulator().PrepareMs1(
            ctx.Models, ctx.MinCharge, ctx.MaxCharge, ctx.WidthModel, ctx.ScanTimes, scanNoise: scanNoise).Ms1Scans;
    }

    /// <summary>
    /// Keeps one model per species, where a species is a mass within 0.05 Da eluting within half a
    /// minute. Isobaric localization variants and isoforms with identical sequences were each fitted
    /// independently to the same peaks, so each one already claims the whole observed signal.
    /// </summary>
    private static ProteoformModel[] DeduplicateByMass(ProteoformModel[] models)
    {
        var kept = new List<ProteoformModel>();
        foreach (var m in models.OrderByDescending(m => m.Abundance))
        {
            if (kept.Any(k => Math.Abs(k.MonoisotopicMass - m.MonoisotopicMass) < 0.05
                              && Math.Abs(k.RtProfile.Mu - m.RtProfile.Mu) < 0.5))
                continue;
            kept.Add(m);
        }

        return kept.ToArray();
    }

    [Test]
    [Explicit("Re-simulates selected rep2 fract7 scans and scores them against the real raw")]
    public static void CompareSimulatedScansToReal()
    {
        var real = ReadRealTargets();
        var (models, minZ, maxZ, width) = LoadModels(ModelsPath);
        var context = new SimulationContext(models, minZ, maxZ, width, real);
        Console.WriteLine($"models {models.Length}, charges {minZ}-{maxZ}, width {width}");
        Directory.CreateDirectory(OutputDirectory);

        var results = new List<(Variant Variant, MsDataScan[] Scans)>();
        foreach (var variant in Variants())
        {
            var sw = System.Diagnostics.Stopwatch.StartNew();
            var scans = variant.Simulate(context);
            Assert.That(scans.Length, Is.EqualTo(real.Length), $"{variant.Label} returned the wrong number of scans");
            Console.WriteLine($"variant {variant.Label}: {sw.Elapsed.TotalSeconds:F1} s");
            results.Add((variant, scans));
        }

        var table = new StringBuilder();
        table.AppendLine(string.Join('\t', "scan", "rt", "it_ms", "region", "variant",
            "real_n", "sim_n", "n_ratio", "real_matched", "sim_matched", "cosine", "sqrt_cosine",
            "tic_ratio", "rel_ks", "real_below1pct", "sim_below1pct"));

        Console.WriteLine();
        Console.WriteLine($"{"scan",5} {"region",13} {"variant",-16} {"real n",7} {"sim n",7} {"n ratio",7} " +
                          $"{"r match",7} {"s match",7} {"cos",6} {"sqrtcos",7} {"TIC",8} {"KS",5} {"r<1%",5} {"s<1%",5}");

        for (int t = 0; t < real.Length; t++)
        {
            var target = Targets.Single(x => x.Scan == real[t].ScanNumber);
            Console.WriteLine($"---- scan {real[t].ScanNumber}, RT {real[t].RetentionTime:F2}, " +
                              $"IT {real[t].InjectionTimeMs:F3} ms, {target.Regime}");

            WriteSpectrum(Path.Combine(OutputDirectory, $"scan{real[t].ScanNumber}.real.tsv"),
                real[t].Mz, real[t].Intensity);

            foreach (var region in new[] { WholeScan }.Concat(target.Regions))
            {
                foreach (var (variant, scans) in results)
                {
                    var spectrum = scans[t].MassSpectrum;
                    var c = CentroidSpectrumComparison.Compare(
                        real[t].Mz, real[t].Intensity, spectrum.XArray, spectrum.YArray, region.Lo, region.Hi);

                    string label = $"{region.Lo:F1}-{region.Hi:F1}";
                    Console.WriteLine($"{real[t].ScanNumber,5} {label,13} {variant.Label,-16} {c.RealPeaks,7} {c.SimulatedPeaks,7} " +
                                      $"{c.PeakCountRatio,7:F2} {c.RealMatchedFraction,7:F2} {c.SimulatedMatchedFraction,7:F2} " +
                                      $"{c.Cosine,6:F3} {c.SqrtCosine,7:F3} {c.TicRatio,8:G3} {c.RelativeIntensityKs,5:F2} " +
                                      $"{c.RealFractionBelowOnePercent,5:F2} {c.SimulatedFractionBelowOnePercent,5:F2}");

                    table.AppendLine(string.Join('\t',
                        real[t].ScanNumber, F(real[t].RetentionTime), F(real[t].InjectionTimeMs), label, variant.Label,
                        c.RealPeaks, c.SimulatedPeaks, F(c.PeakCountRatio), F(c.RealMatchedFraction),
                        F(c.SimulatedMatchedFraction), F(c.Cosine), F(c.SqrtCosine), F(c.TicRatio),
                        F(c.RelativeIntensityKs), F(c.RealFractionBelowOnePercent), F(c.SimulatedFractionBelowOnePercent)));
                }
            }

            foreach (var (variant, scans) in results)
                WriteSpectrum(Path.Combine(OutputDirectory, $"scan{real[t].ScanNumber}.{variant.Label}.tsv"),
                    scans[t].MassSpectrum.XArray, scans[t].MassSpectrum.YArray);
        }

        string tablePath = Path.Combine(OutputDirectory, "comparison.tsv");
        File.WriteAllText(tablePath, table.ToString());
        Console.WriteLine();
        Console.WriteLine($"metrics: {tablePath}");
        Console.WriteLine($"spectra: {OutputDirectory}\\scan<N>.<variant>.tsv");
    }

    /// <summary>
    /// Per-scan injection time against the reported noise amplitude across the run. Their product
    /// is close to constant whenever AGC limits the fill, which is why one noise level cannot fit
    /// every scan.
    /// </summary>
    [Test]
    [Explicit("Prints per-scan injection time against the reported noise amplitude across the run")]
    public static void NoiseLevelVsInjectionTime()
    {
        using var raw = OpenRaw();
        Console.WriteLine($"{"scan",6} {"RT",7} {"IT ms",8} {"peaks",6} {"noise@650",10} {"noise@850",10} {"noise*IT@650",13}");
        for (int n = raw.RunHeaderEx.FirstSpectrum; n <= raw.RunHeaderEx.LastSpectrum; n++)
        {
            if (n % 40 > 1 || !IsMs1(raw, n)) continue;
            var stream = raw.GetCentroidStream(n, false);
            double it = InjectionTime(raw, n);
            double n650 = MedianNoiseNear(stream, 650), n850 = MedianNoiseNear(stream, 850);
            Console.WriteLine($"{n,6} {raw.RetentionTimeFromScanNumber(n),7:F2} {it,8:F3} {stream.Masses.Length,6} " +
                              $"{n650,10:F0} {n850,10:F0} {n650 * it,13:F0}");
        }
    }

    /// <summary>
    /// Which models put signal on one peak of one scan, and how much, next to what the real scan and
    /// mzLib's own reader of it report there. Use it when a region's TIC ratio is far from 1.
    /// </summary>
    [Test]
    [Explicit("Lists the models contributing to one simulated peak")]
    public static void AttributeSimulatedPeak()
    {
        const int scanNumber = 2534;
        const double targetMz = 861.9786;
        const double ppm = 10;

        var real = ReadRealTargets().Single(t => t.ScanNumber == scanNumber);
        double tol = targetMz * ppm * 1e-6;
        for (int i = 0; i < real.Mz.Length; i++)
            if (Math.Abs(real.Mz[i] - targetMz) < tol)
                Console.WriteLine($"real centroid stream: {real.Mz[i]:F4} {real.Intensity[i]:E3}");

        var reader = Readers.MsDataFileReader.GetDataFile(RealPath);
        reader.InitiateDynamicConnection();
        var viaMzLib = reader.GetOneBasedScanFromDynamicConnection(scanNumber);
        reader.CloseDynamicConnection();
        for (int i = 0; i < viaMzLib.MassSpectrum.XArray.Length; i++)
            if (Math.Abs(viaMzLib.MassSpectrum.XArray[i] - targetMz) < tol)
                Console.WriteLine($"mzLib reader:         {viaMzLib.MassSpectrum.XArray[i]:F4} {viaMzLib.MassSpectrum.YArray[i]:E3}");

        var (models, minZ, maxZ, width) = LoadModels(ModelsPath);
        var simulator = new Simulator();
        var contributions = new List<(ProteoformModel Model, double Intensity)>();
        foreach (var model in models.Where(m => Math.Abs(m.RtProfile.Mu - real.RetentionTime) < 2))
        {
            var scan = simulator.SimulateCentroided(new[] { model }, minZ, maxZ, width, new[] { real.RetentionTime }).Scans[0];
            double sum = 0;
            for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
                if (Math.Abs(scan.MassSpectrum.XArray[i] - targetMz) < tol)
                    sum += scan.MassSpectrum.YArray[i];
            if (sum > 0) contributions.Add((model, sum));
        }

        Console.WriteLine($"{contributions.Count} models contribute, summing to {contributions.Sum(c => c.Intensity):E3}");
        foreach (var (m, intensity) in contributions.OrderByDescending(c => c.Intensity).Take(25))
            Console.WriteLine($"  {intensity,10:E2}  M {m.MonoisotopicMass,10:F2}  A {m.Abundance,9:E2}  " +
                              $"rt {m.RtProfile.Mu:F2}/{m.RtProfile.Sigma:F2}  z {((GaussianChargeDistribution)m.ChargeDistribution).MuZ:F1}  " +
                              $"{m.Identifier![..Math.Min(60, m.Identifier.Length)]}");
    }

    /// <summary>
    /// How many noise-like peaks each scan carries, and where, against injection time and TIC. A
    /// noise-like peak has no charge assigned and S/N under 10, which excludes nearly all analyte.
    /// </summary>
    [Test]
    [Explicit("Prints the noise-peak density of real scans across the run")]
    public static void NoiseDensityAcrossRun()
    {
        double[] edges = { 600, 700, 800, 900, 1000, 1200, 1400, 2000 };
        using var raw = OpenRaw();
        Console.WriteLine($"{"scan",6} {"RT",6} {"IT",7} {"TIC",9} {"peaks",6} {"noise",6} " +
                          string.Join(" ", edges.Zip(edges.Skip(1), (a, b) => $"{a:F0}-{b:F0}".PadLeft(9))));
        for (int n = raw.RunHeaderEx.FirstSpectrum; n <= raw.RunHeaderEx.LastSpectrum; n++)
        {
            if (n % 60 > 1 || !IsMs1(raw, n)) continue;
            var stream = raw.GetCentroidStream(n, false);
            var noise = Enumerable.Range(0, stream.Masses.Length)
                .Where(i => stream.Charges[i] == 0 && stream.Intensities[i] < 10 * stream.Noises[i])
                .Select(i => stream.Masses[i]).ToArray();
            var bins = edges.Zip(edges.Skip(1), (a, b) => noise.Count(m => m >= a && m < b));
            double tic = raw.GetScanStatsForScanNumber(n).TIC;
            Console.WriteLine($"{n,6} {raw.RetentionTimeFromScanNumber(n),6:F1} {InjectionTime(raw, n),7:F2} {tic,9:E1} " +
                              $"{stream.Masses.Length,6} {noise.Length,6} " + string.Join(" ", bins.Select(b => $"{b,9}")));
        }
    }

    private static RealScan[] ReadRealTargets()
    {
        using var raw = OpenRaw();
        var wanted = Targets.Select(t => t.Scan).ToHashSet();
        var result = new List<RealScan>();
        int ordinal = 0;

        for (int n = raw.RunHeaderEx.FirstSpectrum; n <= raw.RunHeaderEx.LastSpectrum && result.Count < wanted.Count; n++)
        {
            if (!IsMs1(raw, n)) continue;
            ordinal++;
            if (!wanted.Contains(n)) continue;

            var stream = raw.GetCentroidStream(n, false);
            var order = Enumerable.Range(0, stream.Masses.Length).OrderBy(i => stream.Masses[i]).ToArray();
            result.Add(new RealScan(
                n, ordinal, raw.RetentionTimeFromScanNumber(n), InjectionTime(raw, n),
                order.Select(i => stream.Masses[i]).ToArray(),
                order.Select(i => stream.Intensities[i]).ToArray(),
                order.Select(i => stream.Noises[i]).ToArray()));
        }

        Assert.That(result.Count, Is.EqualTo(wanted.Count), "Some target scans are not MS1 scans in the raw file.");
        return result.ToArray();
    }

    /// <summary>The written mzML holds only MS1 scans, so its scan number is the real scan's MS1 ordinal.</summary>
    private static MsDataScan[] ReadFromSimulatedFile(SimulationContext ctx)
    {
        var file = Readers.MsDataFileReader.GetDataFile(SimulatedPath);
        file.InitiateDynamicConnection();
        try
        {
            return ctx.Targets.Select(t =>
            {
                var scan = file.GetOneBasedScanFromDynamicConnection(t.Ms1Ordinal);
                Assert.That(scan.RetentionTime, Is.EqualTo(t.RetentionTime).Within(1e-3),
                    $"Simulated scan {t.Ms1Ordinal} is not at the real scan's retention time.");
                return scan;
            }).ToArray();
        }
        finally
        {
            file.CloseDynamicConnection();
        }
    }

    /// <summary>Rebuilds the simulated models from a sidecar written by <see cref="Simulator.WriteGroundTruth"/>.</summary>
    private static (ProteoformModel[] Models, int MinCharge, int MaxCharge, IPeakWidthModel Width) LoadModels(string path)
    {
        var lines = File.ReadAllLines(path);
        var header = lines[0].Split('\t');
        int Col(string name) => Array.IndexOf(header, name);
        double D(string[] f, string name) => double.Parse(f[Col(name)], CultureInfo.InvariantCulture);

        var models = new List<ProteoformModel>();
        int minZ = 0, maxZ = 0;
        string widthText = "";
        foreach (string line in lines.Skip(1).Where(l => l.Length > 0))
        {
            var f = line.Split('\t');
            models.Add(new ProteoformModel(
                D(f, "MonoisotopicMass"), D(f, "Abundance"),
                new EmgProfile(D(f, "RtMu"), D(f, "RtSigma"), D(f, "RtTau")),
                new GaussianChargeDistribution(D(f, "ChargeMu"), D(f, "ChargeSigma")),
                f[Col("Identifier")]));
            minZ = int.Parse(f[Col("MinCharge")], CultureInfo.InvariantCulture);
            maxZ = int.Parse(f[Col("MaxCharge")], CultureInfo.InvariantCulture);
            widthText = f[Col("PeakWidthModel")];
        }

        var constant = Regex.Match(widthText, @"^Constant\(sigma=([0-9.eE+-]+)\)$");
        Assert.That(constant.Success, Is.True,
            $"Only a constant peak width can be rebuilt from the sidecar; it records '{widthText}'.");
        var width = new ConstantPeakWidth(double.Parse(constant.Groups[1].Value, CultureInfo.InvariantCulture));
        return (models.ToArray(), minZ, maxZ, width);
    }

    private static IRawDataPlus OpenRaw()
    {
        var raw = RawFileReaderAdapter.FileFactory(RealPath);
        Assert.That(raw.IsOpen && !raw.IsError, Is.True, $"Could not open {RealPath}");
        raw.SelectInstrument(Device.MS, 1);
        return raw;
    }

    private static bool IsMs1(IRawDataPlus raw, int scan) =>
        raw.GetFilterForScanNumber(scan)?.MSOrder == MSOrderType.Ms;

    private static double InjectionTime(IRawDataPlus raw, int scan)
    {
        var trailer = raw.GetTrailerExtraInformation(scan);
        for (int i = 0; i < trailer.Length; i++)
            if (trailer.Labels[i].StartsWith("Ion Injection Time", StringComparison.Ordinal))
                return double.Parse(trailer.Values[i], CultureInfo.InvariantCulture);
        return double.NaN;
    }

    private static double MedianNoiseNear(RealScan scan, double mz, double halfWidth = 50)
    {
        var values = Enumerable.Range(0, scan.Mz.Length)
            .Where(i => Math.Abs(scan.Mz[i] - mz) < halfWidth)
            .Select(i => scan.Noise[i])
            .OrderBy(x => x)
            .ToArray();
        return values.Length > 0 ? values[values.Length / 2] : double.NaN;
    }

    private static double MedianNoiseNear(CentroidStream stream, double mz, double halfWidth = 50)
    {
        var values = Enumerable.Range(0, stream.Masses.Length)
            .Where(i => Math.Abs(stream.Masses[i] - mz) < halfWidth)
            .Select(i => stream.Noises[i])
            .OrderBy(x => x)
            .ToArray();
        return values.Length > 0 ? values[values.Length / 2] : double.NaN;
    }

    private static void WriteSpectrum(string path, double[] mz, double[] intensity)
    {
        var sb = new StringBuilder("mz\tintensity\n");
        for (int i = 0; i < mz.Length; i++)
            sb.Append(F(mz[i])).Append('\t').Append(F(intensity[i])).Append('\n');
        File.WriteAllText(path, sb.ToString());
    }

    private static string F(double v) => v.ToString("R", CultureInfo.InvariantCulture);
}
