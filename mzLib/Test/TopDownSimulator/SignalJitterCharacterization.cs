using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;
using TopDownSimulator.Model;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoFisher.CommonCore.RawFileReader;

namespace Test.TopDownSimulator;

/// <summary>
/// Measures how much a real peak's m/z and intensity wobble from scan to scan, which is the
/// difference between the simulator's perfectly smooth analytical XICs and real data.
/// </summary>
/// <remarks>
/// <para>
/// Both estimators are built to cancel the elution profile rather than to model it:
/// </para>
/// <list type="bullet">
/// <item><b>m/z jitter</b> is each observation's deviation from its own isotopologue's median m/z
/// across scans, split into the part shared by every species in a scan (calibration drift) and the
/// part left over (per-peak centroiding error).</item>
/// <item><b>Intensity jitter</b> is the scan-to-scan scatter of the <i>ratio</i> between two
/// adjacent isotopologues of the same species. That ratio is fixed by chemistry, so the elution
/// profile divides out exactly and what remains is measurement noise.</item>
/// </list>
/// <para>Diagnostic only: prints, asserts nothing, hence <c>[Explicit]</c>.</para>
/// </remarks>
[TestFixture]
public class SignalJitterCharacterization
{
    private const string RawPath = @"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw";
    private const string ResultPath =
        @"D:\JurkatTopdown\Frac7_GPTMD_Search\Task2-TopDownSearch\Individual File Results\02-18-20_jurkat_td_rep2_fract7_Proteoforms.psmtsv";

    private const double PpmTolerance = 15.0;
    private const double RtHalfWidth = 0.35;
    private const int MaxProteoforms = 400;

    /// <summary>One observation of one isotopologue in one scan.</summary>
    private readonly record struct Observation(
        int SeriesId, int ScanNumber, int IsotopologueIndex,
        double TheoreticalMz, double ObservedMz, double Intensity, double Noise)
    {
        public double PpmDeviation => (ObservedMz - TheoreticalMz) / TheoreticalMz * 1e6;
        public double SignalToNoise => Noise > 0 ? Intensity / Noise : double.NaN;
    }

    [Test]
    [Explicit("Reads the local Jurkat raw plus MetaMorpheus results and prints jitter statistics")]
    public static void CharacterizeSignalJitter()
    {
        if (!File.Exists(RawPath) || !File.Exists(ResultPath))
        {
            Console.WriteLine($"SKIPPED - raw exists {File.Exists(RawPath)}, results exist {File.Exists(ResultPath)}");
            return;
        }

        var targets = LoadTargets();
        Console.WriteLine($"Tracking {targets.Count} proteoforms");

        var observations = CollectObservations(targets);
        Console.WriteLine($"Collected {observations.Count} isotopologue observations " +
                          $"across {observations.Select(o => o.ScanNumber).Distinct().Count()} scans");

        ReportCommonModeDrift(observations);
        ReportMzJitter(observations);
        ReportIntensityJitter(observations);
    }

    private sealed record Target(double MonoisotopicMass, int Charge, double RetentionTime);

    private static List<Target> LoadTargets()
    {
        var file = new PsmFromTsvFile(ResultPath, new SpectrumMatchParsingParameters
        {
            ParseMatchedFragmentIons = false,
        });
        file.LoadResults();

        return file.Results
            .Where(p => p is not null && p.MonoisotopicMass > 0 && p.RetentionTime >= 0)
            .Where(p => !double.IsNaN(p.QValue) && p.QValue <= 0.01)
            .Where(p => p.PrecursorCharge > 0)
            .GroupBy(p => (p.FullSequence, p.PrecursorCharge))
            .Select(g => g.OrderByDescending(p => p.Score).First())
            .OrderByDescending(p => p.Score)
            .Take(MaxProteoforms)
            .Select(p => new Target(p.MonoisotopicMass, p.PrecursorCharge, p.RetentionTime))
            .ToList();
    }

    /// <summary>
    /// Walks the scans once, matching every tracked isotopologue against the centroids in that scan.
    /// </summary>
    private static List<Observation> CollectObservations(List<Target> targets)
    {
        // Theoretical isotopologue m/z per target, restricted to the strongest few so that the
        // statistics are not dominated by envelope tails sitting in the noise.
        var series = new List<(int SeriesId, double Mz, int IsoIndex, double RtStart, double RtEnd)>();
        for (int t = 0; t < targets.Count; t++)
        {
            var target = targets[t];
            var kernel = new IsotopeEnvelopeKernel(target.MonoisotopicMass);
            var centroids = kernel.CentroidMzs(target.Charge);
            if (centroids.Length < 3)
                continue;

            // The middle of the envelope, where the isotopologues are brightest.
            int start = Math.Max(0, centroids.Length / 2 - 3);
            int end = Math.Min(centroids.Length, start + 6);
            for (int i = start; i < end; i++)
            {
                series.Add((t, centroids[i], i,
                    target.RetentionTime - RtHalfWidth, target.RetentionTime + RtHalfWidth));
            }
        }

        var observations = new List<Observation>();
        var rawFile = RawFileReaderAdapter.FileFactory(RawPath);
        try
        {
            rawFile.SelectInstrument(Device.MS, 1);
            int first = rawFile.RunHeaderEx.FirstSpectrum;
            int last = rawFile.RunHeaderEx.LastSpectrum;

            for (int n = first; n <= last; n++)
            {
                var filter = rawFile.GetFilterForScanNumber(n);
                if (filter is null || filter.MSOrder != ThermoFisher.CommonCore.Data.FilterEnums.MSOrderType.Ms)
                    continue;

                double rt = rawFile.RetentionTimeFromScanNumber(n);
                var active = series.Where(s => rt >= s.RtStart && rt <= s.RtEnd).ToArray();
                if (active.Length == 0)
                    continue;

                var stream = rawFile.GetCentroidStream(n, false);
                if (stream?.Masses is null || stream.Intensities is null)
                    continue;

                var mzs = stream.Masses;
                var intensities = stream.Intensities;
                var noises = stream.Noises;

                var order = Enumerable.Range(0, mzs.Length).ToArray();
                Array.Sort(order.Select(i => mzs[i]).ToArray(), order);
                var sortedMz = order.Select(i => mzs[i]).ToArray();

                foreach (var s in active)
                {
                    int match = NearestWithin(sortedMz, s.Mz, PpmTolerance);
                    if (match < 0)
                        continue;

                    int original = order[match];
                    observations.Add(new Observation(
                        SeriesId: s.SeriesId,
                        ScanNumber: n,
                        IsotopologueIndex: s.IsoIndex,
                        TheoreticalMz: s.Mz,
                        ObservedMz: sortedMz[match],
                        Intensity: intensities[original],
                        Noise: noises is { Length: > 0 } ? noises[original] : double.NaN));
                }
            }
        }
        finally
        {
            rawFile.Dispose();
        }

        return observations;
    }

    /// <summary>
    /// The part of the m/z error shared by every peak in a scan. If this dominates, jitter is
    /// calibration drift and correlates across a spectrum; if it does not, jitter is independent
    /// per-peak centroiding error. The two stress a feature finder very differently.
    /// </summary>
    private static void ReportCommonModeDrift(List<Observation> observations)
    {
        var perScan = observations
            .GroupBy(o => o.ScanNumber)
            .Where(g => g.Count() >= 5)
            .Select(g => Median(g.Select(o => o.PpmDeviation).ToArray()))
            .ToArray();

        if (perScan.Length == 0)
        {
            Console.WriteLine("common-mode: too few observations");
            return;
        }

        Console.WriteLine();
        Console.WriteLine("---- per-scan common-mode m/z offset (ppm) ----");
        Console.WriteLine($"scans: {perScan.Length}, mean {perScan.Average():F3}, sd {StandardDeviation(perScan):F3}");
        Report("common-mode ppm", perScan);
    }

    /// <summary>
    /// Per-peak m/z scatter after the per-scan common mode is removed, against S/N. Centroiding
    /// error should fall as roughly 1/(S/N), so this is where that gets confirmed or refuted.
    /// </summary>
    private static void ReportMzJitter(List<Observation> observations)
    {
        var commonModeByScan = observations
            .GroupBy(o => o.ScanNumber)
            .Where(g => g.Count() >= 5)
            .ToDictionary(g => g.Key, g => Median(g.Select(o => o.PpmDeviation).ToArray()));

        // Deviation from the series' own median, so a systematically mis-assigned monoisotopic mass
        // shifts the whole series and cancels rather than inflating the scatter.
        var residuals = new List<(double Residual, double SignalToNoise, double Mz)>();
        foreach (var group in observations.GroupBy(o => (o.SeriesId, o.IsotopologueIndex)))
        {
            var members = group.Where(o => commonModeByScan.ContainsKey(o.ScanNumber)).ToArray();
            if (members.Length < 5)
                continue;

            var corrected = members.Select(o => o.PpmDeviation - commonModeByScan[o.ScanNumber]).ToArray();
            double median = Median(corrected);
            for (int i = 0; i < members.Length; i++)
                residuals.Add((corrected[i] - median, members[i].SignalToNoise, members[i].ObservedMz));
        }

        if (residuals.Count == 0)
        {
            Console.WriteLine("m/z jitter: too few observations");
            return;
        }

        Console.WriteLine();
        Console.WriteLine("---- per-peak m/z residual after common-mode removal (ppm) ----");
        Console.WriteLine($"n {residuals.Count}, sd {StandardDeviation(residuals.Select(r => r.Residual).ToArray()):F3} ppm");
        Report("residual ppm", residuals.Select(r => r.Residual).ToArray());

        Console.WriteLine();
        Console.WriteLine($"{"S/N bin",14} {"n",8} {"sd (ppm)",10} {"sd x S/N",10}");
        foreach (var (lo, hi) in new[] { (2.0, 4.0), (4.0, 8.0), (8.0, 16.0), (16.0, 32.0), (32.0, 64.0), (64.0, 1e9) })
        {
            var inBin = residuals.Where(r => r.SignalToNoise >= lo && r.SignalToNoise < hi)
                                 .Select(r => r.Residual).ToArray();
            if (inBin.Length < 50) continue;

            double sd = StandardDeviation(inBin);
            double midpoint = hi > 1e8 ? 128 : Math.Sqrt(lo * hi);
            Console.WriteLine($"{lo,6:F0}-{hi,6:F0} {inBin.Length,8} {sd,10:F3} {sd * midpoint,10:F2}");
        }
    }

    /// <summary>
    /// Scan-to-scan scatter of the ratio between adjacent isotopologues of the same species. The
    /// true ratio is fixed by chemistry, so the elution profile divides out and what is left is
    /// measurement noise. Reported against the fainter peak of the pair, since that is what limits
    /// the pair.
    /// </summary>
    private static void ReportIntensityJitter(List<Observation> observations)
    {
        var bySeriesScan = observations
            .GroupBy(o => (o.SeriesId, o.ScanNumber))
            .ToDictionary(g => g.Key, g => g.ToDictionary(o => o.IsotopologueIndex, o => o));

        // ln of the adjacent-isotopologue ratio, grouped by (series, isotopologue pair).
        var byPair = new Dictionary<(int, int), List<(double LogRatio, double MinIntensity, double MinSn)>>();
        foreach (var ((seriesId, _), byIso) in bySeriesScan)
        {
            foreach (var (iso, lower) in byIso)
            {
                if (!byIso.TryGetValue(iso + 1, out var upper))
                    continue;
                if (!(lower.Intensity > 0) || !(upper.Intensity > 0))
                    continue;

                var key = (seriesId, iso);
                if (!byPair.TryGetValue(key, out var list))
                    byPair[key] = list = new List<(double, double, double)>();

                list.Add((Math.Log(lower.Intensity / upper.Intensity),
                    Math.Min(lower.Intensity, upper.Intensity),
                    Math.Min(lower.SignalToNoise, upper.SignalToNoise)));
            }
        }

        var scatter = new List<(double Residual, double MinIntensity, double MinSn)>();
        foreach (var samples in byPair.Values)
        {
            if (samples.Count < 5) continue;
            double median = Median(samples.Select(s => s.LogRatio).ToArray());
            foreach (var s in samples)
                scatter.Add((s.LogRatio - median, s.MinIntensity, s.MinSn));
        }

        if (scatter.Count == 0)
        {
            Console.WriteLine("intensity jitter: too few observations");
            return;
        }

        Console.WriteLine();
        Console.WriteLine("---- scan-to-scan scatter of adjacent-isotopologue log ratio ----");
        Console.WriteLine($"n {scatter.Count}, sd(ln ratio) {StandardDeviation(scatter.Select(s => s.Residual).ToArray()):F4}");
        Console.WriteLine("Per-peak relative intensity SD is this divided by sqrt(2), if both peaks jitter alike.");

        Console.WriteLine();
        Console.WriteLine($"{"S/N bin",14} {"n",8} {"sd(ln r)",10} {"per-peak CV",12} {"CV x sqrt(S/N)",15}");
        foreach (var (lo, hi) in new[] { (2.0, 4.0), (4.0, 8.0), (8.0, 16.0), (16.0, 32.0), (32.0, 64.0), (64.0, 256.0), (256.0, 1e9) })
        {
            var inBin = scatter.Where(s => s.MinSn >= lo && s.MinSn < hi).Select(s => s.Residual).ToArray();
            if (inBin.Length < 50) continue;

            double sd = StandardDeviation(inBin);
            double cv = sd / Math.Sqrt(2);
            double midpoint = hi > 1e8 ? 512 : Math.Sqrt(lo * hi);
            Console.WriteLine($"{lo,6:F0}-{hi,6:F0} {inBin.Length,8} {sd,10:F4} {cv,12:F4} {cv * Math.Sqrt(midpoint),15:F3}");
        }
    }

    private static int NearestWithin(double[] sortedMz, double target, double ppmTolerance)
    {
        double tolerance = target * ppmTolerance * 1e-6;
        int idx = Array.BinarySearch(sortedMz, target);
        if (idx >= 0) return idx;
        idx = ~idx;

        int best = -1;
        double bestDistance = double.MaxValue;
        for (int i = idx - 1; i <= idx; i++)
        {
            if (i < 0 || i >= sortedMz.Length) continue;
            double distance = Math.Abs(sortedMz[i] - target);
            if (distance <= tolerance && distance < bestDistance)
            {
                bestDistance = distance;
                best = i;
            }
        }

        return best;
    }

    private static void Report(string name, double[] values)
    {
        var sorted = (double[])values.Clone();
        Array.Sort(sorted);
        Console.WriteLine($"{name,-16}: " + string.Join("  ",
            new[] { 0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99 }
                .Select(q => $"p{q * 100:F0}={Quantile(sorted, q):F3}")));
    }

    private static double Quantile(double[] sorted, double q)
    {
        if (sorted.Length == 0) return double.NaN;
        int idx = (int)Math.Round(q * (sorted.Length - 1));
        return sorted[Math.Clamp(idx, 0, sorted.Length - 1)];
    }

    private static double Median(double[] values)
    {
        if (values.Length == 0) return double.NaN;
        var sorted = (double[])values.Clone();
        Array.Sort(sorted);
        return Quantile(sorted, 0.5);
    }

    private static double StandardDeviation(double[] values)
    {
        if (values.Length < 2) return double.NaN;
        double mean = values.Average();
        double sum = 0;
        foreach (double v in values)
            sum += (v - mean) * (v - mean);
        return Math.Sqrt(sum / (values.Length - 1));
    }
}
