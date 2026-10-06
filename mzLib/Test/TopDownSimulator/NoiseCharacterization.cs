using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using NUnit.Framework;
using ThermoFisher.CommonCore.Data.Business;
using ThermoFisher.CommonCore.Data.Interfaces;
using ThermoFisher.CommonCore.RawFileReader;

namespace Test.TopDownSimulator;

/// <summary>
/// Measures what the noise in a real top-down run actually looks like, so the simulator's noise
/// model can be built against numbers rather than intuition.
/// </summary>
/// <remarks>
/// Reads the Thermo centroid stream directly rather than going through
/// <see cref="Readers.ThermoRawFileReader"/>, because the per-peak <c>Noises</c>,
/// <c>Baselines</c> and <c>Resolutions</c> columns are what this analysis is about and mzLib
/// discards all three (<c>ThermoRawFileReader.cs</c>, <c>noiseData: null</c>). Going direct also
/// means only the scans inside the RT windows get read, instead of all ~9000.
/// <para>Everything here is diagnostic: it prints and asserts nothing, hence <c>[Explicit]</c>.</para>
/// </remarks>
[TestFixture]
public class NoiseCharacterization
{
    /// <summary>The late-gradient window the user flagged, where analyte has mostly stopped eluting.</summary>
    private const double NoiseRtStart = 50.0;
    private const double NoiseRtEnd = 55.0;

    /// <summary>A busy window, for contrast.</summary>
    private const double BusyRtStart = 31.0;
    private const double BusyRtEnd = 35.0;

    /// <summary>
    /// Candidate locations per logical file, tried in order. The top-level copies of rep1 fract7
    /// and rep2 fract5 are held open by another process on this machine, so the duplicates under
    /// Rep2_Raw are listed as fallbacks.
    /// </summary>
    private static readonly (string Label, string[] Candidates)[] Files =
    {
        ("rep2_fract5", new[]
        {
            @"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract5.raw",
            @"D:\JurkatTopdown\Rep2_Raw\02-18-20_jurkat_td_rep2_fract5.raw",
        }),
        ("rep2_fract7", new[]
        {
            @"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw",
            @"D:\JurkatTopdown\Rep2_Raw\02-18-20_jurkat_td_rep2_fract7.raw",
        }),
        ("rep1_fract6", new[]
        {
            @"D:\JurkatTopdown\02-18-20_jurkat_td_rep1_fract6.raw",
        }),
        ("rep1_fract7", new[]
        {
            @"D:\JurkatTopdown\02-18-20_jurkat_td_rep1_fract7.raw",
        }),
    };

    /// <summary>One centroid, with the columns mzLib throws away.</summary>
    private readonly record struct Centroid(
        double Mz, double Intensity, double Noise, double Baseline, double Resolution, int Charge)
    {
        public double SignalToNoise => Noise > 0 ? Intensity / Noise : double.NaN;
    }

    private sealed record Ms1Scan(int ScanNumber, double RetentionTime, Centroid[] Peaks);

    [Test]
    [Explicit("Reads local Jurkat .raw files and prints noise statistics")]
    public static void CharacterizeNoise()
    {
        foreach (var (label, candidates) in Files)
        {
            string? path = candidates.FirstOrDefault(CanOpen);
            Console.WriteLine();
            Console.WriteLine("================================================================");
            Console.WriteLine(label);
            Console.WriteLine("================================================================");

            if (path is null)
            {
                Console.WriteLine("SKIPPED - no candidate path is readable:");
                foreach (string c in candidates)
                    Console.WriteLine($"    {c} (exists={File.Exists(c)})");
                continue;
            }

            Console.WriteLine($"path: {path}");
            try
            {
                Characterize(label, path);
            }
            catch (Exception ex)
            {
                Console.WriteLine($"FAILED: {ex.GetType().Name}: {ex.Message}");
            }
        }
    }

    /// <summary>
    /// Measures how much the noise amplitude varies <i>at a fixed m/z</i>, which is the dispersion a
    /// model built from per-bin medians alone would miss. Reads the TSV that
    /// <see cref="CharacterizeNoise"/> dumps, so run that first.
    /// </summary>
    [Test]
    [Explicit("Analyses the peak dumps written by CharacterizeNoise")]
    public static void CharacterizeWithinMzNoiseDispersion()
    {
        // Narrow enough that the level curve is nearly flat across the slice, so what is left is
        // scan-to-scan and peak-to-peak variation rather than the m/z trend.
        const double sliceWidth = 5.0;

        foreach (string dumpPath in Directory.GetFiles(ScratchDir(), "*.peaks.tsv"))
        {
            Console.WriteLine();
            Console.WriteLine($"---- {Path.GetFileName(dumpPath)} ----");

            var noiseBySlice = new Dictionary<int, List<double>>();
            var snBySlice = new Dictionary<int, List<double>>();

            using (var reader = new StreamReader(dumpPath))
            {
                reader.ReadLine();
                string? line;
                while ((line = reader.ReadLine()) is not null)
                {
                    var f = line.Split('\t');
                    if (f.Length < 5) continue;
                    if (!double.TryParse(f[2], NumberStyles.Float, CultureInfo.InvariantCulture, out double mz)) continue;
                    if (!double.TryParse(f[3], NumberStyles.Float, CultureInfo.InvariantCulture, out double intensity)) continue;
                    if (!double.TryParse(f[4], NumberStyles.Float, CultureInfo.InvariantCulture, out double noise)) continue;
                    if (!(noise > 0)) continue;

                    int slice = (int)(mz / sliceWidth);
                    if (!noiseBySlice.TryGetValue(slice, out var noises))
                        noiseBySlice[slice] = noises = new List<double>();
                    noises.Add(noise);

                    if (!snBySlice.TryGetValue(slice, out var sns))
                        snBySlice[slice] = sns = new List<double>();
                    sns.Add(intensity / noise);
                }
            }

            // Pool each slice's values divided by that slice's own median, which removes the m/z
            // trend and leaves the shape of the within-m/z variation.
            var pooledNoise = new List<double>();
            var pooledSn = new List<double>();
            foreach (var (slice, noises) in noiseBySlice)
            {
                if (noises.Count < 200) continue;
                double median = Median(noises.ToArray());
                if (!(median > 0)) continue;
                foreach (double n in noises) pooledNoise.Add(n / median);
                foreach (double sn in snBySlice[slice]) pooledSn.Add(sn);
            }

            Report("noise/median", pooledNoise.ToArray());
            Report("S/N", pooledSn.ToArray());
            Console.WriteLine($"slices used: {noiseBySlice.Count(kv => kv.Value.Count >= 200)}, pooled n: {pooledNoise.Count}");
        }
    }

    /// <summary>
    /// Compares a written simulated mzML against the real run it was built from, bin by bin. This is
    /// the check that the model survived rasterization, reduction, injection, merging, and the mzML
    /// writer — the unit tests only establish that the sampler itself is right.
    /// </summary>
    [Test]
    [Explicit("Compares a simulated mzML written by ExportRep2NoisySimulationForInspection to the real raw")]
    public static void CompareSimulatedMzmlToRealFile()
    {
        const double rtStart = 50.0;
        const double rtEnd = 55.0;

        string realPath = @"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw";
        string simulatedPath = @"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.rt45-56.noisy.simulated.mzML";

        if (!CanOpen(realPath) || !File.Exists(simulatedPath))
        {
            Console.WriteLine($"SKIPPED - real readable: {CanOpen(realPath)}, simulated exists: {File.Exists(simulatedPath)}");
            return;
        }

        var rawFile = RawFileReaderAdapter.FileFactory(realPath);
        Ms1Scan[] real;
        try
        {
            rawFile.SelectInstrument(Device.MS, 1);
            real = ReadMs1Window(rawFile, rtStart, rtEnd);
        }
        finally
        {
            rawFile.Dispose();
        }

        var simulatedFile = Readers.MsDataFileReader.GetDataFile(simulatedPath);
        simulatedFile.LoadAllStaticData();
        var simulated = simulatedFile.GetAllScansList()
            .Where(s => s.MsnOrder == 1 && s.RetentionTime >= rtStart && s.RetentionTime <= rtEnd)
            .Select(s => new Ms1Scan(s.OneBasedScanNumber, s.RetentionTime,
                s.MassSpectrum.XArray
                    .Select((mz, i) => new Centroid(mz, s.MassSpectrum.YArray[i], double.NaN, double.NaN, double.NaN, 0))
                    .ToArray()))
            .ToArray();

        Console.WriteLine($"real scans {real.Length} ({real.Average(s => s.Peaks.Length):F0} peaks/scan), " +
                          $"simulated scans {simulated.Length} ({simulated.Average(s => s.Peaks.Length):F0} peaks/scan)");

        Console.WriteLine();
        Console.WriteLine("peaks/scan and median intensity by 100-m/z bin");
        Console.WriteLine($"{"m/z",12} {"real n",10} {"sim n",10} {"n ratio",9} {"real I",12} {"sim I",12} {"I ratio",9}");

        for (int bin = 6; bin < 20; bin++)
        {
            double lo = bin * 100.0, hi = lo + 100.0;
            var realBin = real.SelectMany(s => s.Peaks).Where(p => p.Mz >= lo && p.Mz < hi).ToArray();
            var simBin = simulated.SelectMany(s => s.Peaks).Where(p => p.Mz >= lo && p.Mz < hi).ToArray();

            double realN = realBin.Length / (double)Math.Max(1, real.Length);
            double simN = simBin.Length / (double)Math.Max(1, simulated.Length);
            double realI = Median(realBin.Select(p => p.Intensity).DefaultIfEmpty(double.NaN).ToArray());
            double simI = Median(simBin.Select(p => p.Intensity).DefaultIfEmpty(double.NaN).ToArray());

            Console.WriteLine($"[{lo,5:F0}-{hi,5:F0}) {realN,10:F1} {simN,10:F1} {simN / realN,9:F2} " +
                              $"{realI,12:E2} {simI,12:E2} {simI / realI,9:F2}");
        }
    }

    /// <summary>
    /// Reported peaks per scan against retention time, over the whole run. The noise model uses one
    /// density everywhere, calibrated on the 50-55 min window; this is the measurement that says how
    /// wrong that is elsewhere.
    /// </summary>
    [Test]
    [Explicit("Prints the peaks-per-scan profile across a whole run")]
    public static void CharacterizeDensityAcrossRun()
    {
        const double binMinutes = 5.0;

        foreach (var (label, candidates) in Files)
        {
            string? path = candidates.FirstOrDefault(CanOpen);
            if (path is null)
            {
                Console.WriteLine($"{label}: SKIPPED (not readable)");
                continue;
            }

            var rawFile = RawFileReaderAdapter.FileFactory(path);
            try
            {
                rawFile.SelectInstrument(Device.MS, 1);
                var countsByBin = new SortedDictionary<int, List<int>>();

                for (int n = rawFile.RunHeaderEx.FirstSpectrum; n <= rawFile.RunHeaderEx.LastSpectrum; n++)
                {
                    var filter = rawFile.GetFilterForScanNumber(n);
                    if (filter is null || filter.MSOrder != ThermoFisher.CommonCore.Data.FilterEnums.MSOrderType.Ms)
                        continue;

                    var stream = rawFile.GetCentroidStream(n, false);
                    if (stream?.Masses is null)
                        continue;

                    int bin = (int)(rawFile.RetentionTimeFromScanNumber(n) / binMinutes);
                    if (!countsByBin.TryGetValue(bin, out var list))
                        countsByBin[bin] = list = new List<int>();
                    list.Add(stream.Masses.Length);
                }

                Console.WriteLine();
                Console.WriteLine($"---- {label}: peaks/scan vs RT ----");
                double overall = countsByBin.SelectMany(kv => kv.Value).Average();
                Console.WriteLine($"whole-run mean: {overall:F0} peaks/scan");
                foreach (var (bin, counts) in countsByBin)
                {
                    double median = Median(counts.Select(c => (double)c).ToArray());
                    Console.WriteLine($"  {bin * binMinutes,5:F0}-{(bin + 1) * binMinutes,5:F0} min  " +
                                      $"scans {counts.Count,5}  median {median,8:F0}  " +
                                      new string('#', (int)(median / 500)));
                }
            }
            finally
            {
                rawFile.Dispose();
            }
        }
    }

    private static void Characterize(string label, string path)
    {
        var rawFile = RawFileReaderAdapter.FileFactory(path);
        try
        {
            if (rawFile.IsError || !rawFile.IsOpen)
            {
                Console.WriteLine("could not open raw file");
                return;
            }

            rawFile.SelectInstrument(Device.MS, 1);
            int first = rawFile.RunHeaderEx.FirstSpectrum;
            int last = rawFile.RunHeaderEx.LastSpectrum;
            Console.WriteLine($"spectra: {first}-{last}, " +
                              $"RT {rawFile.RetentionTimeFromScanNumber(first):F2}-{rawFile.RetentionTimeFromScanNumber(last):F2} min");

            var busy = ReadMs1Window(rawFile, BusyRtStart, BusyRtEnd);
            var quiet = ReadMs1Window(rawFile, NoiseRtStart, NoiseRtEnd);

            ReportWindow(busy, $"BUSY {BusyRtStart}-{BusyRtEnd} min");
            ReportWindow(quiet, $"QUIET {NoiseRtStart}-{NoiseRtEnd} min");

            DumpPeaks(quiet, Path.Combine(ScratchDir(), $"{label}.rt50-55.peaks.tsv"));
        }
        finally
        {
            rawFile.Dispose();
        }
    }

    private static Ms1Scan[] ReadMs1Window(IRawDataPlus rawFile, double rtStart, double rtEnd)
    {
        int first = rawFile.ScanNumberFromRetentionTime(rtStart);
        int last = rawFile.ScanNumberFromRetentionTime(rtEnd);
        var scans = new List<Ms1Scan>();

        for (int n = first; n <= last; n++)
        {
            var filter = rawFile.GetFilterForScanNumber(n);
            if (filter is null || filter.MSOrder != ThermoFisher.CommonCore.Data.FilterEnums.MSOrderType.Ms)
                continue;

            var stream = rawFile.GetCentroidStream(n, false);
            if (stream?.Masses is null || stream.Intensities is null)
                continue;

            var peaks = new Centroid[stream.Masses.Length];
            for (int i = 0; i < peaks.Length; i++)
            {
                peaks[i] = new Centroid(
                    Mz: stream.Masses[i],
                    Intensity: stream.Intensities[i],
                    Noise: stream.Noises is { Length: > 0 } ? stream.Noises[i] : double.NaN,
                    Baseline: stream.Baselines is { Length: > 0 } ? stream.Baselines[i] : double.NaN,
                    Resolution: stream.Resolutions is { Length: > 0 } ? stream.Resolutions[i] : double.NaN,
                    Charge: stream.Charges is { Length: > 0 } ? (int)stream.Charges[i] : 0);
            }

            Array.Sort(peaks, (a, b) => a.Mz.CompareTo(b.Mz));
            scans.Add(new Ms1Scan(n, rawFile.RetentionTimeFromScanNumber(n), peaks));
        }

        return scans.ToArray();
    }

    private static void ReportWindow(Ms1Scan[] window, string label)
    {
        Console.WriteLine();
        Console.WriteLine($"---- {label} ----");
        if (window.Length == 0)
        {
            Console.WriteLine("no MS1 scans in window");
            return;
        }

        var counts = window.Select(s => (double)s.Peaks.Length).ToArray();
        Console.WriteLine($"MS1 scans: {window.Length}");
        Console.WriteLine($"peaks/scan: min {counts.Min():F0}, median {Median(counts):F0}, mean {counts.Average():F0}, max {counts.Max():F0}");

        var all = window.SelectMany(s => s.Peaks).ToArray();
        Report("intensity", all.Select(p => p.Intensity).ToArray());
        Report("noise", all.Select(p => p.Noise).Where(double.IsFinite).ToArray());
        Report("baseline", all.Select(p => p.Baseline).Where(double.IsFinite).ToArray());
        Report("S/N", all.Select(p => p.SignalToNoise).Where(double.IsFinite).ToArray());

        ReportByMzBin(window);
        ReportPersistence(window);
        ReportIsotopePartners(window);
        ReportChargeAssignments(all);
    }

    private static void Report(string name, double[] values)
    {
        if (values.Length == 0)
        {
            Console.WriteLine($"{name,-10}: (absent)");
            return;
        }

        var sorted = (double[])values.Clone();
        Array.Sort(sorted);
        Console.WriteLine($"{name,-10}: " + string.Join("  ",
            new[] { 0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99 }
                .Select(q => $"p{q * 100:F0}={Quantile(sorted, q):G4}")));
    }

    /// <summary>
    /// Peak density, noise level, and resolution per 100-m/z bin. An FT noise floor is roughly
    /// uniform in frequency, and frequency maps to m/z as f ∝ (m/z)^-1/2, so noise peak density
    /// should fall as (m/z)^-1.5 rather than staying flat. Resolution is reported because it pins
    /// the same exponent from the other side and calibrates OrbitrapPeakWidth.K.
    /// </summary>
    private static void ReportByMzBin(Ms1Scan[] window)
    {
        const double binWidth = 100.0;
        var byBin = new SortedDictionary<int, List<Centroid>>();

        foreach (var scan in window)
        foreach (var p in scan.Peaks)
        {
            int bin = (int)(p.Mz / binWidth);
            if (!byBin.TryGetValue(bin, out var list))
                byBin[bin] = list = new List<Centroid>();
            list.Add(p);
        }

        Console.WriteLine("by 100-m/z bin:  n/scan | median intensity | median noise | median S/N | median resolution");
        foreach (var (bin, peaks) in byBin)
        {
            Console.WriteLine($"  [{bin * binWidth,5:F0}-{(bin + 1) * binWidth,5:F0})  " +
                              $"{peaks.Count / (double)window.Length,8:F1}  " +
                              $"{Median(peaks.Select(p => p.Intensity).ToArray()),16:E2}  " +
                              $"{Median(peaks.Select(p => p.Noise).Where(double.IsFinite).DefaultIfEmpty(double.NaN).ToArray()),12:E2}  " +
                              $"{Median(peaks.Select(p => p.SignalToNoise).Where(double.IsFinite).DefaultIfEmpty(double.NaN).ToArray()),10:F2}  " +
                              $"{Median(peaks.Select(p => p.Resolution).Where(double.IsFinite).DefaultIfEmpty(double.NaN).ToArray()),17:F0}");
        }
    }

    /// <summary>
    /// Fraction of peaks that reappear at the same m/z in the next MS1 scan. Random FT noise does
    /// not persist; chemical background does. This is the single number that decides whether the
    /// noise model needs persistent species or only per-scan random draws.
    /// </summary>
    private static void ReportPersistence(Ms1Scan[] window)
    {
        const double ppmTolerance = 10.0;

        foreach (double quantile in new[] { 0.0, 0.5, 0.9 })
        {
            double matched = 0, total = 0;
            for (int s = 0; s + 1 < window.Length; s++)
            {
                var cur = window[s].Peaks;
                var next = window[s + 1].Peaks;
                if (next.Length == 0 || cur.Length == 0) continue;

                double floor = IntensityQuantile(cur, quantile);
                var nextMz = next.Select(p => p.Mz).ToArray();
                foreach (var p in cur)
                {
                    if (p.Intensity < floor) continue;
                    total++;
                    if (HasPeakWithin(nextMz, p.Mz, ppmTolerance))
                        matched++;
                }
            }

            string tag = quantile == 0 ? "all peaks" : $"top {(1 - quantile) * 100:F0}%";
            Console.WriteLine($"persistence to next scan ({tag}, {ppmTolerance:F0} ppm): " +
                              $"{(total > 0 ? matched / total : 0):P1}  (n={total:F0})");
        }
    }

    /// <summary>
    /// Fraction of peaks with a neighbour at +1.00335/z. Peaks with a z=1 or z=2 partner are
    /// small-molecule chemical background, which stresses a deconvolution algorithm quite
    /// differently from an unstructured noise floor.
    /// </summary>
    private static void ReportIsotopePartners(Ms1Scan[] window)
    {
        const double ppmTolerance = 10.0;
        const double neutronMass = 1.00335;

        var sampled = window.Where((_, i) => i % 5 == 0).ToArray();
        var partnerCounts = new int[8];
        int total = 0;

        foreach (var scan in sampled)
        {
            var mz = scan.Peaks.Select(p => p.Mz).ToArray();
            foreach (double m in mz)
            {
                total++;
                for (int z = 1; z <= 7; z++)
                {
                    if (HasPeakWithin(mz, m + neutronMass / z, ppmTolerance))
                        partnerCounts[z]++;
                }
            }
        }

        if (total == 0) return;
        Console.WriteLine("has a +1.00335/z neighbour: " + string.Join("  ",
            Enumerable.Range(1, 7).Select(z => $"z={z}:{partnerCounts[z] / (double)total:P1}")));
    }

    /// <summary>
    /// The instrument's own charge calls. A peak the instrument could not assign a charge to is
    /// one it could not find an isotope partner for, which is the operational definition of an
    /// unstructured noise peak.
    /// </summary>
    private static void ReportChargeAssignments(Centroid[] all)
    {
        if (all.Length == 0) return;
        var byCharge = all.GroupBy(p => Math.Min(p.Charge, 10))
                          .OrderBy(g => g.Key)
                          .ToArray();

        Console.WriteLine("instrument charge calls: " + string.Join("  ",
            byCharge.Select(g => $"z={(g.Key == 10 ? "10+" : g.Key.ToString())}:{g.Count() / (double)all.Length:P1}")));
    }

    private static void DumpPeaks(Ms1Scan[] window, string outPath)
    {
        using var writer = new StreamWriter(outPath);
        writer.WriteLine("ScanNumber\tRetentionTime\tMz\tIntensity\tNoise\tBaseline\tResolution\tCharge");
        foreach (var scan in window)
        foreach (var p in scan.Peaks)
        {
            writer.WriteLine(string.Join('\t',
                scan.ScanNumber.ToString(CultureInfo.InvariantCulture),
                R(scan.RetentionTime), R(p.Mz), R(p.Intensity), R(p.Noise), R(p.Baseline), R(p.Resolution),
                p.Charge.ToString(CultureInfo.InvariantCulture)));
        }

        Console.WriteLine($"peak dump: {outPath}");
    }

    private static string R(double v) => v.ToString("R", CultureInfo.InvariantCulture);

    private static bool CanOpen(string path)
    {
        if (!File.Exists(path)) return false;
        try
        {
            using var stream = File.OpenRead(path);
            return true;
        }
        catch (IOException)
        {
            return false;
        }
    }

    private static bool HasPeakWithin(double[] sortedMz, double target, double ppmTolerance)
    {
        double tolerance = target * ppmTolerance * 1e-6;
        int idx = Array.BinarySearch(sortedMz, target);
        if (idx >= 0) return true;
        idx = ~idx;

        if (idx < sortedMz.Length && Math.Abs(sortedMz[idx] - target) <= tolerance) return true;
        if (idx > 0 && Math.Abs(sortedMz[idx - 1] - target) <= tolerance) return true;
        return false;
    }

    private static double IntensityQuantile(Centroid[] peaks, double quantile)
    {
        if (quantile <= 0 || peaks.Length == 0) return 0;
        var sorted = peaks.Select(p => p.Intensity).ToArray();
        Array.Sort(sorted);
        return Quantile(sorted, quantile);
    }

    /// <summary>Quantile of an already-sorted array.</summary>
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

    private static string ScratchDir()
    {
        string dir = Path.Combine(Path.GetTempPath(), "TopDownSimulatorNoise");
        Directory.CreateDirectory(dir);
        return dir;
    }
}
