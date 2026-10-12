using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using TopDownSimulator.Extraction;
using TopDownSimulator.Noise;

namespace Test.TopDownSimulator;

/// <summary>
/// How the noise differs between runs, and whether a run's identifications predict it. Writes one
/// table of per-scan noise descriptors and one of identifications per run, for analysis outside the
/// test (MZLIB_TOPDOWN_SIM_NOISE_RUNS_DIR, default D:\JurkatTopdown\noise-across-runs).
/// </summary>
[TestFixture]
public class NoiseAcrossRuns
{
    private const string Dir = @"D:\JurkatTopdown";
    private const string Search = @"D:\JurkatTopdown\Rep2_Raw\MM114_Search_50_50\Task1-SearchTask\Individual File Results";

    /// <summary>The ten rep2 fractions share one search; the rep1 runs are replicates of fract6 and fract7.</summary>
    private static readonly string[] Stems =
    {
        "02-17-20_jurkat_td_rep2_fract1", "02-17-20_jurkat_td_rep2_fract2", "02-17-20_jurkat_td_rep2_fract3",
        "02-17-20_jurkat_td_rep2_fract4", "02-18-20_jurkat_td_rep2_fract5", "02-18-20_jurkat_td_rep2_fract6",
        "02-18-20_jurkat_td_rep2_fract7", "02-18-20_jurkat_td_rep2_fract8", "02-18-20_jurkat_td_rep2_fract9",
        "02-18-20_jurkat_td_rep2_fract10", "02-18-20_jurkat_td_rep1_fract6", "02-18-20_jurkat_td_rep1_fract7",
    };

    [Test]
    [Explicit("Writes per-scan noise descriptors and identifications for twelve Jurkat runs; reads ~15 GB of raw data")]
    public static void ExportNoiseDescriptors()
    {
        string outDir = Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_NOISE_RUNS_DIR") ?? Path.Combine(Dir, "noise-across-runs");
        Directory.CreateDirectory(outDir);
        var only = Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_NOISE_RUNS_ONLY");
        var inv = CultureInfo.InvariantCulture;
        int bins = NoiseFloorModel.BinCount;

        foreach (string stem in Stems)
        {
            if (only is not null && !stem.Contains(only)) continue;
            string scansOut = Path.Combine(outDir, stem + ".scans.tsv");
            if (File.Exists(scansOut)) { Console.WriteLine($"{stem}: exists, skipped"); continue; }

            var watch = System.Diagnostics.Stopwatch.StartNew();
            var file = MsDataFileReader.GetDataFile(RawPath(stem));
            file.LoadAllStaticData();
            var all = file.GetAllScansList().OrderBy(s => s.OneBasedScanNumber).ToArray();
            var ms1 = all.Where(s => s.MsnOrder == 1).ToArray();

            // MS2 scans acquired after each MS1 scan, before the next: the instrument's own view of the load.
            var ms2After = new int[ms1.Length];
            int k = -1;
            foreach (var s in all)
            {
                if (s.MsnOrder == 1) k++;
                else if (k >= 0) ms2After[k]++;
            }

            var noise = ScanNoiseConditions.FromSourceScans(ms1, new NoiseFloorModel(NoiseFloorModel.JurkatNoiseLevelAtReferenceMz));
            var sb = new StringBuilder();
            sb.Append("scan\trt\tit\tms2\tpeaks\ttic\tnBright\tticBright\tnDim\tticDim");
            for (int b = 0; b < bins; b++) sb.Append($"\tdens{b}");
            for (int b = 0; b < bins; b++) sb.Append($"\tlvl{b}");
            for (int b = 0; b < bins; b++) sb.Append($"\tbright{b}");
            sb.AppendLine();

            for (int i = 0; i < ms1.Length; i++)
            {
                var s = ms1[i];
                var model = noise[i];
                var mz = s.MassSpectrum.XArray;
                var y = s.MassSpectrum.YArray;
                int nBright = 0, nDim = 0;
                double ticBright = 0, ticDim = 0;
                var brightTicPerBin = new double[bins];
                for (int j = 0; j < mz.Length; j++)
                {
                    if (y[j] >= ScanNoiseConditions.DefaultMaxSignalToNoise * model.NoiseLevelAt(mz[j]))
                    {
                        nBright++;
                        ticBright += y[j];
                        int bin = (int)Math.Floor((mz[j] - NoiseFloorModel.BinStart(0)) / 100.0);
                        if (bin >= 0 && bin < bins) brightTicPerBin[bin] += y[j];
                    }
                    else
                    {
                        nDim++;
                        ticDim += y[j];
                    }
                }

                sb.Append(string.Join('\t',
                    s.OneBasedScanNumber.ToString(inv), s.RetentionTime.ToString("F4", inv),
                    (s.InjectionTime ?? double.NaN).ToString("G5", inv), ms2After[i].ToString(inv),
                    mz.Length.ToString(inv), s.MassSpectrum.SumOfAllY.ToString("G6", inv),
                    nBright.ToString(inv), ticBright.ToString("G6", inv), nDim.ToString(inv), ticDim.ToString("G6", inv)));
                for (int b = 0; b < bins; b++) sb.Append('\t').Append(model.PeaksPerScanInBin(b).ToString(inv));
                for (int b = 0; b < bins; b++) sb.Append('\t').Append(model.NoiseLevelAt(NoiseFloorModel.BinStart(b) + 50).ToString("G5", inv));
                for (int b = 0; b < bins; b++) sb.Append('\t').Append(brightTicPerBin[b].ToString("G5", inv));
                sb.AppendLine();
            }

            File.WriteAllText(scansOut, sb.ToString());

            string psm = Path.Combine(Search, stem + "_Proteoforms.psmtsv");
            if (stem.Contains("rep2") && File.Exists(psm))
            {
                var ids = new MmResultLoader().LoadQualified(psm, stem, 0.01);
                var idb = new StringBuilder("ms2Scan\tprecursorScan\trt\tmass\tprecursorMass\tcharge\tprecursorIntensity\tscore\taccession\tsequence\n");
                foreach (var r in ids.OrderBy(r => r.RetentionTime))
                    idb.AppendLine(string.Join('\t', r.Ms2ScanNumber.ToString(inv), r.PrecursorScanNumber.ToString(inv),
                        r.RetentionTime.ToString("F4", inv), r.MonoisotopicMass.ToString("F4", inv),
                        (r.PrecursorMass ?? double.NaN).ToString("F4", inv), r.PrecursorCharge.ToString(inv),
                        (r.PrecursorIntensity ?? double.NaN).ToString("G6", inv), r.Score.ToString("F2", inv),
                        r.Accession ?? "", r.FullSequence));
                File.WriteAllText(Path.Combine(outDir, stem + ".ids.tsv"), idb.ToString());
                Console.WriteLine($"{stem}: {ms1.Length} MS1 scans, {ids.Count} IDs at q <= 0.01, {watch.Elapsed.TotalSeconds:F0} s");
            }
            else
            {
                Console.WriteLine($"{stem}: {ms1.Length} MS1 scans, no IDs from the shared search, {watch.Elapsed.TotalSeconds:F0} s");
            }

            file = null;
            GC.Collect();
        }
    }

    private static string RawPath(string stem)
    {
        string? copies = Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_RAW_COPY_DIR");
        foreach (string candidate in new[] { $@"{Dir}\Rep2_Raw\{stem}.raw", $@"{Dir}\{stem}.raw", $@"{copies}\{stem}.raw" })
        {
            try
            {
                using var stream = File.Open(candidate, FileMode.Open, FileAccess.Read, FileShare.Read);
                return candidate;
            }
            catch (IOException) { }
            catch (UnauthorizedAccessException) { }
        }

        throw new FileNotFoundException($"No readable copy of {stem}.raw");
    }
}
