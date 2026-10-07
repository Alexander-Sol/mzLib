using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using TopDownSimulator.Comparison;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;
using TopDownSimulator.Noise;
using TopDownSimulator.Prediction;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

/// <summary>
/// Simulates a run from its identifications alone and scores the result against that run's raw
/// file, which the simulation itself never reads.
/// </summary>
/// <remarks>
/// <para>
/// The held-out run is Jurkat rep1 fract7, with rep2 fract7 as training data: same instrument,
/// method and fraction, different identifications. Its priors (<see cref="IdOnlyPriors"/>) are
/// learned from rep2's fitted models.
/// </para>
/// <para>
/// The scan grid and the per-scan noise come from a template run. <c>train</c> uses rep2, which is
/// the honest ID-only setting. <c>self</c> uses rep1's own scans, which leaks the held-out run's
/// noise and so isolates how much of the gap is the signal model and how much the noise template.
/// </para>
/// <para>
/// Run in order: <see cref="TrainPriors"/>, then <see cref="SimulateHeldOutFromIdsOnly"/> with each
/// template, then <see cref="CompareHeldOut"/>. The upper bound is rep1 fitted to its own raw
/// file, <c>AnalysisExample.ExportRep1FullNoisySimulation</c>, with MZLIB_TOPDOWN_SIM_OUTPUT_TAG=.v2.
/// </para>
/// </remarks>
[TestFixture]
public class IdOnlySimulation
{
    private const string Dir = @"D:\JurkatTopdown";
    private const string TrainStem = "02-18-20_jurkat_td_rep2_fract7";
    private const string TestStem = "02-18-20_jurkat_td_rep1_fract7";
    private const double QValue = 0.01;

    private static string PsmPath(string stem) =>
        $@"{Dir}\Frac7_GPTMD_Search\Task2-TopDownSearch\Individual File Results\{stem}_Proteoforms.psmtsv";

    /// <summary>
    /// The first copy of the raw file that can be opened. The top-level copies are sometimes held
    /// open by a viewer, and identical copies sit under Rep2_Raw.
    /// </summary>
    private static string RawPath(string stem)
    {
        foreach (string candidate in new[] { $@"{Dir}\{stem}.raw", $@"{Dir}\Rep2_Raw\{stem}.raw" })
        {
            try
            {
                using var stream = File.Open(candidate, FileMode.Open, FileAccess.Read, FileShare.Read);
                return candidate;
            }
            catch (IOException)
            {
            }
        }

        throw new IOException($"No readable copy of {stem}.raw.");
    }
    private static string FittedSidecar(string stem) => $@"{Dir}\{stem}.full.v2.noisy.simulated.groundtruth.tsv";
    private static string FittedMzml(string stem) => $@"{Dir}\{stem}.full.v2.noisy.simulated.mzML";
    private static string PriorsPath => $@"{Dir}\idonly-priors.{TrainStem}.json";
    private static string IdOnlyMzml(string template) => $@"{Dir}\{TestStem}.idonly-{template}.noisy.simulated.mzML";

    /// <summary>The retention times the per-scan comparison looks at, one per regime.</summary>
    private static readonly (double Rt, string Regime)[] TargetTimes =
    {
        (26.4, "pre-elution"), (35.1, "busy"), (41.1, "busy, histones"), (53.8, "wash"),
    };

    /// <summary>
    /// The template that supplies the scan grid and noise: <c>train</c> (default) or <c>self</c>.
    /// Set with MZLIB_TOPDOWN_SIM_IDONLY_TEMPLATE.
    /// </summary>
    private static string Template =>
        Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_IDONLY_TEMPLATE")?.Trim().ToLowerInvariant() == "self"
            ? "self"
            : "train";

    [Test]
    [Explicit("Learns ID-only priors from rep2 fract7's fitted models")]
    public static void TrainPriors()
    {
        var species = SpeciesGrouper.Group(new MmResultLoader().LoadQualified(PsmPath(TrainStem), TrainStem, QValue));
        var (models, _, _, width) = Simulator.ReadGroundTruth(FittedSidecar(TrainStem));
        var priors = IdOnlyPriors.Fit(species, models, $"{TrainStem}, {species.Length} species, width {width}");
        priors.Save(PriorsPath);

        Console.WriteLine($"species {species.Length}, fitted models {models.Length}, shape samples {priors.Shapes.Length}");
        Console.WriteLine($"rt apex offset {priors.RtApexOffset:F4} min");
        Console.WriteLine($"log10 abundance = {priors.LogAbundanceIntercept:F3} + {priors.LogAbundanceSlope:F4} log10(sum precursor intensity), fallback {priors.FallbackLogAbundance:F3}");
        Console.WriteLine($"charge mu = {priors.ChargeMuIntercept:F3} + {priors.ChargeMuChargeSlope:F4} mean charge + {priors.ChargeMuBasicSlope:F4} K+R");
        Console.WriteLine($"wrote {PriorsPath}");
    }

    [Test]
    [Explicit("Simulates rep1 fract7 from its identifications alone (template from MZLIB_TOPDOWN_SIM_IDONLY_TEMPLATE)")]
    public static void SimulateHeldOutFromIdsOnly()
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var priors = IdOnlyPriors.Load(PriorsPath);
        var species = SpeciesGrouper.Group(new MmResultLoader().LoadQualified(PsmPath(TestStem), TestStem, QValue));
        var models = priors.Predict(species);
        int minZ = species.Min(s => s.MinCharge), maxZ = species.Max(s => s.MaxCharge);

        // The width law is an instrument setting, so it is taken from the training run.
        var (_, _, _, width) = Simulator.ReadGroundTruth(FittedSidecar(TrainStem));

        string templateStem = Template == "self" ? TestStem : TrainStem;
        var templateScans = ReadMs1(RawPath(templateStem));
        var scanNoise = ScanNoiseConditions.FromSourceScans(templateScans, new NoiseFloorModel(663.0));
        double[] scanTimes = templateScans.Select(s => s.RetentionTime).ToArray();

        Console.WriteLine($"species {species.Length}, charges {minZ}-{maxZ}, template {templateStem} ({templateScans.Length} MS1 scans), loaded in {sw.Elapsed}");

        var export = new Simulator().WriteMzml(
            models, minZ, maxZ, width, scanTimes, IdOnlyMzml(Template),
            noise: new NoiseFloorModel(663.0), scanNoise: scanNoise);

        Console.WriteLine($"wrote {export.MzmlPath}: {export.ScanCount} scans, {export.PeakCount / (double)export.ScanCount:F0} peaks/scan, " +
                          $"{export.FeatureCount} features, signal {export.Noise!.SignalPeaks}, total {sw.Elapsed}");
    }

    [Test]
    [Explicit("Scores the ID-only rep1 simulations and the fitted upper bound against rep1's raw file")]
    public static void CompareHeldOut()
    {
        CompareParameters();

        var real = ReadMs1(RawPath(TestStem));
        var simulations = new List<(string Label, MsDataScan[] Scans)>
        {
            ("fitted (upper bound)", ReadMs1(FittedMzml(TestStem))),
        };
        foreach (string template in new[] { "self", "train" })
            if (File.Exists(IdOnlyMzml(template)))
                simulations.Add(($"ids only, {template} template", ReadMs1(IdOnlyMzml(template))));

        Console.WriteLine();
        Console.WriteLine("Whole run, simulation interpolated onto the real scan times:");
        Console.WriteLine($"  {"simulation",-28} {"log TIC r",10} {"peaks/scan r",13} {"TIC ratio",10} {"peak ratio",11}");
        foreach (var (label, scans) in simulations)
        {
            var matched = real.Select(r => Nearest(scans, r.RetentionTime)).ToArray();
            double ticR = Pearson(real.Select(r => Math.Log10(1 + r.MassSpectrum.SumOfAllY)), matched.Select(s => Math.Log10(1 + s.MassSpectrum.SumOfAllY)));
            double nR = Pearson(real.Select(r => (double)r.MassSpectrum.XArray.Length), matched.Select(s => (double)s.MassSpectrum.XArray.Length));
            double ticRatio = matched.Sum(s => s.MassSpectrum.SumOfAllY) / real.Sum(r => r.MassSpectrum.SumOfAllY);
            double nRatio = matched.Sum(s => (double)s.MassSpectrum.XArray.Length) / real.Sum(r => (double)r.MassSpectrum.XArray.Length);
            Console.WriteLine($"  {label,-28} {ticR,10:F3} {nR,13:F3} {ticRatio,10:F2} {nRatio,11:F2}");
        }

        Console.WriteLine();
        Console.WriteLine("Per scan (m/z 600-2000):");
        Console.WriteLine($"  {"RT",5} {"regime",-15} {"simulation",-28} {"real n",7} {"sim n",7} {"s match",7} {"cos",6} {"sqrtcos",7} {"TIC",7} {"KS",5}");
        foreach (var (rt, regime) in TargetTimes)
        {
            var r = Nearest(real, rt);
            foreach (var (label, scans) in simulations)
            {
                var s = Nearest(scans, rt);
                var c = CentroidSpectrumComparison.Compare(r.MassSpectrum.XArray, r.MassSpectrum.YArray,
                    s.MassSpectrum.XArray, s.MassSpectrum.YArray, 600, 2000);
                Console.WriteLine($"  {r.RetentionTime,5:F1} {regime,-15} {label,-28} {c.RealPeaks,7} {c.SimulatedPeaks,7} " +
                                  $"{c.SimulatedMatchedFraction,7:F2} {c.Cosine,6:F3} {c.SqrtCosine,7:F3} {c.TicRatio,7:F2} {c.RelativeIntensityKs,5:F2}");
            }
        }
    }

    /// <summary>How far each predicted parameter lands from rep1's own fit, species by species.</summary>
    private static void CompareParameters()
    {
        var priors = IdOnlyPriors.Load(PriorsPath);
        var species = SpeciesGrouper.Group(new MmResultLoader().LoadQualified(PsmPath(TestStem), TestStem, QValue));
        var predicted = priors.Predict(species).ToDictionary(m => m.Identifier!);
        var (fitted, _, _, _) = Simulator.ReadGroundTruth(FittedSidecar(TestStem));

        var pairs = fitted
            .Where(f => f.Abundance > 0 && f.Identifier is not null && predicted.ContainsKey(f.Identifier))
            .Select(f => (Fit: f, Pred: predicted[f.Identifier!]))
            .ToArray();

        Console.WriteLine($"Predicted vs fitted parameters for rep1 ({pairs.Length} species with a positive fitted abundance):");
        Report("rt mu error (min)", pairs.Select(p => p.Pred.RtProfile.Mu - p.Fit.RtProfile.Mu));
        Report("log10 abundance error", pairs.Select(p => Math.Log10(p.Pred.Abundance / p.Fit.Abundance)));
        Report("charge mu error", pairs.Select(p =>
            ((GaussianChargeDistribution)p.Pred.ChargeDistribution).MuZ - ((GaussianChargeDistribution)p.Fit.ChargeDistribution).MuZ));
        Console.WriteLine($"  log10 abundance: Pearson r {Pearson(pairs.Select(p => Math.Log10(p.Fit.Abundance)), pairs.Select(p => Math.Log10(p.Pred.Abundance))):F3}");
    }

    private static MsDataScan[] ReadMs1(string path)
    {
        var file = MsDataFileReader.GetDataFile(path);
        file.LoadAllStaticData();
        return file.GetAllScansList().Where(s => s.MsnOrder == 1).OrderBy(s => s.RetentionTime).ToArray();
    }

    private static MsDataScan Nearest(MsDataScan[] scans, double rt)
    {
        int lo = 0, hi = scans.Length - 1;
        while (lo < hi)
        {
            int mid = (lo + hi) / 2;
            if (scans[mid].RetentionTime < rt) lo = mid + 1;
            else hi = mid;
        }

        return lo > 0 && Math.Abs(scans[lo - 1].RetentionTime - rt) < Math.Abs(scans[lo].RetentionTime - rt)
            ? scans[lo - 1]
            : scans[lo];
    }

    private static double Pearson(IEnumerable<double> xs, IEnumerable<double> ys)
    {
        var x = xs.ToArray();
        var y = ys.ToArray();
        double mx = x.Average(), my = y.Average(), sxy = 0, sxx = 0, syy = 0;
        for (int i = 0; i < x.Length; i++)
        {
            sxy += (x[i] - mx) * (y[i] - my);
            sxx += (x[i] - mx) * (x[i] - mx);
            syy += (y[i] - my) * (y[i] - my);
        }

        return sxy / Math.Sqrt(sxx * syy);
    }

    private static void Report(string label, IEnumerable<double> values)
    {
        var v = values.Where(double.IsFinite).OrderBy(x => x).ToArray();
        double Q(double q) => v[(int)Math.Min(v.Length - 1, Math.Floor(q * v.Length))];
        Console.WriteLine($"  {label,-24} p5 {Q(0.05),7:F3}  p25 {Q(0.25),7:F3}  p50 {Q(0.5),7:F3}  p75 {Q(0.75),7:F3}  p95 {Q(0.95),7:F3}  " +
                          $"mean |err| {v.Average(Math.Abs),6:F3}");
    }
}
