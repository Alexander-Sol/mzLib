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
    /// <summary>
    /// The output tag of the fitted exports to train and score against, ".v2" by default. Set with
    /// MZLIB_TOPDOWN_SIM_FIT_TAG. Priors and ID-only exports from tags other than ".v2" carry the tag
    /// in their names, so earlier results stay readable.
    /// </summary>
    private static string FitTag =>
        Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_FIT_TAG")?.Trim() is { Length: > 0 } tag ? tag : ".v2";

    private static string DerivedTag => FitTag == ".v2" ? "" : FitTag;
    private static string FittedSidecar(string stem) => $@"{Dir}\{stem}.full{FitTag}.noisy.simulated.groundtruth.tsv";
    private static string FittedMzml(string stem) => $@"{Dir}\{stem}.full{FitTag}.noisy.simulated.mzML";
    private static string PriorsPath => $@"{Dir}\idonly-priors.{TrainStem}{DerivedTag}.json";
    private static string IdOnlyMzml(string template, ChargePredictor charge) =>
        $@"{Dir}\{TestStem}.idonly-{template}{ChargeSuffix(charge)}{DerivedTag}.noisy.simulated.mzML";

    private static string ChargeSuffix(ChargePredictor charge) => charge switch
    {
        ChargePredictor.SequenceOnly => "-seqcharge",
        ChargePredictor.ResidueWeights => "-rescharge",
        _ => "",
    };

    private static readonly ChargePredictor[] ChargePredictors =
        { ChargePredictor.ObservedCharges, ChargePredictor.SequenceOnly, ChargePredictor.ResidueWeights };

    /// <summary>
    /// Where the charge centre comes from: mass and a weight per charged residue (default), or with
    /// MZLIB_TOPDOWN_SIM_IDONLY_CHARGE=sequence, mass and net basic residues, or with =observed, the
    /// identifications' precursor charges.
    /// </summary>
    private static ChargePredictor Charge =>
        Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_IDONLY_CHARGE")?.Trim().ToLowerInvariant() switch
        {
            "sequence" => ChargePredictor.SequenceOnly,
            "observed" => ChargePredictor.ObservedCharges,
            _ => ChargePredictor.ResidueWeights,
        };

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
        Console.WriteLine($"sequence charge mu = {priors.SequenceChargeIntercept:F3} + {priors.SequenceChargeSqrtMassSlope:F4} sqrt mass + {priors.SequenceChargeNetBasicSlope:F4} (KRH - DE)");
        if (priors.ResidueCharge is { } r)
            Console.WriteLine($"residue charge mu = {r.Intercept:F3} + {r.SqrtMassSlope:F4} sqrt mass" +
                              string.Concat(r.Residues.Select((c, i) => $" {r.Weights[i]:+0.0000;-0.0000} {c}")) + $", residual sd {r.ResidualSd:F3}");
        if (priors.ChargeSigmaIntercept is { } sa)
            Console.WriteLine($"charge sigma = {sa:F3} + {priors.ChargeSigmaMuSlope:F4} charge mu");
        Console.WriteLine($"median charge sigma {priors.Shapes.Select(s => s.ChargeSigma).OrderBy(x => x).ElementAt(priors.Shapes.Length / 2):F3}");
        Console.WriteLine($"wrote {PriorsPath}");
    }

    [Test]
    [Explicit("Simulates rep1 fract7 from its identifications alone (template from MZLIB_TOPDOWN_SIM_IDONLY_TEMPLATE)")]
    public static void SimulateHeldOutFromIdsOnly()
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var priors = IdOnlyPriors.Load(PriorsPath);
        var species = SpeciesGrouper.Group(new MmResultLoader().LoadQualified(PsmPath(TestStem), TestStem, QValue));
        var models = priors.Predict(species, chargePredictor: Charge);

        // Without precursor charges the range has to come from the predicted distributions too.
        int minZ, maxZ;
        if (Charge == ChargePredictor.ObservedCharges)
        {
            minZ = species.Min(s => s.MinCharge);
            maxZ = species.Max(s => s.MaxCharge);
        }
        else
        {
            var charges = models.Select(m => (GaussianChargeDistribution)m.ChargeDistribution).ToArray();
            minZ = Math.Max(2, (int)Math.Floor(charges.Min(c => c.MuZ - 3 * c.SigmaZ)));
            maxZ = Math.Min(80, (int)Math.Ceiling(charges.Max(c => c.MuZ + 3 * c.SigmaZ)));
        }

        // The width law is an instrument setting, so it is taken from the training run.
        var (_, _, _, width) = Simulator.ReadGroundTruth(FittedSidecar(TrainStem));

        string templateStem = Template == "self" ? TestStem : TrainStem;
        var templateScans = ReadMs1(RawPath(templateStem));
        var scanNoise = ScanNoiseConditions.FromSourceScans(templateScans, new NoiseFloorModel(663.0));
        double[] scanTimes = templateScans.Select(s => s.RetentionTime).ToArray();

        Console.WriteLine($"species {species.Length}, charge predictor {Charge}, charges {minZ}-{maxZ}, template {templateStem} ({templateScans.Length} MS1 scans), loaded in {sw.Elapsed}");

        var export = new Simulator().WriteMzml(
            models, minZ, maxZ, width, scanTimes, IdOnlyMzml(Template, Charge),
            noise: new NoiseFloorModel(663.0), scanNoise: scanNoise);

        Console.WriteLine($"wrote {export.MzmlPath}: {export.ScanCount} scans, {export.PeakCount / (double)export.ScanCount:F0} peaks/scan, " +
                          $"{export.FeatureCount} features, signal {export.Noise!.SignalPeaks}, total {sw.Elapsed}");
    }

    [Test]
    [Explicit("Scores the ID-only rep1 simulations and the fitted upper bound against rep1's raw file")]
    public static void CompareHeldOut()
    {
        foreach (var charge in ChargePredictors)
            CompareParameters(charge);

        var real = ReadMs1(RawPath(TestStem));
        var simulations = new List<(string Label, MsDataScan[] Scans)>
        {
            ("fitted (upper bound)", ReadMs1(FittedMzml(TestStem))),
        };
        foreach (var charge in ChargePredictors)
        foreach (string template in new[] { "self", "train" })
            if (File.Exists(IdOnlyMzml(template, charge)))
                simulations.Add(($"ids, {template}{ChargeSuffix(charge)}", ReadMs1(IdOnlyMzml(template, charge))));

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

    /// <summary>
    /// Scores each model's charge distribution against the charge profile observed in rep1's raw
    /// data, species by species. Whole-scan metrics cannot see this: they are dominated by noise
    /// and by the tallest peaks.
    /// </summary>
    /// <remarks>
    /// For each species the extractor sums its isotopologue intensities per charge over the anchor's
    /// retention time ± 0.25 min, across charges 5 to 40. That observed profile is compared with
    /// f(z) of the fitted model and of both ID-only predictions, by cosine and by whether the most
    /// intense charge lands within one of the observed one. Random peaks inside the extraction windows
    /// add a floor to every observed profile, the same for every model.
    /// </remarks>
    [Test]
    [Explicit("Scores fitted and predicted charge distributions against rep1's observed charge profiles")]
    public static void CompareChargeProfiles()
    {
        const int minZ = 5, maxZ = 40;
        var priors = IdOnlyPriors.Load(PriorsPath);
        var species = SpeciesGrouper.Group(new MmResultLoader().LoadQualified(PsmPath(TestStem), TestStem, QValue));
        var bySpecies = species.ToDictionary(s => s.Anchor.Identifier);
        var observedCharge = priors.Predict(species).ToDictionary(m => m.Identifier!);
        var sequenceCharge = priors.Predict(species, chargePredictor: ChargePredictor.SequenceOnly).ToDictionary(m => m.Identifier!);
        var residueCharge = priors.Predict(species, chargePredictor: ChargePredictor.ResidueWeights).ToDictionary(m => m.Identifier!);
        var (fitted, _, _, _) = Simulator.ReadGroundTruth(FittedSidecar(TestStem));
        Console.WriteLine($"fitted {FittedSidecar(TestStem)}, priors {PriorsPath}" +
                          (priors.ResidueCharge is null ? " (no residue model; residue weights = sequence only)" : ""));

        var extractor = new GroundTruthExtractor(ReadMs1(RawPath(TestStem)), ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
        string[] labels = { "fitted to rep1 raw", "ids, observed charges", "ids, sequence only", "ids, residue weights" };
        var apexOffsets = labels.Append("mean member precursor charge").ToDictionary(l => l, _ => new List<int>());
        var scores = labels.ToDictionary(l => l, _ => new List<(double Cosine, bool ApexWithinOne)>());

        foreach (var fit in fitted.Where(f => f.Abundance > 0 && f.Identifier is not null && bySpecies.ContainsKey(f.Identifier)))
        {
            var anchor = bySpecies[fit.Identifier!].Anchor;
            var truth = extractor.Extract(anchor.MonoisotopicMass, anchor.RetentionTime, 0.25, minZ, maxZ);
            var observed = truth.ChargeXics.Select(xic => xic.Sum()).ToArray();
            if (observed.Sum() <= 0) continue;
            int observedApex = minZ + Array.IndexOf(observed, observed.Max());

            foreach (var (label, model) in new[]
                     {
                         ("fitted to rep1 raw", fit),
                         ("ids, observed charges", observedCharge[fit.Identifier!]),
                         ("ids, sequence only", sequenceCharge[fit.Identifier!]),
                         ("ids, residue weights", residueCharge[fit.Identifier!]),
                     })
            {
                var predicted = Enumerable.Range(minZ, maxZ - minZ + 1).Select(z => model.ChargeDistribution.Evaluate(z)).ToArray();
                int predictedApex = minZ + Array.IndexOf(predicted, predicted.Max());
                scores[label].Add((SpectralAngle.ComputeCosineSimilarity(observed, predicted), Math.Abs(predictedApex - observedApex) <= 1));
                apexOffsets[label].Add(predictedApex - observedApex);
            }

            apexOffsets["mean member precursor charge"].Add(
                (int)Math.Round(IdOnlyPriors.MeanPrecursorCharge(bySpecies[fit.Identifier!])) - observedApex);
        }

        Console.WriteLine("Most intense charge minus the observed one (negative = predicted too low):");
        foreach (var (label, offsets) in apexOffsets)
        {
            var o = offsets.OrderBy(x => x).ToArray();
            Console.WriteLine($"  {label,-30} p10 {o[o.Length / 10],3}  p25 {o[o.Length / 4],3}  median {o[o.Length / 2],3}  p75 {o[3 * o.Length / 4],3}  p90 {o[9 * o.Length / 10],3}");
        }

        Console.WriteLine($"Charge profiles of {scores.First().Value.Count} rep1 species, charges {minZ}-{maxZ}:");
        Console.WriteLine($"  {"model",-24} {"median cos",11} {"p25 cos",8} {"apex within 1",14}");
        foreach (var (label, list) in scores)
        {
            var cos = list.Select(x => x.Cosine).OrderBy(x => x).ToArray();
            Console.WriteLine($"  {label,-24} {cos[cos.Length / 2],11:F3} {cos[cos.Length / 4],8:F3} {list.Count(x => x.ApexWithinOne) / (double)list.Count,14:P0}");
        }
    }

    /// <summary>
    /// What the charge fitter sees over a wide window: each species' per-charge XIC apex, as the
    /// fitter takes it, its floor relative to its maximum, and what the trimmed fit makes of it.
    /// </summary>
    [Test]
    [Explicit("Prints per-charge apex profiles of rep1 species over charges 5-40")]
    public static void ChargeApexProfiles()
    {
        const int minZ = 5, maxZ = 40;
        var species = SpeciesGrouper.Group(new MmResultLoader().LoadQualified(PsmPath(TestStem), TestStem, QValue));
        var extractor = new GroundTruthExtractor(ReadMs1(RawPath(TestStem)), ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
        var floors = new List<double>();
        int shown = 0;
        foreach (var s in species.Where((_, i) => i % 5 == 0))
        {
            var truth = extractor.Extract(s.Anchor.MonoisotopicMass, s.Anchor.RetentionTime, 0.25, minZ, maxZ);
            var apex = truth.ChargeXics.Select(x => x.Max()).ToArray();
            double max = apex.Max();
            if (max <= 0) continue;
            double floor = apex.OrderBy(a => a).ElementAt(apex.Length / 5) / max;
            floors.Add(floor);
            if (shown++ >= 12) continue;
            var fit = new global::TopDownSimulator.Fitting.ChargeDistributionFitter(trimFraction: 0.05).Fit(truth);
            Console.WriteLine($"{s.Anchor.MonoisotopicMass,9:F1} members {s.MinCharge + 2}-{s.MaxCharge - 2}  fit mu {fit.Distribution.MuZ,5:F1} sigma {fit.Distribution.SigmaZ,4:F1} used {fit.ChargesUsed,2}  p20/max {floor:F3}");
            Console.WriteLine("   " + string.Join(" ", apex.Select(a => ((int)Math.Round(99 * a / max)).ToString().PadLeft(2))));
        }

        var f = floors.OrderBy(x => x).ToArray();
        Console.WriteLine($"p20 apex / max over {f.Length} species: p25 {f[f.Length / 4]:F3}  median {f[f.Length / 2]:F3}  p75 {f[3 * f.Length / 4]:F3}");
    }

    /// <summary>How far each predicted parameter lands from rep1's own fit, species by species.</summary>
    private static void CompareParameters(ChargePredictor chargePredictor)
    {
        var priors = IdOnlyPriors.Load(PriorsPath);
        var species = SpeciesGrouper.Group(new MmResultLoader().LoadQualified(PsmPath(TestStem), TestStem, QValue));
        var predicted = priors.Predict(species, chargePredictor: chargePredictor).ToDictionary(m => m.Identifier!);
        var (fitted, _, _, _) = Simulator.ReadGroundTruth(FittedSidecar(TestStem));

        var pairs = fitted
            .Where(f => f.Abundance > 0 && f.Identifier is not null && predicted.ContainsKey(f.Identifier))
            .Select(f => (Fit: f, Pred: predicted[f.Identifier!]))
            .ToArray();

        Console.WriteLine($"Predicted ({chargePredictor}) vs fitted parameters for rep1 ({pairs.Length} species with a positive fitted abundance):");
        Report("rt mu error (min)", pairs.Select(p => p.Pred.RtProfile.Mu - p.Fit.RtProfile.Mu));
        Report("log10 abundance error", pairs.Select(p => Math.Log10(p.Pred.Abundance / p.Fit.Abundance)));
        Report("charge mu error", pairs.Select(p =>
            ((GaussianChargeDistribution)p.Pred.ChargeDistribution).MuZ - ((GaussianChargeDistribution)p.Fit.ChargeDistribution).MuZ));
        Console.WriteLine($"  log10 abundance: Pearson r {Pearson(pairs.Select(p => Math.Log10(p.Fit.Abundance)), pairs.Select(p => Math.Log10(p.Pred.Abundance))):F3}");

        // The species that dominate the signal-dominated scan, where errors decide its cosine.
        Console.WriteLine("  brightest species eluting at 41.1 min (fitted | predicted): mass, log10 A, mu_z, sigma_z, rt mu");
        foreach (var (fit, pred) in pairs.Where(p => Math.Abs(p.Fit.RtProfile.Mu - 41.1) < 0.4)
                     .OrderByDescending(p => Math.Max(p.Fit.Abundance, p.Pred.Abundance)).Take(8))
        {
            var fz = (GaussianChargeDistribution)fit.ChargeDistribution;
            var pz = (GaussianChargeDistribution)pred.ChargeDistribution;
            Console.WriteLine($"    {fit.MonoisotopicMass,9:F1}  {Math.Log10(fit.Abundance),5:F2} | {Math.Log10(pred.Abundance),5:F2}  " +
                              $"{fz.MuZ,5:F1} | {pz.MuZ,5:F1}  {fz.SigmaZ,4:F1} | {pz.SigmaZ,4:F1}  {fit.RtProfile.Mu,6:F2} | {pred.RtProfile.Mu,6:F2}");
        }
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
