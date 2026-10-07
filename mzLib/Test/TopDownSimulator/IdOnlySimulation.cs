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
    private const double QValue = 0.01;

    /// <summary>A run to simulate from its identifications, and where its raw file and fitted exports live.</summary>
    private sealed record HeldOutRun(string Stem, string Directory, string PsmPath);

    /// <summary>
    /// The held-out run: rep1 fract7 (default; same fraction, GPTMD search) or, with
    /// MZLIB_TOPDOWN_SIM_IDONLY_HELDOUT=rep2-fract6 or rep2-fract5, another fraction of the training
    /// replicate, identified by MM114_Search_50_50 (Classic deconvolution, no GPTMD). Their fitted
    /// upper bounds come from <c>AnalysisExample.ExportRep2OtherFractionFullNoisySimulation</c>.
    /// </summary>
    private static HeldOutRun HeldOut =>
        Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_IDONLY_HELDOUT")?.Trim().ToLowerInvariant() switch
        {
            "rep2-fract6" => OtherFraction("fract6"),
            "rep2-fract5" => OtherFraction("fract5"),
            _ => new HeldOutRun("02-18-20_jurkat_td_rep1_fract7", Dir, PsmPath("02-18-20_jurkat_td_rep1_fract7")),
        };

    private static HeldOutRun OtherFraction(string fraction)
    {
        string stem = $"02-18-20_jurkat_td_rep2_{fraction}";
        return new HeldOutRun(stem, $@"{Dir}\Rep2_Raw",
            $@"{Dir}\Rep2_Raw\MM114_Search_50_50\Task1-SearchTask\Individual File Results\{stem}_Proteoforms.psmtsv");
    }

    private static string TestStem => HeldOut.Stem;

    private static string PsmPath(string stem) =>
        stem == TrainStem || stem == "02-18-20_jurkat_td_rep1_fract7"
            ? $@"{Dir}\Frac7_GPTMD_Search\Task2-TopDownSearch\Individual File Results\{stem}_Proteoforms.psmtsv"
            : HeldOut.PsmPath;

    private static string RunDir(string stem) => stem == TrainStem ? Dir : HeldOut.Directory;

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
    private static string FittedSidecar(string stem) => $@"{RunDir(stem)}\{stem}.full{FitTag}.noisy.simulated.groundtruth.tsv";
    private static string FittedMzml(string stem) => $@"{RunDir(stem)}\{stem}.full{FitTag}.noisy.simulated.mzML";
    private static string PriorsPath => $@"{Dir}\idonly-priors.{TrainStem}{DerivedTag}.json";
    private static string IdOnlyMzml(string template, ChargePredictor charge) =>
        $@"{HeldOut.Directory}\{TestStem}.idonly-{template}{ChargeSuffix(charge)}{DerivedTag}.noisy.simulated.mzML";

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
        foreach (string template in new[] { "self", "train", "notemplate" })
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

        CompareSignalPeaks(real, simulations);
        CompareSpeciesEnvelopes(real, simulations);
    }

    /// <summary>
    /// Per scan, only the peaks at S/N ≥ 10 over m/z 600–2000, at the target times and summarised over
    /// every 10th scan from 30 to 50 min. The noise level is the real scan's, from its injection time
    /// (<see cref="ScanNoiseConditions"/>), and the simulated scan is cut at the same threshold.
    /// </summary>
    private static void CompareSignalPeaks(MsDataScan[] real, List<(string Label, MsDataScan[] Scans)> simulations)
    {
        var noiseByScan = ScanNoiseConditions.FromSourceScans(real, new NoiseFloorModel(663.0), conditionDensity: false)
            .Select((m, i) => (real[i], m))
            .ToDictionary(p => p.Item1, p => p.m);

        CentroidSpectrumComparison Compare(MsDataScan r, MsDataScan s)
        {
            var noise = noiseByScan[r];
            var (rMz, rI) = SignalPeaks.Select(r.MassSpectrum.XArray, r.MassSpectrum.YArray, noise, 10, 600, 2000);
            var (sMz, sI) = SignalPeaks.Select(s.MassSpectrum.XArray, s.MassSpectrum.YArray, noise, 10, 600, 2000);
            return CentroidSpectrumComparison.Compare(rMz, rI, sMz, sI);
        }

        Console.WriteLine();
        Console.WriteLine("Per scan, peaks at S/N >= 10 only (m/z 600-2000):");
        Console.WriteLine($"  {"RT",5} {"regime",-15} {"simulation",-28} {"real n",7} {"sim n",7} {"r match",7} {"s match",7} {"cos",6} {"sqrtcos",7} {"TIC",6}");
        foreach (var (rt, regime) in TargetTimes)
        {
            var r = Nearest(real, rt);
            foreach (var (label, scans) in simulations)
            {
                var c = Compare(r, Nearest(scans, rt));
                Console.WriteLine($"  {r.RetentionTime,5:F1} {regime,-15} {label,-28} {c.RealPeaks,7} {c.SimulatedPeaks,7} " +
                                  $"{c.RealMatchedFraction,7:F2} {c.SimulatedMatchedFraction,7:F2} {c.Cosine,6:F3} {c.SqrtCosine,7:F3} {c.TicRatio,6:F2}");
            }
        }

        var sample = real.Where(r => r.RetentionTime is >= 30 and <= 50).Where((_, i) => i % 10 == 0).ToArray();
        Console.WriteLine($"  Medians over {sample.Length} scans, 30-50 min:");
        Console.WriteLine($"  {"simulation",-28} {"count ratio",11} {"r match",7} {"s match",7} {"cos",6} {"sqrtcos",7} {"TIC",6}");
        foreach (var (label, scans) in simulations)
        {
            var cs = sample.Select(r => Compare(r, Nearest(scans, r.RetentionTime))).ToArray();
            Console.WriteLine($"  {label,-28} {Median(cs.Select(c => c.PeakCountRatio)),11:F2} {Median(cs.Select(c => c.RealMatchedFraction)),7:F2} " +
                              $"{Median(cs.Select(c => c.SimulatedMatchedFraction)),7:F2} {Median(cs.Select(c => c.Cosine)),6:F3} " +
                              $"{Median(cs.Select(c => c.SqrtCosine)),7:F3} {Median(cs.Select(c => c.TicRatio)),6:F2}");
        }
    }

    /// <summary>
    /// Per species, the real and simulated envelopes around the apex rep1's own fit places it at:
    /// isotopologue intensities per charge, summed over apex ± 0.1 min, over every charge at which the
    /// species lands inside the scan window. Judges the signal model species by species, where the
    /// whole-scan metrics only see the tallest peaks and the noise.
    /// </summary>
    private static void CompareSpeciesEnvelopes(MsDataScan[] real, List<(string Label, MsDataScan[] Scans)> simulations)
    {
        const double halfWidth = 0.1;
        var species = SpeciesGrouper.Group(new MmResultLoader().LoadQualified(PsmPath(TestStem), TestStem, QValue));
        var anchors = species.Select(s => s.Anchor.Identifier).ToHashSet();
        var (fitted, _, _, _) = Simulator.ReadGroundTruth(FittedSidecar(TestStem));
        var realExtractor = new GroundTruthExtractor(real, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
        var targets = fitted
            .Where(f => f.Abundance > 0 && f.Identifier is not null && anchors.Contains(f.Identifier))
            .Select(f => (f.MonoisotopicMass, Apex: f.RtProfile.Mu, Charges: realExtractor.ObservableCharges(f.MonoisotopicMass)))
            .ToArray();
        var realTruths = targets.AsParallel().AsOrdered()
            .Select(t => realExtractor.Extract(t.MonoisotopicMass, t.Apex, halfWidth, t.Charges.Min, t.Charges.Max))
            .ToArray();

        Console.WriteLine();
        Console.WriteLine($"Per species ({targets.Length}), envelope summed over the fitted apex ± {halfWidth} min, charges inside the scan window:");
        Console.WriteLine($"  {"simulation",-28} {"n",4} {"median cos",10} {"p25 cos",8} {"charge cos",10} {"apex within 1",13} {"log10 ratio p25/p50/p75",24} {"mean |log10 ratio|",18}");
        foreach (var (label, scans) in simulations)
        {
            var extractor = new GroundTruthExtractor(scans, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
            var results = targets.AsParallel().AsOrdered()
                .Select((t, i) => SpeciesEnvelopeComparison.Compare(realTruths[i],
                    extractor.Extract(t.MonoisotopicMass, t.Apex, halfWidth, t.Charges.Min, t.Charges.Max)))
                .Where(c => c is not null)
                .Select(c => c!)
                .ToArray();
            var cos = results.Select(c => c.Cosine).OrderBy(x => x).ToArray();
            var ratio = results.Select(c => Math.Log10(c.IntensityRatio)).OrderBy(x => x).ToArray();
            Console.WriteLine($"  {label,-28} {results.Length,4} {cos[cos.Length / 2],10:F3} {cos[cos.Length / 4],8:F3} " +
                              $"{Median(results.Select(c => c.ChargeProfileCosine)),10:F3} " +
                              $"{results.Count(c => Math.Abs(c.RealApexCharge - c.SimulatedApexCharge) <= 1) / (double)results.Length,13:P0} " +
                              $"{$"{ratio[ratio.Length / 4]:F2} / {ratio[ratio.Length / 2]:F2} / {ratio[3 * ratio.Length / 4]:F2}",24} {ratio.Average(Math.Abs),18:F3}");
        }
    }

    private static double Median(IEnumerable<double> values)
    {
        var v = values.Where(double.IsFinite).OrderBy(x => x).ToArray();
        return v.Length == 0 ? double.NaN : v[v.Length / 2];
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

    private static string AcquisitionProfilePath => $@"{Dir}\acquisition-profile.{TrainStem}.json";
    private static string UnidentifiedPath => $@"{Dir}\unidentified.{TrainStem}{DerivedTag}.json";

    /// <summary>
    /// Calibrates the unidentified-analyte component on rep2: renders rep2's fitted models plus
    /// <see cref="UnidentifiedAnalytes"/> drawn from them, without noise, on every 5th MS1 scan from
    /// 30 to 50 min, and picks the count and abundance shift whose summed TIC and number of peaks at
    /// S/N ≥ 10 (against each real scan's injection-time noise) best match the real scans.
    /// Also learns and saves the <see cref="AcquisitionProfile"/>.
    /// </summary>
    [Test]
    [Explicit("Calibrates unidentified analytes and learns the acquisition profile on rep2 fract7")]
    public static void CalibrateTemplateFreeModel()
    {
        var real = ReadMs1(RawPath(TrainStem));
        var profile = AcquisitionProfile.Learn(real, TrainStem);
        profile.Save(AcquisitionProfilePath);
        Console.WriteLine($"AGC: IT = min({profile.Agc.MaxInjectionTime:F0}, {profile.Agc.ChargeTarget:G3} / (TIC + {profile.Agc.BackgroundTic:G2})), load below S/N 10 = {profile.DimLoadRatio:F2} x load above; " +
                          $"{profile.ScanTimes().Length} scan times from the profile against {real.Length} real; wrote {AcquisitionProfilePath}");

        var (fitted, minZ, maxZ, width) = Simulator.ReadGroundTruth(FittedSidecar(TrainStem));
        var sample = real.Where(s => s.RetentionTime is >= 30 and <= 50).Where((_, i) => i % 5 == 0).ToArray();
        double[] times = sample.Select(s => s.RetentionTime).ToArray();
        var noise = sample.Select(s => new NoiseFloorModel(NoiseFloorModel.JurkatNoiseTimesInjectionTime / (s.InjectionTime ?? 50))).ToArray();
        int RealBright(int i) => SignalPeaks.Select(sample[i].MassSpectrum.XArray, sample[i].MassSpectrum.YArray, noise[i], 10, 600, 2000).Mz.Length;
        var realBright = Enumerable.Range(0, sample.Length).Select(RealBright).ToArray();

        // Where the real TIC sits: above or below S/N 10, against what the noise model would put there.
        var realBrightTic = Enumerable.Range(0, sample.Length).Select(i =>
            SignalPeaks.Select(sample[i].MassSpectrum.XArray, sample[i].MassSpectrum.YArray, noise[i], 10, 600, 2000).Intensity.Sum()).ToArray();
        var sourceNoise = ScanNoiseConditions.FromSourceScans(sample, new NoiseFloorModel());
        var rng = new Random(1);
        var modelNoiseTic = sourceNoise.Select(m =>
        {
            var peaks = new List<SimulatedPeak>();
            m.SampleScan(rng, peaks);
            return peaks.Sum(p => p.Intensity);
        }).ToArray();
        Console.WriteLine($"real TIC median {Median(sample.Select(s => s.MassSpectrum.SumOfAllY)):G3}; at S/N >= 10 {Median(realBrightTic):G3}; " +
                          $"noise model sampled on each scan's own density {Median(modelNoiseTic):G3}; real bright peaks {Median(realBright.Select(b => (double)b)):F0}");

        var simulator = new Simulator();
        (double LogTic, double LogCount) Score(int count, double shift)
        {
            var background = new UnidentifiedAnalytes(count, shift).Draw(fitted, real[0].RetentionTime, real[^1].RetentionTime, seed: 1);
            var scans = simulator.SimulateCentroided(fitted.Concat(background).ToArray(), minZ, maxZ, width, times).Scans;
            var selected = Enumerable.Range(0, sample.Length)
                .Select(i => SignalPeaks.Select(scans[i].MassSpectrum.XArray, scans[i].MassSpectrum.YArray, noise[i], 10, 600, 2000)).ToArray();
            var tic = Enumerable.Range(0, sample.Length).Select(i => Math.Log10((1 + selected[i].Intensity.Sum()) / (1 + realBrightTic[i]))).ToArray();
            var bright = Enumerable.Range(0, sample.Length).Select(i => Math.Log10((1 + selected[i].Mz.Length) / (1.0 + realBright[i]))).ToArray();
            return (Median(tic), Median(bright));
        }

        // Below S/N 10 the acquisition profile's density curve already carries the load: it was
        // measured as every peak under S/N 10, unidentified analytes included. So the unidentified
        // component is calibrated on the bright part of the spectrum only.
        Console.WriteLine($"{sample.Length} scans, 30-50 min; median log10(sim/real) of the TIC and the number of peaks at S/N >= 10:");
        var best = (Count: 0, Shift: 0.0, Error: double.PositiveInfinity);
        foreach (int count in new[] { 0, 500, 1000, 2000, 4000, 8000 })
        foreach (double shift in count == 0 ? new[] { 0.0 } : new[] { -2.0, -1.5, -1.0, -0.5, 0.0 })
        {
            var (tic, bright) = Score(count, shift);
            double error = tic * tic + bright * bright;
            Console.WriteLine($"  count {count,5} shift {shift,5:F1}:  TIC {tic,6:F3}  bright peaks {bright,6:F3}");
            if (error < best.Error) best = (count, shift, error);
        }

        var chosen = new UnidentifiedAnalytes(best.Count, best.Shift);
        File.WriteAllText(UnidentifiedPath, System.Text.Json.JsonSerializer.Serialize(chosen));
        Console.WriteLine($"chosen {chosen}; wrote {UnidentifiedPath}");
    }

    /// <summary>
    /// Simulates the held-out run from its identifications with no template run at all: scan times
    /// and noise density from the rep2 acquisition profile, noise amplitude from AGC on the simulated
    /// ion load, and unidentified analytes calibrated on rep2. Run
    /// <see cref="CalibrateTemplateFreeModel"/> first. Written as the <c>notemplate</c> template.
    /// </summary>
    [Test]
    [Explicit("Simulates the held-out run from IDs only, without any template run")]
    public static void SimulateHeldOutWithoutTemplate()
    {
        var sw = System.Diagnostics.Stopwatch.StartNew();
        var priors = IdOnlyPriors.Load(PriorsPath);
        var profile = AcquisitionProfile.Load(AcquisitionProfilePath);
        var unidentified = System.Text.Json.JsonSerializer.Deserialize<UnidentifiedAnalytes>(File.ReadAllText(UnidentifiedPath))!;
        var species = SpeciesGrouper.Group(new MmResultLoader().LoadQualified(PsmPath(TestStem), TestStem, QValue));
        var models = priors.Predict(species, chargePredictor: Charge);
        var (fitted, _, _, width) = Simulator.ReadGroundTruth(FittedSidecar(TrainStem));
        double[] scanTimes = profile.ScanTimes();
        var background = unidentified.Draw(fitted, scanTimes[0], scanTimes[^1], seed: 1);

        var charges = models.Concat(background).Select(m => (GaussianChargeDistribution)m.ChargeDistribution).ToArray();
        int minZ = Math.Max(2, (int)Math.Floor(charges.Min(c => c.MuZ - 3 * c.SigmaZ)));
        int maxZ = Math.Min(80, (int)Math.Ceiling(charges.Max(c => c.MuZ + 3 * c.SigmaZ)));
        Console.WriteLine($"species {species.Length}, unidentified {background.Length} ({unidentified}), charges {minZ}-{maxZ}, {scanTimes.Length} scans");

        var template = new NoiseFloorModel();
        var export = new Simulator().WriteMzml(
            models, minZ, maxZ, width, scanTimes, IdOnlyMzml("notemplate", Charge),
            background: background, scanNoiseFromSignal: scans => profile.NoiseModels(scans, template));

        Console.WriteLine($"wrote {export.MzmlPath}: {export.ScanCount} scans, {export.PeakCount / (double)export.ScanCount:F0} peaks/scan, " +
                          $"{export.FeatureCount} features, total {sw.Elapsed}");
    }

    /// <summary>
    /// Scores the sparse scans of the busy region: the real MS1 scans from 30 to 50 min with the
    /// fewest peaks at S/N ≥ 10 (but at least 20), where a few species carry clear signal. For each,
    /// the nearest scan of every simulation is compared on its S/N ≥ 10 peaks (threshold from the real
    /// scan's injection time), and the identified species eluting there are listed.
    /// </summary>
    [Test]
    [Explicit("Scores the sparsest real scans of 30-50 min against each simulation")]
    public static void CompareSparseScans()
    {
        var real = ReadMs1(RawPath(TestStem));
        var simulations = new List<(string Label, MsDataScan[] Scans)> { ("fitted", ReadMs1(FittedMzml(TestStem))) };
        foreach (string template in new[] { "train", "notemplate" })
            if (File.Exists(IdOnlyMzml(template, Charge)))
                simulations.Add(($"ids, {template}", ReadMs1(IdOnlyMzml(template, Charge))));
        var (fitted, _, _, _) = Simulator.ReadGroundTruth(FittedSidecar(TestStem));

        NoiseFloorModel NoiseOf(MsDataScan s) => new(NoiseFloorModel.JurkatNoiseTimesInjectionTime / (s.InjectionTime ?? 50));
        (double[] Mz, double[] I) Bright(MsDataScan s, NoiseFloorModel n) =>
            SignalPeaks.Select(s.MassSpectrum.XArray, s.MassSpectrum.YArray, n, 10, 600, 2000);

        // Sparse but identified: at least 30 % of the real bright peaks are matched by the fitted
        // simulation, so the signal there is species the simulation knows about. One per minute.
        var fittedScans = simulations[0].Scans;
        var sparse = real.Where(s => s.RetentionTime is >= 30 and <= 50)
            .Select(s =>
            {
                var n = NoiseOf(s);
                var (rMz, rI) = Bright(s, n);
                var (fMz, fI) = Bright(Nearest(fittedScans, s.RetentionTime), n);
                return (Scan: s, Count: rMz.Length, Explained: rMz.Length > 0 ? CentroidSpectrumComparison.Compare(rMz, rI, fMz, fI).RealMatchedFraction : 0);
            })
            .Where(p => p.Count >= 20 && p.Explained >= 0.3)
            .GroupBy(p => (int)p.Scan.RetentionTime)
            .Select(g => g.OrderBy(p => p.Count).First())
            .OrderBy(p => p.Count)
            .Take(8)
            .OrderBy(p => p.Scan.RetentionTime)
            .Select(p => (p.Scan, p.Count))
            .ToArray();

        foreach (var (r, count) in sparse)
        {
            var noise = NoiseOf(r);
            var (rMz, rI) = Bright(r, noise);
            var eluting = fitted.Where(f => f.Abundance > 0 && Math.Abs(f.RtProfile.Mu - r.RetentionTime) < 2 * f.RtProfile.Sigma + 0.1)
                .OrderByDescending(f => f.Abundance * f.RtProfile.Evaluate(r.RetentionTime)).Take(4)
                .Select(f => $"{f.MonoisotopicMass:F0} Da").ToArray();
            Console.WriteLine($"RT {r.RetentionTime:F2} (real scan {r.OneBasedScanNumber}, IT {r.InjectionTime:F1} ms): {count} peaks at S/N >= 10, " +
                              $"brightest fitted species eluting: {(eluting.Length == 0 ? "none" : string.Join(", ", eluting))}");
            foreach (var (label, scans) in simulations)
            {
                var s = Nearest(scans, r.RetentionTime);
                var (sMz, sI) = Bright(s, noise);
                var c = CentroidSpectrumComparison.Compare(rMz, rI, sMz, sI);
                var whole = CentroidSpectrumComparison.Compare(r.MassSpectrum.XArray, r.MassSpectrum.YArray, s.MassSpectrum.XArray, s.MassSpectrum.YArray, 600, 2000);
                Console.WriteLine($"    {label,-16} scan {s.OneBasedScanNumber,5}: bright {c.SimulatedPeaks,4} vs {c.RealPeaks,4}, real matched {c.RealMatchedFraction,4:P0}, " +
                                  $"sim matched {c.SimulatedMatchedFraction,4:P0}, bright cos {c.Cosine:F2}, bright TIC {c.TicRatio:F2}x, whole-scan cos {whole.Cosine:F2}");
            }
        }
    }

    /// <summary>
    /// Splits the no-template simulation's ion load by source on every 10th held-out MS1 scan from 30
    /// to 50 min: the rendered signal (identified, unidentified), the injection time AGC gives it, and
    /// the noise that follows, against the real scan's TIC, bright TIC and injection time.
    /// </summary>
    [Test]
    [Explicit("Diagnoses the no-template simulation's ion load against the held-out run")]
    public static void DiagnoseTemplateFreeLoad()
    {
        var real = ReadMs1(RawPath(TestStem)).Where(s => s.RetentionTime is >= 30 and <= 50).Where((_, i) => i % 10 == 0).ToArray();
        var priors = IdOnlyPriors.Load(PriorsPath);
        var profile = AcquisitionProfile.Load(AcquisitionProfilePath);
        var unidentified = System.Text.Json.JsonSerializer.Deserialize<UnidentifiedAnalytes>(File.ReadAllText(UnidentifiedPath))!;
        var species = SpeciesGrouper.Group(new MmResultLoader().LoadQualified(PsmPath(TestStem), TestStem, QValue));
        var models = priors.Predict(species, chargePredictor: Charge);
        var (fitted, _, _, width) = Simulator.ReadGroundTruth(FittedSidecar(TrainStem));
        var background = unidentified.Draw(fitted, profile.FirstScanTime, profile.LastScanTime, seed: 1);
        double[] times = real.Select(s => s.RetentionTime).ToArray();

        var simulator = new Simulator();
        var identifiedScans = simulator.SimulateCentroided(models, 2, 47, width, times).Scans;
        var backgroundScans = simulator.SimulateCentroided(background, 2, 47, width, times).Scans;
        var allScans = simulator.SimulateCentroided(models.Concat(background).ToArray(), 2, 47, width, times).Scans;
        var noise = profile.NoiseModels(allScans);
        var rng = new Random(1);

        Console.WriteLine($"  {"RT",5} {"real TIC",9} {"real >10",9} {"real IT",7} | {"ident",9} {"unident",9} {"sim IT",7} {"noise TIC",9} {"sim >10",9}");
        for (int i = 0; i < real.Length; i++)
        {
            var realNoise = new NoiseFloorModel(NoiseFloorModel.JurkatNoiseTimesInjectionTime / (real[i].InjectionTime ?? 50));
            double realBright = SignalPeaks.Select(real[i].MassSpectrum.XArray, real[i].MassSpectrum.YArray, realNoise, 10).Intensity.Sum();
            double simBright = SignalPeaks.Select(allScans[i].MassSpectrum.XArray, allScans[i].MassSpectrum.YArray, realNoise, 10).Intensity.Sum();
            var peaks = new List<SimulatedPeak>();
            noise[i].SampleScan(rng, peaks);
            double simIt = NoiseFloorModel.JurkatNoiseTimesInjectionTime / noise[i].NoiseLevelAtReferenceMz;
            Console.WriteLine($"  {real[i].RetentionTime,5:F1} {real[i].MassSpectrum.SumOfAllY,9:G3} {realBright,9:G3} {real[i].InjectionTime,7:F1} | " +
                              $"{identifiedScans[i].MassSpectrum.SumOfAllY,9:G3} {backgroundScans[i].MassSpectrum.SumOfAllY,9:G3} {simIt,7:F1} {peaks.Sum(p => p.Intensity),9:G3} {simBright,9:G3}");
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
