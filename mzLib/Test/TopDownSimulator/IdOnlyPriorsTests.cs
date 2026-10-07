using System;
using System.IO;
using System.Linq;
using NUnit.Framework;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;
using TopDownSimulator.Prediction;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

[TestFixture]
public class IdOnlyPriorsTests
{
    private static MmResultRecord Record(
        string id, double mass, double rt, int charge, double score = 10, double? intensity = 1e6,
        string sequence = "PEPTIDEK") =>
        new("run", 1, 2, charge, mass, rt, score, sequence, "P1", id, intensity);

    [Test]
    public void GrouperMergesChargeVariantsAndIsotopeErrorsButNotDistinctSpecies()
    {
        var records = new[]
        {
            Record("a", 10000.00, 20.0, 9, score: 30),
            Record("b", 10000.01, 20.2, 12, score: 20),        // same species, other charge
            Record("c", 10001.01, 20.1, 10, score: 10),        // off by one isotope
            Record("d", 10000.50, 20.0, 9, score: 25),         // half a dalton away: distinct
            Record("e", 10000.00, 25.0, 9, score: 5),          // elutes far later: distinct
        };

        var species = SpeciesGrouper.Group(records);

        Assert.That(species, Has.Length.EqualTo(3));
        var first = species.Single(s => s.Anchor.Identifier == "a");
        Assert.That(first.Members.Select(m => m.Identifier), Is.EquivalentTo(new[] { "a", "b", "c" }));
        Assert.That(first.MinCharge, Is.EqualTo(7));
        Assert.That(first.MaxCharge, Is.EqualTo(14));
    }

    [Test]
    public void PrecursorGroupingMergesInterpretationsOfOneEnvelope()
    {
        // Three interpretations of the envelope measured at ~13766, one real neighbour +29 Da away.
        var records = new[]
        {
            Record("a", 13766.50, 41.0, 18, score: 30) with { PrecursorMass = 13766.5 },
            Record("b", 13781.52, 41.1, 17, score: 20) with { PrecursorMass = 13765.5 },   // claims +15 Da, measured here
            Record("c", 13748.52, 41.2, 19, score: 15) with { PrecursorMass = 13768.5 },   // claims -18 Da, measured here
            Record("d", 13795.52, 41.0, 18, score: 25) with { PrecursorMass = 13795.5 },   // a separate envelope
        };

        var byTheory = SpeciesGrouper.Group(records);
        var byEnvelope = SpeciesGrouper.Group(records, precursorMassWindow: SpeciesGrouper.DefaultPrecursorMassWindow);

        Assert.That(byTheory, Has.Length.EqualTo(4));
        Assert.That(byEnvelope, Has.Length.EqualTo(2));
        Assert.That(byEnvelope.Single(s => s.Anchor.Identifier == "a").Members.Select(m => m.Identifier),
            Is.EquivalentTo(new[] { "a", "b", "c" }));
    }

    [Test]
    public void FitRecoversTheRegressionsItWasTrainedOn()
    {
        // Abundance = 10 x summed intensity, charge centre = mean charge + 0.1 K+R, apex 0.05 min late.
        var rng = new Random(1);
        var species = Enumerable.Range(0, 60).Select(i =>
        {
            int charge = 8 + i % 10;
            var r = Record($"id{i}", 8000 + 200 * i, 10 + i * 0.5, charge,
                intensity: Math.Pow(10, 5 + rng.NextDouble() * 2), sequence: new string('K', i % 7) + "PEPTIDE");
            return new IdentifiedSpecies(r, new[] { r }, charge - 2, charge + 2);
        }).ToArray();

        var fitted = species.Select(s => new ProteoformModel(
            s.Anchor.MonoisotopicMass,
            10 * s.Anchor.PrecursorIntensity!.Value,
            new EmgProfile(s.Anchor.RetentionTime + 0.05, 0.1, 0.02),
            new GaussianChargeDistribution(s.Anchor.PrecursorCharge + 0.1 * IdOnlyPriors.BasicResidueCount(s), 1.3),
            s.Anchor.Identifier)).ToArray();

        var priors = IdOnlyPriors.Fit(species, fitted, "synthetic");

        Assert.That(priors.RtApexOffset, Is.EqualTo(0.05).Within(1e-9));
        Assert.That(priors.LogAbundanceSlope, Is.EqualTo(1).Within(1e-9));
        Assert.That(priors.LogAbundanceIntercept, Is.EqualTo(1).Within(1e-9));
        Assert.That(priors.ChargeMuChargeSlope, Is.EqualTo(1).Within(1e-9));
        Assert.That(priors.ChargeMuBasicSlope, Is.EqualTo(0.1).Within(1e-9));

        var predicted = priors.Predict(species);
        for (int i = 0; i < species.Length; i++)
        {
            Assert.That(predicted[i].Abundance, Is.EqualTo(fitted[i].Abundance).Within(1e-6).Percent);
            Assert.That(predicted[i].RtProfile.Mu, Is.EqualTo(fitted[i].RtProfile.Mu).Within(1e-9));
        }
    }

    [Test]
    public void ResidueModelRecoversAWeightPerResidue()
    {
        // mu_z = 1 + 0.05 sqrt(mass) + 0.3 K + 0.2 R + 0.1 H - 0.15 D - 0.05 E, exactly.
        var rng = new Random(3);
        var species = Enumerable.Range(0, 80).Select(i =>
        {
            string seq = new string('K', rng.Next(12)) + new string('R', rng.Next(8)) + new string('H', rng.Next(4))
                         + new string('D', rng.Next(10)) + new string('E', rng.Next(10)) + "PEPTIDE";
            var r = Record($"id{i}", 5000 + 300 * i, 10, 10, sequence: seq);
            return new IdentifiedSpecies(r, new[] { r }, 8, 12);
        }).ToArray();
        double Truth(IdentifiedSpecies s) => 1 + 0.05 * Math.Sqrt(s.Anchor.MonoisotopicMass)
            + 0.3 * s.Anchor.FullSequence.Count(c => c == 'K') + 0.2 * s.Anchor.FullSequence.Count(c => c == 'R')
            + 0.1 * s.Anchor.FullSequence.Count(c => c == 'H') - 0.15 * s.Anchor.FullSequence.Count(c => c == 'D')
            - 0.05 * (s.Anchor.FullSequence.Count(c => c == 'E') - 2);   // "PEPTIDE" carries two E and one D
        var targets = species.Select(Truth).ToArray();

        var model = ResidueChargeModel.Fit(species, targets);

        Assert.That(model.Residues, Is.EqualTo("KRHDE"));
        Assert.That(model.SqrtMassSlope, Is.EqualTo(0.05).Within(1e-8));
        Assert.That(model.Weights, Is.EqualTo(new[] { 0.3, 0.2, 0.1, -0.15, -0.05 }).Within(1e-8));
        Assert.That(model.ResidualSd, Is.LessThan(1e-8));
        Assert.That(model.Predict(species[7]), Is.EqualTo(targets[7]).Within(1e-8));
    }

    [Test]
    public void SequenceChargeRangeCoversThePredictedEnvelope()
    {
        var shapes = new[] { new ShapeSample(0.1, 0, 2.0) };
        var priors = new IdOnlyPriors(0, 0, 1, 5, 0, 1, 0, shapes, "test",
            SequenceChargeIntercept: 10, SequenceChargeSqrtMassSlope: 0, SequenceChargeNetBasicSlope: 0);
        var r = Record("x", 9000, 10, 9);
        var range = priors.SequenceChargeRange(new IdentifiedSpecies(r, new[] { r }, 7, 11));
        Assert.That(range, Is.EqualTo((4, 16)));
    }

    [Test]
    public void PredictionIsStablePerSpeciesWhateverElseIsInTheSet()
    {
        var shapes = Enumerable.Range(0, 50).Select(i => new ShapeSample(0.05 + 0.002 * i, 0, 1 + 0.02 * i)).ToArray();
        var priors = new IdOnlyPriors(0, 0, 1, 5, 0, 1, 0, shapes, "test");
        var a = new IdentifiedSpecies(Record("a", 9000, 10, 9), new[] { Record("a", 9000, 10, 9) }, 7, 11);
        var b = new IdentifiedSpecies(Record("b", 9500, 12, 10), new[] { Record("b", 9500, 12, 10) }, 8, 12);

        var alone = priors.Predict(new[] { a })[0];
        var together = priors.Predict(new[] { b, a })[1];

        Assert.That(together, Is.EqualTo(alone));
    }

    [Test]
    public void MissingPrecursorIntensityFallsBackToThePopulationMedian()
    {
        var priors = new IdOnlyPriors(0, 0, 1, 5.5, 0, 1, 0, new[] { new ShapeSample(0.1, 0, 1) }, "test");
        var r = Record("x", 9000, 10, 9, intensity: null);
        var model = priors.Predict(new[] { new IdentifiedSpecies(r, new[] { r }, 7, 11) })[0];
        Assert.That(model.Abundance, Is.EqualTo(Math.Pow(10, 5.5)).Within(1e-6).Percent);
    }

    [Test]
    public void BasicResidueCountIgnoresModificationText()
    {
        var r = Record("x", 9000, 10, 9, sequence: "PK[Common Biological:Acetylation on K]RAR[UniProt:Methyl]K");
        Assert.That(IdOnlyPriors.BasicResidueCount(new IdentifiedSpecies(r, new[] { r }, 7, 11)), Is.EqualTo(4));
    }

    [Test]
    public void PriorsRoundTripThroughJson()
    {
        var priors = new IdOnlyPriors(-0.01, 0.1, 0.8, 5.1, 1.7, 0.83, 0.06,
            new[] { new ShapeSample(0.11, 0.03, 1.4), new ShapeSample(0.09, 0, 2.1) }, "unit test");
        string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, "priors-roundtrip.json");
        try
        {
            priors.Save(path);
            var loaded = IdOnlyPriors.Load(path);
            Assert.That(loaded with { Shapes = priors.Shapes }, Is.EqualTo(priors));
            Assert.That(loaded.Shapes, Is.EqualTo(priors.Shapes));
        }
        finally
        {
            File.Delete(path);
        }
    }

    [Test]
    public void GroundTruthSidecarRoundTrips()
    {
        var models = new[]
        {
            new ProteoformModel(10000, 1.5e6, new EmgProfile(20, 0.2, 0.05), new GaussianChargeDistribution(8.3, 1.1), "A"),
            new ProteoformModel(12500, 4e5, new EmgProfile(21, 0.15, 0), new GaussianChargeDistribution(10.2, 1.4), "B"),
        };
        var width = OrbitrapPeakWidth.FromSigmaAt(1000, 0.01);
        string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, "sidecar-roundtrip.tsv");
        try
        {
            Simulator.WriteGroundTruth(models, 6, 14, width, path);
            var (read, minZ, maxZ, readWidth) = Simulator.ReadGroundTruth(path);

            Assert.That(read, Is.EqualTo(models));
            Assert.That((minZ, maxZ), Is.EqualTo((6, 14)));
            Assert.That(readWidth, Is.EqualTo(width));
        }
        finally
        {
            File.Delete(path);
        }
    }
}
