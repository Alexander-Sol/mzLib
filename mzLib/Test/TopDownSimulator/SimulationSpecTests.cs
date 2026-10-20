using System;
using System.IO;
using System.Linq;
using NUnit.Framework;
using TopDownSimulator.Model;
using TopDownSimulator.Noise;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

/// <summary>
/// SimSpec v1 round-trip and guard tests. No external data: everything runs against temp files.
/// The contract these pin down lives in ProteoStar's agent_info/TopDown-SimSpec-Contract-v1.md.
/// </summary>
[TestFixture]
public class SimulationSpecTests
{
    private static ProteoformModel[] TwoModels() => new[]
    {
        new ProteoformModel(5000.25, 1.0e9,
            new EmgProfile(20.0, 0.22, 0.08),
            new GaussianChargeDistribution(7.5, 1.2), "P1|SMALL"),
        new ProteoformModel(35000.75, 4.2e8,
            new EmgProfile(24.5, 0.30, 0.12),
            new GaussianChargeDistribution(30.0, 3.5), "P2|LARGE"),
    };

    private static double[] ScanTimes(int n = 12, double start = 19.5, double step = 0.1) =>
        Enumerable.Range(0, n).Select(i => start + i * step).ToArray();

    private string _tempPath;

    [SetUp]
    public void SetUp() =>
        _tempPath = Path.Combine(TestContext.CurrentContext.TestDirectory, $"simspec-{Guid.NewGuid():N}.json");

    [TearDown]
    public void TearDown()
    {
        if (File.Exists(_tempPath)) File.Delete(_tempPath);
    }

    [Test]
    public void OrbitrapSpecRoundTripsThroughDisk()
    {
        var models = TwoModels();
        double[] scanTimes = ScanTimes();
        var reduction = new ScanReductionOptions { RelativeIntensityThreshold = 1e-4, MinimumIntensity = 2.5 };

        var written = SimSpec.FromModels(
            models, minCharge: 3, maxCharge: 40, new OrbitrapPeakWidth(1.2e-7), scanTimes, reduction,
            noise: SimNoise.OrbitrapDefaults(7),
            generator: SimGenerator.Now("deadbeef"),
            sourceRun: new SimSourceRun("Jurkat_rep2.raw", "unit test"));
        written.Write(_tempPath);

        var read = SimSpec.Read(_tempPath);

        Assert.That(read.SchemaVersion, Is.EqualTo(SimSpec.CurrentSchemaVersion));
        Assert.That(read.ChargeRange.Min, Is.EqualTo(3));
        Assert.That(read.ChargeRange.Max, Is.EqualTo(40));
        Assert.That(read.PeakWidthModel.Kind, Is.EqualTo("Orbitrap"));
        Assert.That(read.PeakWidthModel.K, Is.EqualTo(1.2e-7).Within(1e-18));
        Assert.That(read.ScanTimes, Is.EqualTo(scanTimes).AsCollection);
        Assert.That(read.Reduction!.RelativeIntensityThreshold, Is.EqualTo(1e-4));
        Assert.That(read.Reduction!.AbsoluteIntensityFloor, Is.EqualTo(2.5));
        Assert.That(read.Noise!.Seed, Is.EqualTo(7));
        Assert.That(read.Noise!.Jitter.Enabled, Is.True);
        Assert.That(read.Noise!.Floor.DensityScale, Is.EqualTo(1.0));
        Assert.That(read.Noise!.Floor.SignalToNoise.Threshold, Is.EqualTo(1.94));
        Assert.That(read.Proteoforms, Has.Count.EqualTo(2));

        var p2 = read.Proteoforms[1];
        Assert.That(p2.Identifier, Is.EqualTo("P2|LARGE"));
        Assert.That(p2.MonoisotopicMass, Is.EqualTo(35000.75));
        Assert.That(p2.Abundance, Is.EqualTo(4.2e8));
        Assert.That(p2.Rt.Tau, Is.EqualTo(0.12));
        Assert.That(p2.Charge.MuZ, Is.EqualTo(30.0));
        Assert.That(p2.Charge.SigmaZ, Is.EqualTo(3.5));
    }

    [Test]
    public void ToSimulationInputsReconstructsTheModels()
    {
        var models = TwoModels();
        var spec = SimSpec.FromModels(models, 3, 40, new ConstantPeakWidth(0.01), ScanTimes());

        var inputs = spec.ToSimulationInputs();

        Assert.That(inputs.WidthModel, Is.TypeOf<ConstantPeakWidth>());
        Assert.That(((ConstantPeakWidth)inputs.WidthModel).SigmaMz, Is.EqualTo(0.01));
        Assert.That(inputs.Proteoforms, Has.Length.EqualTo(2));

        var original = models[0];
        var rebuilt = inputs.Proteoforms[0];
        Assert.That(rebuilt.MonoisotopicMass, Is.EqualTo(original.MonoisotopicMass));
        Assert.That(rebuilt.Abundance, Is.EqualTo(original.Abundance));
        Assert.That(rebuilt.RtProfile, Is.EqualTo(original.RtProfile));           // record equality
        Assert.That(rebuilt.ChargeDistribution, Is.EqualTo(original.ChargeDistribution));
        Assert.That(rebuilt.Identifier, Is.EqualTo(original.Identifier));
    }

    [Test]
    public void ReadRejectsAnUnknownSchemaVersion()
    {
        File.WriteAllText(_tempPath,
            "{\"schemaVersion\":999,\"instrument\":{\"polarity\":\"Positive\",\"analyzer\":\"Orbitrap\",\"msLevel\":1}," +
            "\"chargeRange\":{\"min\":3,\"max\":40},\"peakWidthModel\":{\"kind\":\"Constant\",\"sigma\":0.01}," +
            "\"scanTimes\":[1.0],\"proteoforms\":[]}");

        Assert.That(() => SimSpec.Read(_tempPath), Throws.InstanceOf<InvalidDataException>());
    }

    [Test]
    public void NoiseConfigRoundTripsAndReconstructsTheModels()
    {
        var floor = new NoiseFloorModel(
            noiseLevelAtReferenceMz: 1500, densityScale: 0.5,
            signalToNoise: new SignalToNoiseDistribution(2.0, 1.2), levelDispersion: 0.3);
        var noise = SimNoise.FromModels(floor, PeakJitterModel.Orbitrap, seed: 42);

        var spec = SimSpec.FromModels(
            TwoModels(), 3, 40, new OrbitrapPeakWidth(1.2e-7), ScanTimes(), noise: noise);
        spec.Write(_tempPath);

        var read = SimSpec.Read(_tempPath).Noise!;
        Assert.That(read.Seed, Is.EqualTo(42));

        var rebuiltFloor = read.ToFloorModel();
        Assert.That(rebuiltFloor.NoiseLevelAtReferenceMz, Is.EqualTo(1500));
        Assert.That(rebuiltFloor.DensityScale, Is.EqualTo(0.5));
        Assert.That(rebuiltFloor.LevelDispersion, Is.EqualTo(0.3));
        Assert.That(rebuiltFloor.SignalToNoise.Threshold, Is.EqualTo(2.0));
        Assert.That(rebuiltFloor.SignalToNoise.MeanExcess, Is.EqualTo(1.2));

        var rebuiltJitter = read.ToJitterModel();
        Assert.That(rebuiltJitter, Is.Not.Null);
        Assert.That(rebuiltJitter!.MzPpmAtUnitSignalToNoise,
            Is.EqualTo(PeakJitterModel.Orbitrap.MzPpmAtUnitSignalToNoise));
        Assert.That(rebuiltJitter.CommonModeMzPpm, Is.EqualTo(PeakJitterModel.Orbitrap.CommonModeMzPpm));
        Assert.That(rebuiltJitter.DropPeaksFallingBelowFloor, Is.True);
    }

    [Test]
    public void DisabledJitterReconstructsAsNull()
    {
        var noise = SimNoise.FromModels(new NoiseFloorModel(), PeakJitterModel.None, seed: 0);

        Assert.That(noise.Jitter.Enabled, Is.False);
        Assert.That(noise.ToJitterModel(), Is.Null);
    }

    [Test]
    public void NonGaussianChargeDistributionIsRejectedUpFront()
    {
        var model = new ProteoformModel(5000.0, 1e9,
            new EmgProfile(20.0, 0.2, 0.08), new DummyChargeDistribution());

        Assert.That(
            () => SimSpec.FromModels(new[] { model }, 3, 40, new ConstantPeakWidth(0.01), ScanTimes()),
            Throws.InstanceOf<NotSupportedException>());
    }

    private sealed class DummyChargeDistribution : IChargeStateDistribution
    {
        public double Evaluate(int z) => 1.0;
    }
}
