using System;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using TopDownSimulator.Model;
using TopDownSimulator.Noise;
using TopDownSimulator.Prediction;

namespace Test.TopDownSimulator;

[TestFixture]
public class TemplateFreeModelTests
{
    private static MsDataScan Scan(int number, double rt, double tic, double injectionTime)
    {
        // Two peaks carrying the TIC, plus enough low ones in the wash-like region to have a density.
        var mz = new[] { 700.0, 900.0 };
        var intensity = new[] { tic / 2, tic / 2 };
        return new MsDataScan(new MzSpectrum(mz, intensity, false), number, 1, true, Polarity.Positive, rt,
            new MzRange(600, 2000), "", MZAnalyzerType.Orbitrap, tic, injectionTime, null, $"scan={number}");
    }

    [Test]
    public void AgcRecoversTheTargetAndTheCap()
    {
        // IT = min(50, 1e9 / TIC) exactly.
        var scans = Enumerable.Range(0, 40).Select(i =>
        {
            double tic = Math.Pow(10, 6 + i * 0.1);
            return Scan(i + 1, i * 0.03, tic, Math.Min(50, 1e9 / tic));
        }).ToArray();

        var agc = AutomaticGainControl.Fit(scans);

        Assert.That(agc.MaxInjectionTime, Is.EqualTo(50));
        Assert.That(agc.ChargeTarget, Is.EqualTo(1e9).Within(1).Percent);
        Assert.That(agc.InjectionTime(1e6), Is.EqualTo(50));
        Assert.That(agc.InjectionTime(1e10), Is.EqualTo(0.1).Within(2).Percent);
    }

    [Test]
    public void ProfileStepsScanTimesByEachBinsInterval()
    {
        // 0.02 min apart for the first two minutes, 0.05 min apart after.
        var times = Enumerable.Range(0, 100).Select(i => i * 0.02)
            .Concat(Enumerable.Range(1, 40).Select(i => 2 + i * 0.05)).ToArray();
        var scans = times.Select((t, i) => Scan(i + 1, t, 1e8, 10)).ToArray();

        var profile = AcquisitionProfile.Learn(scans, "test");
        var generated = profile.ScanTimes();

        Assert.That(profile.Ms1Interval[0], Is.EqualTo(0.02).Within(1e-9));
        Assert.That(profile.Ms1Interval[^1], Is.EqualTo(0.05).Within(1e-9));
        Assert.That(generated.Length, Is.EqualTo(times.Length).Within(2));
        Assert.That(generated[0], Is.EqualTo(0));
    }

    [Test]
    public void ProfileNoiseFollowsTheSignalLoadThroughAgc()
    {
        var times = Enumerable.Range(0, 60).Select(i => i * 0.05).ToArray();
        var scans = times.Select((t, i) => Scan(i + 1, t, i < 30 ? 1e7 : 1e10, i < 30 ? 50 : 0.1)).ToArray();
        var profile = AcquisitionProfile.Learn(scans, "test");

        var models = profile.NoiseModels(new[] { Scan(1, 0.5, 1e5, 0), Scan(2, 0.6, 1e10, 0) });

        // Quiet scan: capped at 50 ms. Loaded scan: IT = target / 1e10.
        Assert.That(models[0].NoiseLevelAtReferenceMz, Is.EqualTo(NoiseFloorModel.JurkatNoiseTimesInjectionTime / 50).Within(1e-6));
        Assert.That(models[1].NoiseLevelAtReferenceMz, Is.GreaterThan(10 * models[0].NoiseLevelAtReferenceMz));
    }

    [Test]
    public void UnidentifiedAnalytesAreShiftedCopiesOfTheFittedPopulation()
    {
        var fitted = new[]
        {
            new ProteoformModel(10000, 1e6, new EmgProfile(30, 0.1, 0), new GaussianChargeDistribution(12, 3), "a"),
            new ProteoformModel(20000, 1e7, new EmgProfile(40, 0.12, 0.05), new GaussianChargeDistribution(20, 4), "b"),
        };

        var drawn = new UnidentifiedAnalytes(500, LogAbundanceShift: -1, LogAbundanceSd: 0).Draw(fitted, 25, 45, seed: 3);

        Assert.That(drawn, Has.Length.EqualTo(500));
        Assert.That(drawn.All(m => m.Identifier!.StartsWith("unidentified:")));
        Assert.That(drawn.All(m => m.RtProfile.Mu is >= 25 and <= 45));
        Assert.That(drawn.All(m => m.Abundance == 1e5 || Math.Abs(m.Abundance - 1e6) < 1e-3));
        var again = new UnidentifiedAnalytes(500, -1, 0).Draw(fitted, 25, 45, seed: 3);
        Assert.That(again, Is.EqualTo(drawn));
    }
}
