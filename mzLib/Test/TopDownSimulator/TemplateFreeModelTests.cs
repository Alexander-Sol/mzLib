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

    [Test]
    public void UnidentifiedLoadScalesAndFollowsTheIdentifiedLoad()
    {
        var fitted = new[]
        {
            new ProteoformModel(10000, 1e6, new EmgProfile(30, 0.1, 0), new GaussianChargeDistribution(12, 3), "a"),
        };
        var identified = new[] { (Rt: 50.0, Intensity: 9.0), (Rt: 60.0, Intensity: 1.0), (Rt: 70.0, Intensity: 0.0) };
        var plain = new UnidentifiedAnalytes(2000, 0, LogAbundanceSd: 0, RtSd: 0.1);

        // Without following, the identified load changes nothing but the abundance.
        var scaled = plain.Draw(fitted, 0, 90, seed: 5, identified, loadScale: 4);
        var unscaled = plain.Draw(fitted, 0, 90, seed: 5);
        Assert.That(scaled.Select(m => m.Abundance / 4), Is.EqualTo(unscaled.Select(m => m.Abundance)).Within(1e-6));
        Assert.That(scaled.Select(m => m.RtProfile.Mu), Is.EqualTo(unscaled.Select(m => m.RtProfile.Mu)));

        // Following entirely, elution centres land on the identifications in proportion to intensity.
        var followed = (plain with { FollowIdentified = 1 }).Draw(fitted, 0, 90, seed: 5, identified);
        int at50 = followed.Count(m => Math.Abs(m.RtProfile.Mu - 50) < 1);
        int at60 = followed.Count(m => Math.Abs(m.RtProfile.Mu - 60) < 1);
        Assert.That(at50 + at60, Is.EqualTo(2000));
        Assert.That(at50 / 2000.0, Is.EqualTo(0.9).Within(0.03));

        // Half following: the rest stay on the training model's elution time.
        var half = (plain with { FollowIdentified = 0.5 }).Draw(fitted, 0, 90, seed: 5, identified);
        Assert.That(half.Count(m => Math.Abs(m.RtProfile.Mu - 30) < 1) / 2000.0, Is.EqualTo(0.5).Within(0.04));
    }

    [Test]
    public void LoadScaleNeedsEnoughIdentifications()
    {
        int n = UnidentifiedAnalytes.MinimumIdentificationsForLoad;
        Assert.That(UnidentifiedAnalytes.LoadScale(Enumerable.Repeat(3.0, n), Enumerable.Repeat(1.0, 2 * n)), Is.EqualTo(1.5).Within(1e-12));
        Assert.That(UnidentifiedAnalytes.LoadScale(Enumerable.Repeat(3.0, n - 1), Enumerable.Repeat(1.0, 2 * n)), Is.Null);
        // Missing or zero intensities do not count as identifications.
        var withGaps = Enumerable.Repeat(3.0, n - 1).Append(0).Append(double.NaN);
        Assert.That(UnidentifiedAnalytes.LoadScale(withGaps, Enumerable.Repeat(1.0, 2 * n)), Is.Null);
    }
}
