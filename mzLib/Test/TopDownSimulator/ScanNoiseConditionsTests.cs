using System;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using TopDownSimulator.Model;
using TopDownSimulator.Noise;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

[TestFixture]
public class ScanNoiseConditionsTests
{
    private static readonly IPeakWidthModel Width = new ConstantPeakWidth(0.012);

    private static MsDataScan Scan(double? injectionTime, double[] mz, double[] intensity, int oneBased = 1) =>
        new(new MzSpectrum(mz, intensity, false), oneBased, 1, true, Polarity.Positive, 30.0,
            new MzRange(600, 2000), "synthetic", MZAnalyzerType.Orbitrap, intensity.Sum(),
            injectionTime, null, $"scan={oneBased}");

    [Test]
    public void AmplitudeIsTheConstantOverInjectionTime()
    {
        var template = new NoiseFloorModel(noiseLevelAtReferenceMz: 500);
        var scans = new[]
        {
            Scan(0.5, Array.Empty<double>(), Array.Empty<double>()),
            Scan(40.0, Array.Empty<double>(), Array.Empty<double>()),
            Scan(null, Array.Empty<double>(), Array.Empty<double>()),
        };

        var models = ScanNoiseConditions.FromSourceScans(scans, template, noiseTimesInjectionTime: 2e4);

        Assert.That(models[0].NoiseLevelAtReferenceMz, Is.EqualTo(4e4).Within(1e-9));
        Assert.That(models[1].NoiseLevelAtReferenceMz, Is.EqualTo(500).Within(1e-9));
        Assert.That(models[2].NoiseLevelAtReferenceMz, Is.EqualTo(500), "No injection time keeps the template's level.");
    }

    [Test]
    public void DensityCountsOnlyPeaksUnderTheSignalToNoiseCutoffInTheirOwnBin()
    {
        var template = new NoiseFloorModel(noiseLevelAtReferenceMz: 1000, levelDispersion: 0);
        double level650 = template.NoiseLevelAt(650);

        // Two quiet peaks and one bright one in 600-700, one quiet peak in 900-1000, one outside.
        var scan = Scan(20.0,
            new[] { 610.0, 650.0, 660.0, 950.0, 2500.0 },
            new[] { 2 * level650, 3 * level650, 50 * level650, 2 * template.NoiseLevelAt(950), 1.0 });

        var model = ScanNoiseConditions.FromSourceScans(new[] { scan }, template, noiseTimesInjectionTime: 2e4)[0];

        Assert.That(model.PeaksPerScanInBin(0), Is.EqualTo(2));
        Assert.That(model.PeaksPerScanInBin(3), Is.EqualTo(1));
        Assert.That(model.ExpectedPeaksPerScan, Is.EqualTo(3));
    }

    [Test]
    public void AmplitudeOnlyConditioningKeepsTheTemplateDensity()
    {
        var template = new NoiseFloorModel();
        var scan = Scan(2.0, new[] { 650.0 }, new[] { 1.0 });

        var model = ScanNoiseConditions.FromSourceScans(new[] { scan }, template, conditionDensity: false)[0];

        Assert.That(model.ExpectedPeaksPerScan, Is.EqualTo(template.ExpectedPeaksPerScan).Within(1e-9));
        Assert.That(model.NoiseLevelAtReferenceMz, Is.EqualTo(NoiseFloorModel.JurkatNoiseTimesInjectionTime / 2.0));
    }

    [Test]
    public void InjectorUsesEachScansOwnModel()
    {
        var empty = new double[NoiseFloorModel.BinCount];
        var dense = Enumerable.Repeat(200.0, NoiseFloorModel.BinCount).ToArray();
        var quiet = new NoiseFloorModel(1000, levelDispersion: 0, peaksPerScanPerBin: empty);
        var loud = new NoiseFloorModel(1e5, levelDispersion: 0, peaksPerScanPerBin: dense);

        var scans = new[]
        {
            Scan(1, Array.Empty<double>(), Array.Empty<double>(), 1),
            Scan(1, Array.Empty<double>(), Array.Empty<double>(), 2),
        };

        var (result, summary) = new NoiseInjector(new[] { quiet, loud }, Width, seed: 3).Apply(scans);

        Assert.That(result[0].MassSpectrum.XArray, Is.Empty, "An empty density table must add no noise.");
        Assert.That(result[1].MassSpectrum.XArray.Length, Is.EqualTo(200 * NoiseFloorModel.BinCount).Within(400));
        Assert.That(result[1].MassSpectrum.YArray.Min(), Is.GreaterThanOrEqualTo(loud.NoiseLevelAt(600) * 1.5));
        Assert.That(summary.NoisePeaks, Is.EqualTo(result[1].MassSpectrum.XArray.Length));
    }

    [Test]
    public void InjectorRejectsAScanCountItWasNotBuiltFor()
    {
        var injector = new NoiseInjector(new[] { new NoiseFloorModel() }, Width);
        var scans = new[] { Scan(1, new[] { 700.0 }, new[] { 1.0 }, 1), Scan(1, new[] { 700.0 }, new[] { 1.0 }, 2) };
        Assert.Throws<ArgumentException>(() => injector.Apply(scans));
    }

    [Test]
    public void PerScanFloorsDecideWhichSignalSurvivesInEachScan()
    {
        // The same apex-height signal in two scans, against a quiet floor and one 1000x higher.
        var model = new ProteoformModel(10000, 1e6, new EmgProfile(30.0, 0.5, 0), new GaussianChargeDistribution(9, 0.8));
        var none = new double[NoiseFloorModel.BinCount];
        var quiet = new NoiseFloorModel(10, levelDispersion: 0, peaksPerScanPerBin: none);
        var loud = new NoiseFloorModel(1e4, levelDispersion: 0, peaksPerScanPerBin: none);

        var prepared = new Simulator().PrepareMs1(
            new[] { model }, 7, 11, Width, new[] { 30.0, 30.0 },
            jitter: PeakJitterModel.None, scanNoise: new[] { quiet, loud });

        int quietPeaks = prepared.Ms1Scans[0].MassSpectrum.XArray.Length;
        int loudPeaks = prepared.Ms1Scans[1].MassSpectrum.XArray.Length;
        Assert.That(quietPeaks, Is.GreaterThan(loudPeaks));
        Assert.That(prepared.Floors[0].At(1000), Is.LessThan(prepared.Floors[1].At(1000)));

        foreach (var (scan, floor) in prepared.Ms1Scans.Zip(prepared.Floors))
            for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
                Assert.That(scan.MassSpectrum.YArray[i], Is.GreaterThanOrEqualTo(floor.At(scan.MassSpectrum.XArray[i])));
    }
}
