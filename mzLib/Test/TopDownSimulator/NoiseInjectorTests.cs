using System;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using TopDownSimulator.Model;
using TopDownSimulator.Noise;

namespace Test.TopDownSimulator;

[TestFixture]
public class NoiseInjectorTests
{
    private static readonly IPeakWidthModel Width = OrbitrapPeakWidth.FromSigmaAt(1000, 0.01);

    private static MsDataScan BuildScan(int oneBased, double rt, double[] mz, double[] intensity) =>
        new(
            massSpectrum: new MzSpectrum(mz, intensity, false),
            oneBasedScanNumber: oneBased,
            msnOrder: 1,
            isCentroid: true,
            polarity: Polarity.Positive,
            retentionTime: rt,
            scanWindowRange: new MzRange(600, 2000),
            scanFilter: "synthetic",
            mzAnalyzer: MZAnalyzerType.Orbitrap,
            totalIonCurrent: intensity.Sum(),
            injectionTime: 1.0,
            noiseData: null,
            nativeId: $"scan={oneBased}");

    /// <summary>Signal peaks bright enough that no noise draw can be confused for one.</summary>
    private static MsDataScan[] BrightSignalScans(int scanCount = 4) =>
        Enumerable.Range(0, scanCount)
            .Select(s => BuildScan(s + 1, 30.0 + s * 0.01,
                new[] { 800.0, 800.5, 801.0 },
                new[] { 1e9, 8e8, 5e8 }))
            .ToArray();

    [Test]
    public void InjectionAddsRoughlyTheExpectedNumberOfPeaks()
    {
        var model = new NoiseFloorModel(densityScale: 0.01);
        var scans = BrightSignalScans();

        var (result, summary) = new NoiseInjector(model, Width).Apply(scans);

        Assert.That(summary.SignalPeaks, Is.EqualTo(scans.Length * 3));
        Assert.That(summary.NoisePeaks,
            Is.EqualTo(scans.Length * model.ExpectedPeaksPerScan).Within(0.3 * scans.Length * model.ExpectedPeaksPerScan));
        Assert.That(result.Length, Is.EqualTo(scans.Length));

        foreach (var scan in result)
            Assert.That(scan.MassSpectrum.XArray.Length, Is.GreaterThan(3));
    }

    [Test]
    public void SignalPeaksSurviveInjectionIntactWhenNothingLandsOnThem()
    {
        // A sparse noise floor makes a collision with one of the three signal peaks very unlikely.
        var model = new NoiseFloorModel(densityScale: 0.001);
        var scans = BrightSignalScans(1);

        var (result, _) = new NoiseInjector(model, Width).Apply(scans);
        var spectrum = result[0].MassSpectrum;

        foreach (var (mz, intensity) in new[] { (800.0, 1e9), (800.5, 8e8), (801.0, 5e8) })
        {
            int idx = Array.FindIndex(spectrum.XArray, x => Math.Abs(x - mz) < 1e-6);
            Assert.That(idx, Is.GreaterThanOrEqualTo(0), $"signal peak at {mz} is missing");
            Assert.That(spectrum.YArray[idx], Is.EqualTo(intensity).Within(1e-3));
        }
    }

    [Test]
    public void OutputSpectraAreSortedAndTicIsConsistent()
    {
        var model = new NoiseFloorModel(densityScale: 0.05);
        var (result, _) = new NoiseInjector(model, Width).Apply(BrightSignalScans());

        foreach (var scan in result)
        {
            var x = scan.MassSpectrum.XArray;
            for (int i = 1; i < x.Length; i++)
                Assert.That(x[i], Is.GreaterThan(x[i - 1]), "output spectrum is not strictly ascending");

            Assert.That(scan.TotalIonCurrent, Is.EqualTo(scan.MassSpectrum.YArray.Sum()).Within(1e-6 * scan.TotalIonCurrent));
            Assert.That(scan.IsCentroid, Is.True);
        }
    }

    [Test]
    public void ScanMetadataIsPreserved()
    {
        var scans = BrightSignalScans(2);
        var (result, _) = new NoiseInjector(new NoiseFloorModel(densityScale: 0.01), Width).Apply(scans);

        for (int s = 0; s < scans.Length; s++)
        {
            Assert.That(result[s].OneBasedScanNumber, Is.EqualTo(scans[s].OneBasedScanNumber));
            Assert.That(result[s].RetentionTime, Is.EqualTo(scans[s].RetentionTime));
            Assert.That(result[s].MsnOrder, Is.EqualTo(scans[s].MsnOrder));
            Assert.That(result[s].NativeId, Is.EqualTo(scans[s].NativeId));
            Assert.That(result[s].Polarity, Is.EqualTo(scans[s].Polarity));
        }
    }

    [Test]
    public void InjectionIsReproducibleForAGivenSeedAndIndependentOfSeedOtherwise()
    {
        var model = new NoiseFloorModel(densityScale: 0.02);
        var scans = BrightSignalScans(3);

        var (first, _) = new NoiseInjector(model, Width, seed: 11).Apply(scans);
        var (second, _) = new NoiseInjector(model, Width, seed: 11).Apply(scans);
        var (other, _) = new NoiseInjector(model, Width, seed: 12).Apply(scans);

        for (int s = 0; s < scans.Length; s++)
        {
            Assert.That(second[s].MassSpectrum.XArray, Is.EqualTo(first[s].MassSpectrum.XArray));
            Assert.That(second[s].MassSpectrum.YArray, Is.EqualTo(first[s].MassSpectrum.YArray));
        }

        Assert.That(other[0].MassSpectrum.XArray, Is.Not.EqualTo(first[0].MassSpectrum.XArray));
    }

    /// <summary>
    /// Different scans must get different noise, or every XIC in the file would sit on an identical
    /// baseline and a feature finder would see a perfectly flat, perfectly correlated background.
    /// </summary>
    [Test]
    public void DifferentScansGetDifferentNoise()
    {
        var model = new NoiseFloorModel(densityScale: 0.02);
        var (result, _) = new NoiseInjector(model, Width).Apply(BrightSignalScans(2));

        Assert.That(result[1].MassSpectrum.XArray, Is.Not.EqualTo(result[0].MassSpectrum.XArray));
    }

    [Test]
    public void UnresolvablePeaksAreCollapsedIntoOneCentroid()
    {
        // Two peaks a hundredth of a σ apart cannot be reported separately by the instrument.
        double sigma = Width.SigmaAt(800);
        var scans = new[]
        {
            BuildScan(1, 30.0, new[] { 800.0, 800.0 + sigma * 0.01 }, new[] { 100.0, 300.0 }),
        };

        // Density zero isolates the merge behaviour from the sampling.
        var model = new NoiseFloorModel(densityScale: 0);
        var (result, summary) = new NoiseInjector(model, Width).Apply(scans);
        var spectrum = result[0].MassSpectrum;

        Assert.That(summary.NoisePeaks, Is.Zero);
        Assert.That(summary.MergedPeaks, Is.EqualTo(1));
        Assert.That(spectrum.XArray.Length, Is.EqualTo(1));
        Assert.That(spectrum.YArray[0], Is.EqualTo(400.0).Within(1e-9), "merged intensity should be the sum");

        // Intensity-weighted position, so the brighter of the two dominates.
        double expectedMz = (800.0 * 100.0 + (800.0 + sigma * 0.01) * 300.0) / 400.0;
        Assert.That(spectrum.XArray[0], Is.EqualTo(expectedMz).Within(1e-12));
    }

    /// <summary>
    /// A dense run of peaks each just inside the merge tolerance must not chain into one centroid
    /// spanning many σ. At the measured noise density the mean spacing near m/z 700 is about 2σ, so
    /// single-linkage chaining would swallow most of a spectrum.
    /// </summary>
    [Test]
    public void MergingDoesNotChainAcrossALongRunOfClosePeaks()
    {
        double sigma = Width.SigmaAt(800);
        double step = sigma * 0.9;
        var mz = Enumerable.Range(0, 50).Select(i => 800.0 + i * step).ToArray();
        var intensity = Enumerable.Repeat(100.0, mz.Length).ToArray();

        var (result, _) = new NoiseInjector(new NoiseFloorModel(densityScale: 0), Width)
            .Apply(new[] { BuildScan(1, 30.0, mz, intensity) });

        var x = result[0].MassSpectrum.XArray;

        // Chaining would give one peak. Bounded grouping pairs them up, so roughly half survive.
        Assert.That(x.Length, Is.GreaterThan(10), "peaks chained into too few centroids");
        Assert.That(x[^1] - x[0], Is.GreaterThan(30 * step), "the output no longer spans the input range");
    }

    /// <summary>
    /// Groups made purely of noise pass through uncollapsed: the density curve they were drawn from
    /// was measured from centroids a real instrument had already peak-detected, so merging here
    /// would apply that reduction a second time.
    /// </summary>
    [Test]
    public void NoiseOnlyGroupsAreNotCollapsedButGroupsContainingSignalAre()
    {
        // No signal peaks at all, and a density high enough that neighbours land inside the merge
        // tolerance constantly.
        var empty = BuildScan(1, 30.0, Array.Empty<double>(), Array.Empty<double>());
        var model = new NoiseFloorModel();
        var (noiseOnly, noiseSummary) = new NoiseInjector(model, Width).Apply(new[] { empty });

        Assert.That(noiseSummary.MergedPeaks, Is.Zero, "noise was merged against noise");
        Assert.That(noiseOnly[0].MassSpectrum.XArray.Length, Is.EqualTo(noiseSummary.NoisePeaks));

        // With a signal peak present, noise landing on it is absorbed.
        var withSignal = BuildScan(1, 30.0, new[] { 700.0 }, new[] { 1e9 });
        var (_, mixedSummary) = new NoiseInjector(model, Width).Apply(new[] { withSignal });
        Assert.That(mixedSummary.MergedPeaks, Is.GreaterThan(0),
            "noise coincident with a signal peak should have been absorbed");
    }

    [Test]
    public void ResolvablePeaksAreNotCollapsed()
    {
        double sigma = Width.SigmaAt(800);
        var scans = new[]
        {
            BuildScan(1, 30.0, new[] { 800.0, 800.0 + sigma * 3 }, new[] { 100.0, 300.0 }),
        };

        var (result, summary) = new NoiseInjector(new NoiseFloorModel(densityScale: 0), Width).Apply(scans);

        Assert.That(summary.MergedPeaks, Is.Zero);
        Assert.That(result[0].MassSpectrum.XArray.Length, Is.EqualTo(2));
    }

    [Test]
    public void ZeroDensityLeavesScansUnchanged()
    {
        var scans = BrightSignalScans(2);
        var (result, summary) = new NoiseInjector(new NoiseFloorModel(densityScale: 0), Width).Apply(scans);

        Assert.That(summary.NoisePeaks, Is.Zero);
        for (int s = 0; s < scans.Length; s++)
        {
            Assert.That(result[s].MassSpectrum.XArray, Is.EqualTo(scans[s].MassSpectrum.XArray));
            Assert.That(result[s].MassSpectrum.YArray, Is.EqualTo(scans[s].MassSpectrum.YArray));
        }
    }

    [Test]
    public void ConstructorRejectsMissingArguments()
    {
        Assert.Throws<ArgumentNullException>(() => new NoiseInjector(null!, Width));
        Assert.Throws<ArgumentNullException>(() => new NoiseInjector(new NoiseFloorModel(), null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => new NoiseInjector(new NoiseFloorModel(), Width, mergeWithinSigmas: -1));
    }
}
