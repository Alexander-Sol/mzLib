using System;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using TopDownSimulator.Model;
using TopDownSimulator.Noise;

namespace Test.TopDownSimulator;

[TestFixture]
public class PeakJitterTests
{
    private static readonly IPeakWidthModel Width = OrbitrapPeakWidth.FromSigmaAt(1000, 0.01);

    /// <summary>Density zero isolates jitter from the noise floor.</summary>
    private static NoiseFloorModel QuietModel(double level = 1000) =>
        new(noiseLevelAtReferenceMz: level, densityScale: 0);

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

    /// <summary>
    /// One peak repeated across many scans at constant true intensity, which is what lets the
    /// recovered scatter be compared against the model directly.
    /// </summary>
    private static MsDataScan[] RepeatedPeak(double mz, double intensity, int scans = 4000) =>
        Enumerable.Range(0, scans)
            .Select(s => BuildScan(s + 1, 30.0 + s * 0.01, new[] { mz }, new[] { intensity }))
            .ToArray();

    private static double StandardDeviation(double[] values)
    {
        double mean = values.Average();
        return Math.Sqrt(values.Sum(v => (v - mean) * (v - mean)) / (values.Length - 1));
    }

    [Test]
    public void MeasuredLawsReproduceTheMeasuredScatter()
    {
        var jitter = PeakJitterModel.Orbitrap;

        // The six S/N bins measured on rep2 fract7, and the per-peak m/z scatter in each.
        foreach (var (sn, measured) in new[]
                 { (2.83, 5.447), (5.66, 3.842), (11.3, 2.766), (22.6, 1.780), (45.3, 1.579), (128.0, 1.137) })
            Assert.That(jitter.MzPpmSigma(sn), Is.EqualTo(measured).Within(0.25 * measured), $"m/z sigma at S/N {sn}");

        // And the log-intensity scatter in the same bins.
        foreach (var (sn, measured) in new[]
                 { (2.83, 0.6318), (5.66, 0.4622), (11.3, 0.3644), (22.6, 0.2883), (45.3, 0.2368), (128.0, 0.1881) })
            Assert.That(jitter.LogIntensitySigma(sn), Is.EqualTo(measured).Within(0.10 * measured),
                $"log-intensity sigma at S/N {sn}");
    }

    [Test]
    public void JitterFallsWithSignalToNoiseAndFlattensOntoItsFloor()
    {
        var jitter = PeakJitterModel.Orbitrap;

        Assert.That(jitter.MzPpmSigma(4), Is.GreaterThan(jitter.MzPpmSigma(64)));
        Assert.That(jitter.LogIntensitySigma(4), Is.GreaterThan(jitter.LogIntensitySigma(64)));

        Assert.That(jitter.MzPpmSigma(1e12), Is.EqualTo(jitter.MzPpmFloor).Within(1e-6));
        Assert.That(jitter.LogIntensitySigma(1e12), Is.EqualTo(jitter.LogIntensitySigmaFloor).Within(1e-6));

        // Non-positive or non-finite S/N degrades to the floor rather than producing NaN.
        Assert.That(jitter.MzPpmSigma(0), Is.EqualTo(jitter.MzPpmFloor));
        Assert.That(jitter.LogIntensitySigma(double.NaN), Is.EqualTo(jitter.LogIntensitySigmaFloor));
    }

    /// <summary>
    /// The scatter actually produced end to end must match what the law promises, not merely be
    /// non-zero. Uses a bright peak so the floor terms dominate and dropout cannot bias the sample.
    /// </summary>
    [Test]
    public void RecoveredScatterMatchesTheModelledSigma()
    {
        const double mz = 900.0;
        const double level = 1000.0;
        const double intensity = 1e5;

        var model = QuietModel(level);
        // Common mode off, to isolate the per-peak terms.
        var jitter = PeakJitterModel.Orbitrap with { CommonModeMzPpm = 0 };

        // S/N is measured against the level at this m/z, not against the reference-m/z level: the
        // amplitude curve rises about 2.4x between m/z 650 and 900.
        double signalToNoise = intensity / model.NoiseLevelAt(mz);
        Assert.That(signalToNoise, Is.GreaterThan(20), "the test peak needs to be well clear of the floor");

        var (scans, summary) = new NoiseInjector(model, Width, seed: 5, jitter: jitter)
            .Apply(RepeatedPeak(mz, intensity));

        Assert.That(summary.DroppedByJitter, Is.Zero, "a bright peak should never fall below the floor");

        var observedMz = scans.Select(s => s.MassSpectrum.XArray.Single()).ToArray();
        var observedIntensity = scans.Select(s => s.MassSpectrum.YArray.Single()).ToArray();

        double ppmScatter = StandardDeviation(observedMz.Select(x => (x - mz) / mz * 1e6).ToArray());
        double logScatter = StandardDeviation(observedIntensity.Select(y => Math.Log(y / intensity)).ToArray());

        double expectedPpm = jitter.MzPpmSigma(signalToNoise);
        double expectedLog = jitter.LogIntensitySigma(signalToNoise);

        Assert.That(ppmScatter, Is.EqualTo(expectedPpm).Within(0.1 * expectedPpm));
        Assert.That(logScatter, Is.EqualTo(expectedLog).Within(0.1 * expectedLog));
    }

    /// <summary>
    /// Intensity jitter is mean-preserving. A median-preserving draw would inflate expected
    /// intensity by exp(σ²/2), quietly brightening exactly the faint peaks whose detectability the
    /// simulation exists to test.
    /// </summary>
    [Test]
    public void IntensityJitterPreservesTheMean()
    {
        var rng = new Random(17);
        foreach (double sigma in new[] { 0.165, 0.5, 0.84 })
        {
            var factors = Enumerable.Range(0, 200_000)
                .Select(_ => PeakJitterModel.IntensityFactor(rng, sigma))
                .ToArray();

            Assert.That(factors.Average(), Is.EqualTo(1.0).Within(0.02), $"mean at sigma {sigma}");
            // And it really is dispersed, not merely centred.
            Assert.That(StandardDeviation(factors), Is.GreaterThan(0.5 * sigma));
        }

        Assert.That(PeakJitterModel.IntensityFactor(rng, 0), Is.EqualTo(1.0));
    }

    /// <summary>
    /// The common-mode m/z term moves every peak in a scan together, which is what distinguishes
    /// calibration drift from independent centroiding error.
    /// </summary>
    [Test]
    public void CommonModeMzOffsetIsSharedAcrossAScan()
    {
        double[] mz = Enumerable.Range(0, 200).Select(i => 700.0 + i * 3.0).ToArray();
        double[] intensity = Enumerable.Repeat(1e6, mz.Length).ToArray();
        var scans = new[] { BuildScan(1, 30.0, mz, intensity) };

        // Common mode only: per-peak terms zeroed, so all scatter must be the shared offset.
        var sharedOnly = PeakJitterModel.Orbitrap with
        {
            MzPpmAtUnitSignalToNoise = 0,
            MzPpmFloor = 0,
            LogIntensitySigmaAtUnitSignalToNoise = 0,
            LogIntensitySigmaFloor = 0,
        };

        var (result, _) = new NoiseInjector(QuietModel(), Width, seed: 3, jitter: sharedOnly).Apply(scans);
        var deviations = result[0].MassSpectrum.XArray
            .Select((x, i) => (x - mz[i]) / mz[i] * 1e6)
            .ToArray();

        Assert.That(StandardDeviation(deviations), Is.LessThan(1e-6), "peaks did not all shift by the same ppm");
        Assert.That(Math.Abs(deviations[0]), Is.GreaterThan(0), "no common-mode shift was applied at all");
    }

    /// <summary>
    /// A peak just above the detection floor genuinely disappears in a good fraction of scans — its
    /// S/N is about 1.6, so σ_log ≈ 0.84 — while a bright one never does. This XIC dropout is a real
    /// stressor for feature finding, and it comes out of the model rather than being added to it.
    /// </summary>
    [Test]
    public void FaintPeaksDropOutWhileBrightOnesDoNot()
    {
        const double level = 1000.0;
        var model = QuietModel(level);
        double floor = new NoiseRelativeFloor(model).At(900.0);

        var (_, faintSummary) = new NoiseInjector(model, Width, seed: 9, jitter: PeakJitterModel.Orbitrap)
            .Apply(RepeatedPeak(900.0, floor * 1.05, scans: 2000));
        var (_, brightSummary) = new NoiseInjector(model, Width, seed: 9, jitter: PeakJitterModel.Orbitrap)
            .Apply(RepeatedPeak(900.0, floor * 100, scans: 2000));

        double faintDropRate = faintSummary.DroppedByJitter / 2000.0;
        TestContext.Out.WriteLine($"faint dropout {faintDropRate:P1}, bright dropout {brightSummary.DroppedByJitter / 2000.0:P1}");

        Assert.That(faintDropRate, Is.GreaterThan(0.2).And.LessThan(0.8));
        Assert.That(brightSummary.DroppedByJitter, Is.Zero);
    }

    [Test]
    public void DropoutCanBeDisabled()
    {
        const double level = 1000.0;
        var model = QuietModel(level);
        double floor = new NoiseRelativeFloor(model).At(900.0);
        var keepAll = PeakJitterModel.Orbitrap with { DropPeaksFallingBelowFloor = false };

        var (scans, summary) = new NoiseInjector(model, Width, seed: 9, jitter: keepAll)
            .Apply(RepeatedPeak(900.0, floor * 1.05, scans: 500));

        Assert.That(summary.DroppedByJitter, Is.Zero);
        Assert.That(scans.All(s => s.MassSpectrum.XArray.Length == 1));
    }

    [Test]
    public void NoneLeavesPeaksExactlyAlone()
    {
        var scans = RepeatedPeak(900.0, 1e6, scans: 20);
        var (result, summary) = new NoiseInjector(QuietModel(), Width, jitter: PeakJitterModel.None).Apply(scans);

        Assert.That(summary.DroppedByJitter, Is.Zero);
        for (int s = 0; s < scans.Length; s++)
        {
            Assert.That(result[s].MassSpectrum.XArray, Is.EqualTo(scans[s].MassSpectrum.XArray));
            Assert.That(result[s].MassSpectrum.YArray, Is.EqualTo(scans[s].MassSpectrum.YArray));
        }
    }

    [Test]
    public void NoiseFloorPeaksAreNotJittered()
    {
        // Noise is drawn from a distribution measured on already-reported centroids, so its scatter
        // is baked in; jittering it again would double-count. With no signal present, switching
        // jitter on must change nothing.
        var empty = new[] { BuildScan(1, 30.0, Array.Empty<double>(), Array.Empty<double>()) };
        var model = new NoiseFloorModel(densityScale: 0.05);

        var (withoutJitter, _) = new NoiseInjector(model, Width, seed: 4).Apply(empty);
        var (withJitter, _) = new NoiseInjector(model, Width, seed: 4, jitter: PeakJitterModel.Orbitrap).Apply(empty);

        Assert.That(withJitter[0].MassSpectrum.XArray, Is.EqualTo(withoutJitter[0].MassSpectrum.XArray));
        Assert.That(withJitter[0].MassSpectrum.YArray, Is.EqualTo(withoutJitter[0].MassSpectrum.YArray));
    }

    [Test]
    public void JitterIsReproducibleForAGivenSeed()
    {
        var scans = RepeatedPeak(900.0, 1e6, scans: 10);

        var (first, _) = new NoiseInjector(QuietModel(), Width, seed: 21, jitter: PeakJitterModel.Orbitrap).Apply(scans);
        var (second, _) = new NoiseInjector(QuietModel(), Width, seed: 21, jitter: PeakJitterModel.Orbitrap).Apply(scans);
        var (other, _) = new NoiseInjector(QuietModel(), Width, seed: 22, jitter: PeakJitterModel.Orbitrap).Apply(scans);

        Assert.That(second[0].MassSpectrum.XArray, Is.EqualTo(first[0].MassSpectrum.XArray));
        Assert.That(other[0].MassSpectrum.XArray, Is.Not.EqualTo(first[0].MassSpectrum.XArray));

        // Consecutive scans must differ, or every XIC would sit on an identical, perfectly
        // correlated baseline.
        Assert.That(first[1].MassSpectrum.YArray.Single(), Is.Not.EqualTo(first[0].MassSpectrum.YArray.Single()));
    }

    /// <summary>
    /// The point of the exercise: a simulated XIC must not be perfectly smooth any more.
    /// </summary>
    [Test]
    public void SimulatedXicIsNoLongerPerfectlySmooth()
    {
        const double level = 1000.0;
        var model = QuietModel(level);
        var scans = RepeatedPeak(900.0, 30 * level, scans: 200);

        var (clean, _) = new NoiseInjector(model, Width, seed: 2, jitter: PeakJitterModel.None).Apply(scans);
        var (jittered, _) = new NoiseInjector(model, Width, seed: 2, jitter: PeakJitterModel.Orbitrap).Apply(scans);

        double CleanRoughness(MsDataScan[] xs) => StandardDeviation(
            Enumerable.Range(1, xs.Length - 1)
                .Select(i => xs[i].MassSpectrum.YArray.Single() - xs[i - 1].MassSpectrum.YArray.Single())
                .ToArray());

        Assert.That(CleanRoughness(clean), Is.EqualTo(0).Within(1e-9), "the unjittered XIC should be exactly flat");
        Assert.That(CleanRoughness(jittered), Is.GreaterThan(0.1 * 30 * level),
            "the jittered XIC should carry scan-to-scan structure");
    }
}
