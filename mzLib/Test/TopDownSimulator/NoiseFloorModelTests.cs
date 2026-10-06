using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using TopDownSimulator.Model;
using TopDownSimulator.Noise;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

[TestFixture]
public class NoiseFloorModelTests
{
    /// <summary>
    /// The noise amplitude measured in the 50-55 min window of rep1 fract7, the file the fidelity
    /// checks below compare against.
    /// </summary>
    private const double Rep1Fract7NoiseLevel = 715.0;

    private static List<SimulatedPeak> SampleScans(NoiseFloorModel model, int scanCount, int seed = 1)
    {
        var peaks = new List<SimulatedPeak>();
        for (int s = 0; s < scanCount; s++)
            model.SampleScan(NoiseInjector.StreamForScan(seed, s), peaks);
        return peaks;
    }

    private static double Quantile(IEnumerable<double> values, double q)
    {
        var sorted = values.ToArray();
        Array.Sort(sorted);
        int idx = (int)Math.Round(q * (sorted.Length - 1));
        return sorted[Math.Clamp(idx, 0, sorted.Length - 1)];
    }

    [Test]
    public void ExpectedDensityMatchesTheMeasuredPeakCount()
    {
        var model = new NoiseFloorModel();

        // Measured median peaks/scan in the 50-55 min window was 16 434 - 17 554 across four files.
        Assert.That(model.ExpectedPeaksPerScan, Is.EqualTo(16886).Within(0.5));
    }

    [Test]
    public void SampledScanCountIsPoissonAboutTheExpectedDensity()
    {
        var model = new NoiseFloorModel();
        var counts = Enumerable.Range(0, 20)
            .Select(s => SampleScans(model, 1, seed: s).Count)
            .ToArray();

        double expected = model.ExpectedPeaksPerScan;
        double tolerance = 6 * Math.Sqrt(expected);
        Assert.That(counts.Average(), Is.EqualTo(expected).Within(tolerance));
        Assert.That(counts.Min(), Is.GreaterThan(expected - tolerance));
        Assert.That(counts.Max(), Is.LessThan(expected + tolerance));
    }

    [Test]
    public void DensityScaleIsLinearInPeakCount()
    {
        var full = new NoiseFloorModel();
        var tenth = new NoiseFloorModel(densityScale: 0.1);

        Assert.That(tenth.ExpectedPeaksPerScan, Is.EqualTo(full.ExpectedPeaksPerScan * 0.1).Within(1e-9));
        Assert.That(SampleScans(tenth, 10).Count,
            Is.EqualTo(SampleScans(full, 10).Count * 0.1).Within(0.2 * SampleScans(full, 10).Count * 0.1));
    }

    [Test]
    public void SampledPeaksAreSortedAndInsideTheCalibratedRange()
    {
        var model = new NoiseFloorModel();
        var peaks = SampleScans(model, 1);

        Assert.That(peaks, Is.Not.Empty);
        for (int i = 1; i < peaks.Count; i++)
            Assert.That(peaks[i].Mz, Is.GreaterThanOrEqualTo(peaks[i - 1].Mz), $"unsorted at index {i}");

        Assert.That(peaks.Min(p => p.Mz), Is.GreaterThanOrEqualTo(model.MinMz));
        Assert.That(peaks.Max(p => p.Mz), Is.LessThanOrEqualTo(model.MaxMz));
        Assert.That(peaks.All(p => p.Origin == PeakOrigin.FtNoise));
        Assert.That(peaks.All(p => p.Intensity > 0));
    }

    /// <summary>
    /// Noise peak density falls by a factor of ~430 from its maximum near m/z 750 out to m/z 1950.
    /// A model that got this wrong would put its noise where the instrument does not.
    /// </summary>
    [Test]
    public void PeakDensityFallsSteeplyWithMz()
    {
        var model = new NoiseFloorModel();
        var peaks = SampleScans(model, 5);

        double PerScanIn(double lo, double hi) => peaks.Count(p => p.Mz >= lo && p.Mz < hi) / 5.0;

        Assert.That(PerScanIn(700, 800), Is.EqualTo(4303).Within(0.15 * 4303));
        Assert.That(PerScanIn(1000, 1100), Is.EqualTo(1337).Within(0.15 * 1337));
        Assert.That(PerScanIn(1900, 2000), Is.EqualTo(8.2).Within(0.5 * 8.2));

        // The maximum is at 700-800, not at the low edge: density rises from 600 before it falls.
        Assert.That(PerScanIn(700, 800), Is.GreaterThan(PerScanIn(600, 700)));
        Assert.That(PerScanIn(700, 800), Is.GreaterThan(PerScanIn(800, 900)));
    }

    [Test]
    public void NoiseLevelRisesToAPlateauNearMz1450()
    {
        var model = new NoiseFloorModel(noiseLevelAtReferenceMz: 1000);

        Assert.That(model.NoiseLevelAt(650), Is.EqualTo(1000).Within(1e-9));
        Assert.That(model.NoiseLevelAt(1450), Is.EqualTo(4630).Within(1e-9));

        // Held flat outside the calibrated range rather than extrapolated.
        Assert.That(model.NoiseLevelAt(100), Is.EqualTo(model.NoiseLevelAt(650)).Within(1e-9));
        Assert.That(model.NoiseLevelAt(5000), Is.EqualTo(model.NoiseLevelAt(1950)).Within(1e-9));

        // Unimodal, which NoiseRelativeFloor.MinOver relies on.
        var levels = Enumerable.Range(0, 140).Select(i => model.NoiseLevelAt(600 + i * 10)).ToArray();
        int peakIndex = Array.IndexOf(levels, levels.Max());
        for (int i = 1; i <= peakIndex; i++)
            Assert.That(levels[i], Is.GreaterThanOrEqualTo(levels[i - 1]), $"not rising at {600 + i * 10}");
        for (int i = peakIndex + 1; i < levels.Length; i++)
            Assert.That(levels[i], Is.LessThanOrEqualTo(levels[i - 1]), $"not falling at {600 + i * 10}");
    }

    [Test]
    public void SampledSignalToNoiseMatchesItsAnalyticQuantiles()
    {
        var distribution = SignalToNoiseDistribution.Orbitrap;
        var rng = new Random(42);
        var draws = Enumerable.Range(0, 200_000).Select(_ => distribution.Sample(rng)).ToArray();

        foreach (double q in new[] { 0.25, 0.50, 0.75, 0.95, 0.99 })
            Assert.That(Quantile(draws, q), Is.EqualTo(distribution.Quantile(q)).Within(0.03 * distribution.Quantile(q)),
                $"quantile {q}");

        Assert.That(draws.Min(), Is.GreaterThanOrEqualTo(distribution.Threshold));
    }

    /// <summary>
    /// The measured S/N distribution the default was fitted to, from the quiet window of four Jurkat
    /// runs. These quantiles were stable across every file and both RT windows.
    /// </summary>
    [Test]
    public void DefaultSignalToNoiseReproducesTheMeasuredQuantiles()
    {
        var distribution = SignalToNoiseDistribution.Orbitrap;

        Assert.That(distribution.Quantile(0.25), Is.EqualTo(2.09).Within(0.10));
        Assert.That(distribution.Quantile(0.50), Is.EqualTo(2.75).Within(0.10));
        Assert.That(distribution.Quantile(0.75), Is.EqualTo(3.90).Within(0.15));
        Assert.That(distribution.Quantile(0.95), Is.EqualTo(6.68).Within(0.50));
    }

    /// <summary>
    /// The end-to-end fidelity check: the <i>shape</i> of the simulated noise intensity distribution
    /// against the one measured in rep1 fract7's 50-55 min window. This is what the model exists to
    /// reproduce, and it is the assertion that fails first if the density curve, the level curve,
    /// the S/N distribution, or the dispersion drifts.
    /// </summary>
    /// <remarks>
    /// Compared as ratios to the median rather than in absolute intensity, because the absolute
    /// amplitude is the model's one free parameter — it spanned a factor of 3 across the four
    /// calibration files — and is meant to be set per run. The shape is what is claimed to be
    /// transferable, so the shape is what is tested.
    /// </remarks>
    [Test]
    public void SampledIntensityDistributionHasTheMeasuredShape()
    {
        // Measured intensity quantiles, rep1 fract7, 50-55 min.
        const double measuredP25 = 2667, measuredP50 = 4316, measuredP75 = 7113, measuredP95 = 15730;

        var model = new NoiseFloorModel(noiseLevelAtReferenceMz: Rep1Fract7NoiseLevel);
        var intensities = SampleScans(model, 5).Select(p => p.Intensity).ToArray();
        double median = Quantile(intensities, 0.50);

        Assert.Multiple(() =>
        {
            foreach (var (q, measured) in new[]
                     { (0.25, measuredP25), (0.75, measuredP75), (0.95, measuredP95) })
            {
                double simulatedRatio = Quantile(intensities, q) / median;
                double measuredRatio = measured / measuredP50;
                TestContext.Out.WriteLine(
                    $"p{q * 100:F0}/p50: simulated {simulatedRatio:F3}, measured {measuredRatio:F3}, " +
                    $"ratio {simulatedRatio / measuredRatio:F3}");
                Assert.That(simulatedRatio, Is.EqualTo(measuredRatio).Within(0.10 * measuredRatio),
                    $"intensity quantile ratio p{q * 100:F0}/p50");
            }

            // The absolute scale should still land in the right order of magnitude when the measured
            // amplitude for this file is supplied.
            TestContext.Out.WriteLine($"p50: simulated {median:F0}, measured {measuredP50:F0}");
            Assert.That(median, Is.EqualTo(measuredP50).Within(0.25 * measuredP50));
        });
    }

    [Test]
    public void SamplingIsReproducibleForAGivenSeed()
    {
        var model = new NoiseFloorModel();
        var first = SampleScans(model, 3, seed: 7);
        var second = SampleScans(model, 3, seed: 7);
        var different = SampleScans(model, 3, seed: 8);

        Assert.That(second, Is.EqualTo(first));
        Assert.That(different, Is.Not.EqualTo(first));
    }

    [Test]
    public void AdjacentScansAreNotCorrelated()
    {
        var model = new NoiseFloorModel(densityScale: 0.05);
        var a = SampleScans(model, 1, seed: 3).Select(p => p.Mz).ToArray();
        var b = new List<SimulatedPeak>();
        model.SampleScan(NoiseInjector.StreamForScan(3, 1), b);

        // Independent draws, so overlap should be at chance level: with ~850 peaks over 1400 Th and
        // a 10 ppm window, the expected coincidence rate is well under 5 %.
        int shared = a.Count(mz => b.Any(p => Math.Abs(p.Mz - mz) < mz * 1e-5));
        Assert.That(shared / (double)a.Length, Is.LessThan(0.05));
    }

    [Test]
    public void ConstructorRejectsInvalidParameters()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new NoiseFloorModel(noiseLevelAtReferenceMz: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => new NoiseFloorModel(noiseLevelAtReferenceMz: double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => new NoiseFloorModel(densityScale: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SignalToNoiseDistribution(0, 1));
        Assert.Throws<ArgumentOutOfRangeException>(() => new SignalToNoiseDistribution(1, 0));
    }

    [Test]
    public void NoiseRelativeFloorTracksTheNoiseLevel()
    {
        var model = new NoiseFloorModel(noiseLevelAtReferenceMz: 1000);
        var floor = new NoiseRelativeFloor(model);

        Assert.That(floor.SnThreshold, Is.EqualTo(SignalToNoiseDistribution.Orbitrap.Threshold));
        Assert.That(floor.At(650), Is.EqualTo(1000 * floor.SnThreshold).Within(1e-9));
        Assert.That(floor.At(1450), Is.GreaterThan(floor.At(650)));
        Assert.That(floor.MinOver(650, 1450), Is.EqualTo(floor.At(650)).Within(1e-9));
        Assert.That(floor.MinOver(1450, 1950), Is.EqualTo(floor.At(1950)).Within(1e-9));
    }

    /// <summary>
    /// The nominal floor sits in the body of the generated noise rather than underneath all of it.
    /// </summary>
    /// <remarks>
    /// Not every noise peak clears it, and that is correct. A real instrument reports a peak when it
    /// beats a multiple of <i>its own</i> local noise estimate, and that estimate scatters about the
    /// median curve with <see cref="NoiseFloorModel.LevelDispersion"/>. So a peak sitting where the
    /// local noise happens to be low is legitimately dimmer than the threshold times the median
    /// level. What must hold is that the floor lands in the body of the distribution, not above it.
    /// </remarks>
    [Test]
    public void NoiseRelativeFloorSitsWithinTheGeneratedNoise()
    {
        var model = new NoiseFloorModel(noiseLevelAtReferenceMz: 800);
        var floor = new NoiseRelativeFloor(model);

        var peaks = SampleScans(model, 2);
        double aboveFloor = peaks.Count(p => p.Intensity >= floor.At(p.Mz)) / (double)peaks.Count;

        Assert.That(aboveFloor, Is.GreaterThan(0.75), "floor is above the bulk of the noise");
        Assert.That(aboveFloor, Is.LessThan(0.99), "floor is below all of the noise, so dispersion is not being applied");

        // With dispersion switched off, the floor is exactly the lower limit of what is generated.
        var undispersed = new NoiseFloorModel(noiseLevelAtReferenceMz: 800, levelDispersion: 0);
        foreach (var peak in SampleScans(undispersed, 1))
            Assert.That(peak.Intensity, Is.GreaterThanOrEqualTo(new NoiseRelativeFloor(undispersed).At(peak.Mz) * (1 - 1e-9)));
    }
}
