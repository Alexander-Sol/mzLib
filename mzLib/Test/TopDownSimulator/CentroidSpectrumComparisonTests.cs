using NUnit.Framework;
using TopDownSimulator.Comparison;

namespace Test.TopDownSimulator;

[TestFixture]
public class CentroidSpectrumComparisonTests
{
    [Test]
    public void IdenticalSpectraScorePerfectly()
    {
        double[] mz = { 500.0, 500.1, 500.2, 800.0 };
        double[] intensity = { 10, 100, 1000, 5 };

        var c = CentroidSpectrumComparison.Compare(mz, intensity, mz, intensity);

        Assert.That(c.RealPeaks, Is.EqualTo(4));
        Assert.That(c.RealMatchedFraction, Is.EqualTo(1));
        Assert.That(c.SimulatedMatchedFraction, Is.EqualTo(1));
        Assert.That(c.Cosine, Is.EqualTo(1).Within(1e-12));
        Assert.That(c.SqrtCosine, Is.EqualTo(1).Within(1e-12));
        Assert.That(c.TicRatio, Is.EqualTo(1).Within(1e-12));
        Assert.That(c.RelativeIntensityKs, Is.EqualTo(0));
    }

    [Test]
    public void SpuriousFaintPeaksMoveTheSqrtCosineMoreThanTheCosine()
    {
        double[] realMz = { 860.0, 860.0625, 860.125 };
        double[] realI = { 1e6, 2e6, 1e6 };

        // The same envelope with faint peaks inserted between the isotopologues.
        double[] simMz = { 860.0, 860.02, 860.04, 860.0625, 860.08, 860.1, 860.125 };
        double[] simI = { 1e6, 1e4, 1e4, 2e6, 1e4, 1e4, 1e6 };

        var c = CentroidSpectrumComparison.Compare(realMz, realI, simMz, simI);

        Assert.That(c.RealMatchedFraction, Is.EqualTo(1));
        Assert.That(c.SimulatedMatchedFraction, Is.EqualTo(3.0 / 7).Within(1e-12));
        Assert.That(c.Cosine, Is.GreaterThan(0.999));
        Assert.That(c.SqrtCosine, Is.LessThan(c.Cosine));
        Assert.That(c.SimulatedFractionBelowOnePercent, Is.EqualTo(4.0 / 7).Within(1e-12));
        Assert.That(c.RealFractionBelowOnePercent, Is.EqualTo(0));
    }

    [Test]
    public void MatchingRespectsThePpmToleranceAndTheRange()
    {
        double[] realMz = { 1000.0, 1000.5, 2000.0 };
        double[] realI = { 1, 1, 1 };
        double[] simMz = { 1000.005, 1000.52, 2000.0 }; // 5 ppm, 20 ppm, outside range

        var c = CentroidSpectrumComparison.Compare(realMz, realI, simMz, realI, minMz: 900, maxMz: 1500, ppmTolerance: 10);

        Assert.That(c.RealPeaks, Is.EqualTo(2));
        Assert.That(c.SimulatedPeaks, Is.EqualTo(2));
        Assert.That(c.RealMatchedFraction, Is.EqualTo(0.5));
    }

    [Test]
    public void ScaleDoesNotAffectTheShapeMetrics()
    {
        double[] mz = { 500.0, 600.0, 700.0 };
        double[] real = { 1, 10, 100 };
        double[] sim = { 50, 500, 5000 };

        var c = CentroidSpectrumComparison.Compare(mz, real, mz, sim);

        Assert.That(c.Cosine, Is.EqualTo(1).Within(1e-12));
        Assert.That(c.RelativeIntensityKs, Is.EqualTo(0));
        Assert.That(c.TicRatio, Is.EqualTo(50).Within(1e-9));
    }
}
