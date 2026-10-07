using System;
using System.Linq;
using NUnit.Framework;
using TopDownSimulator.Comparison;
using TopDownSimulator.Extraction;
using TopDownSimulator.Noise;

namespace Test.TopDownSimulator;

[TestFixture]
public class SignalComparisonTests
{
    [Test]
    public void SignalPeaksKeepsPeaksAboveTheModelledNoise()
    {
        var noise = new NoiseFloorModel(500.0);
        var mz = new[] { 700.0, 800.0, 900.0, 1000.0 };
        // At S/N 3, 12, 9.9 and 20 against the model's level at each m/z.
        var sn = new[] { 3.0, 12, 9.9, 20 };
        var intensity = mz.Select((m, i) => sn[i] * noise.NoiseLevelAt(m)).ToArray();

        var (keptMz, keptIntensity) = SignalPeaks.Select(mz, intensity, noise, minSignalToNoise: 10);

        Assert.That(keptMz, Is.EqualTo(new[] { 800.0, 1000.0 }));
        Assert.That(keptIntensity, Is.EqualTo(new[] { intensity[1], intensity[3] }));
    }

    [Test]
    public void SignalPeaksRespectsTheMzRange()
    {
        var noise = new NoiseFloorModel(1.0);
        var mz = new[] { 500.0, 700.0, 2500.0 };
        var intensity = mz.Select(m => 100 * noise.NoiseLevelAt(m)).ToArray();
        var (keptMz, _) = SignalPeaks.Select(mz, intensity, noise, 10, 600, 2000);
        Assert.That(keptMz, Is.EqualTo(new[] { 700.0 }));
    }

    private static ProteoformGroundTruth Truth(double[][] perChargeIsotopologue)
    {
        int nC = perChargeIsotopologue.Length;
        return new ProteoformGroundTruth
        {
            MonoisotopicMass = 10000,
            RetentionTimeCenter = 0,
            MinCharge = 8,
            MaxCharge = 8 + nC - 1,
            ZeroBasedScanIndices = new[] { 0, 1 },
            ScanTimes = new[] { 0.0, 0.01 },
            CentroidMzs = perChargeIsotopologue.Select(iso => iso.Select((_, k) => 1000.0 + k).ToArray()).ToArray(),
            // Split each isotopologue over two scans so the comparison has to sum them.
            IsotopologueIntensities = perChargeIsotopologue.Select(iso => iso.Select(v => new[] { 0.25 * v, 0.75 * v }).ToArray()).ToArray(),
            IsotopologuePeakWindows = perChargeIsotopologue.Select(iso => iso.Select(_ => new[] { Array.Empty<PeakSample>(), Array.Empty<PeakSample>() }).ToArray()).ToArray(),
            ChargeXics = perChargeIsotopologue.Select(_ => new double[2]).ToArray(),
            MzWindowHalfWidth = 0.05,
        };
    }

    [Test]
    public void EnvelopeComparisonSeparatesShapeChargeAndScale()
    {
        var real = Truth(new[] { new[] { 1.0, 2, 1 }, new[] { 2.0, 4, 2 }, new[] { 1.0, 2, 1 } });
        var scaled = Truth(new[] { new[] { 3.0, 6, 3 }, new[] { 6.0, 12, 6 }, new[] { 3.0, 6, 3 } });
        var shifted = Truth(new[] { new[] { 0.0, 0, 0 }, new[] { 1.0, 2, 1 }, new[] { 2.0, 4, 2 } });

        var same = SpeciesEnvelopeComparison.Compare(real, scaled)!;
        Assert.That(same.Cosine, Is.EqualTo(1).Within(1e-12));
        Assert.That(same.ChargeProfileCosine, Is.EqualTo(1).Within(1e-12));
        Assert.That(same.IntensityRatio, Is.EqualTo(3).Within(1e-12));
        Assert.That((same.RealApexCharge, same.SimulatedApexCharge), Is.EqualTo((9, 9)));

        var off = SpeciesEnvelopeComparison.Compare(real, shifted)!;
        Assert.That(off.ChargeProfileCosine, Is.LessThan(0.9));
        Assert.That(off.SimulatedApexCharge, Is.EqualTo(10));
    }

    [Test]
    public void EnvelopeComparisonOfAnEmptyExtractionIsNull()
    {
        var real = Truth(new[] { new[] { 1.0, 2, 1 } });
        var empty = Truth(new[] { new[] { 0.0, 0, 0 } });
        Assert.That(SpeciesEnvelopeComparison.Compare(real, empty), Is.Null);
    }
}
