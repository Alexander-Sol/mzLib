using System;
using System.Linq;
using NUnit.Framework;
using TopDownSimulator.Extraction;
using TopDownSimulator.Fitting;
using TopDownSimulator.Model;

namespace Test.TopDownSimulator;

[TestFixture]
public class ChargeDistributionFitterTests
{
    private static ProteoformGroundTruth BuildTruth(int minCharge, double[] perChargeApex)
    {
        int nC = perChargeApex.Length;
        // Encode apex as a constant single-scan XIC per charge.
        int nS = 3;
        var chargeXics = new double[nC][];
        for (int c = 0; c < nC; c++)
        {
            chargeXics[c] = new double[nS];
            chargeXics[c][1] = perChargeApex[c];
        }
        var centroidMzs = new double[nC][];
        var intensities = new double[nC][][];
        var windows = new PeakSample[nC][][][];
        for (int c = 0; c < nC; c++)
        {
            centroidMzs[c] = new[] { 500.0 };
            intensities[c] = new double[1][];
            intensities[c][0] = new double[nS];
            windows[c] = new PeakSample[1][][];
            windows[c][0] = new PeakSample[nS][];
            for (int s = 0; s < nS; s++) windows[c][0][s] = Array.Empty<PeakSample>();
        }

        return new ProteoformGroundTruth
        {
            MonoisotopicMass = 10000,
            RetentionTimeCenter = 0,
            MinCharge = minCharge,
            MaxCharge = minCharge + nC - 1,
            ZeroBasedScanIndices = new[] { 0, 1, 2 },
            ScanTimes = new[] { 0.0, 1.0, 2.0 },
            CentroidMzs = centroidMzs,
            IsotopologueIntensities = intensities,
            IsotopologuePeakWindows = windows,
            ChargeXics = chargeXics,
            MzWindowHalfWidth = 0.05,
        };
    }

    [Test]
    public void SymmetricGaussianOnChargeRecoversMuSigma()
    {
        const double trueMu = 9.0, trueSigma = 1.2;
        int minZ = 5, maxZ = 13;
        int nC = maxZ - minZ + 1;
        var apex = new double[nC];
        for (int c = 0; c < nC; c++)
        {
            int z = minZ + c;
            double d = (z - trueMu) / trueSigma;
            apex[c] = Math.Exp(-0.5 * d * d);
        }
        var truth = BuildTruth(minZ, apex);

        var fit = new ChargeDistributionFitter().Fit(truth);

        Assert.That(fit.Distribution.MuZ, Is.EqualTo(trueMu).Within(1e-6));
        Assert.That(fit.Distribution.SigmaZ, Is.EqualTo(trueSigma).Within(5e-2));
        Assert.That(fit.ChargesUsed, Is.EqualTo(nC));
    }

    private static double[] Gaussian(int minZ, int maxZ, double mu, double sigma) =>
        Enumerable.Range(minZ, maxZ - minZ + 1).Select(z => Math.Exp(-0.5 * Math.Pow((z - mu) / sigma, 2))).ToArray();

    [Test]
    public void TrimmedFitIgnoresNoiseInTheOuterCharges()
    {
        // A wide window: the envelope at 18 ± 2, and a noise floor plus another species' peak far out.
        const int minZ = 5, maxZ = 40;
        var apex = Gaussian(minZ, maxZ, 18.3, 2.0).Select(a => a + 0.01).ToArray();
        apex[34 - minZ] += 0.4;

        var moments = new ChargeDistributionFitter().Fit(BuildTruth(minZ, apex)).Distribution;
        var trimmed = new ChargeDistributionFitter(trimFraction: 0.05).Fit(BuildTruth(minZ, apex)).Distribution;

        Assert.That(Math.Abs(moments.MuZ - 18.3), Is.GreaterThan(1.0));
        Assert.That(trimmed.MuZ, Is.EqualTo(18.3).Within(0.1));
        Assert.That(trimmed.SigmaZ, Is.EqualTo(2.0).Within(0.15));
    }

    [Test]
    public void TrimmedFitStopsAtANeighbouringEnvelope()
    {
        // Another species' envelope at lower charge, overlapping this one's flank above the trim floor.
        var apex = Gaussian(5, 30, 18.0, 2.0).Zip(Gaussian(5, 30, 10.0, 1.2), (a, b) => a + 0.4 * b).ToArray();

        var trimmed = new ChargeDistributionFitter(trimFraction: 0.05).Fit(BuildTruth(5, apex)).Distribution;

        Assert.That(trimmed.MuZ, Is.EqualTo(18.0).Within(0.1));
        Assert.That(trimmed.SigmaZ, Is.EqualTo(2.0).Within(0.1));
    }

    [Test]
    public void TrimmedFitIsExactForATruncatedGaussian()
    {
        // The window ends just above the envelope's centre, as a members-only window can.
        var apex = Gaussian(10, 16, 15.4, 1.8);

        var moments = new ChargeDistributionFitter().Fit(BuildTruth(10, apex)).Distribution;
        var trimmed = new ChargeDistributionFitter(trimFraction: 0.05).Fit(BuildTruth(10, apex)).Distribution;

        Assert.That(moments.MuZ, Is.LessThan(14.8));
        Assert.That(trimmed.MuZ, Is.EqualTo(15.4).Within(1e-6));
        Assert.That(trimmed.SigmaZ, Is.EqualTo(1.8).Within(1e-6));
    }

    [Test]
    public void SingleChargeFallsBackOnSigma()
    {
        var truth = BuildTruth(minCharge: 8, perChargeApex: new[] { 0.0, 100.0, 0.0 });
        var fit = new ChargeDistributionFitter(fallbackSigmaZ: 1.5).Fit(truth);

        Assert.That(fit.Distribution.MuZ, Is.EqualTo(9.0).Within(1e-9));
        Assert.That(fit.Distribution.SigmaZ, Is.EqualTo(1.5));
        Assert.That(fit.ChargesUsed, Is.EqualTo(1));
    }
}
