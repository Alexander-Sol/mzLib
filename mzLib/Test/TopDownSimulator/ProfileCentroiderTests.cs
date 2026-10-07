using System;
using System.Linq;
using Chemistry;
using NUnit.Framework;
using TopDownSimulator.Model;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

[TestFixture]
public class ProfileCentroiderTests
{
    private const double Sigma = 0.012;
    private const double Rt = 20.0;

    private static ProteoformModel Model(double mass, double abundance = 1e6, double muZ = 8.3) =>
        new(mass, abundance, new EmgProfile(Rt, 0.2, 0), new GaussianChargeDistribution(muZ, 1.0));

    private static (double[] Mz, double[] Intensity) Centroid(params ProteoformModel[] models) =>
        new ProfileCentroider(models, 8, 8, new ConstantPeakWidth(Sigma)).Centroid(Rt);

    [Test]
    public void ResolvedEnvelopeGivesOnePeakPerIsotopologueAtItsTheoreticalPosition()
    {
        var model = Model(10000);
        var (mz, intensity) = Centroid(model);
        var kernel = new IsotopeEnvelopeKernel(10000);
        var forward = new ForwardModel(new[] { model }, 8, 8, Sigma);

        // At z=8 the isotopologues are 0.125 Th = 10σ apart, so every one is its own maximum.
        Assert.That(mz, Has.Length.EqualTo(kernel.IsotopologueCount));
        for (int i = 0; i < mz.Length; i++)
        {
            Assert.That(mz[i], Is.EqualTo(kernel.NeutralMass(i).ToMz(8)).Within(1e-6));
            Assert.That(intensity[i], Is.EqualTo(forward.Evaluate(Rt, mz[i])).Within(1e-9).Percent);
        }
    }

    [Test]
    public void ReportedHeightIsTheForwardModelAtALocalMaximum()
    {
        // Two species 0.4 Da apart at z=8: 0.05 Th, about 4σ, so their peaks overlap but stay resolved.
        var models = new[] { Model(10000), Model(10000.4, abundance: 4e5) };
        var forward = new ForwardModel(models, 8, 8, Sigma);
        var (mz, intensity) = Centroid(models);

        for (int i = 0; i < mz.Length; i++)
        {
            double here = forward.Evaluate(Rt, mz[i]);
            Assert.That(intensity[i], Is.EqualTo(here).Within(1e-9).Percent);
            Assert.That(here, Is.GreaterThanOrEqualTo(forward.Evaluate(Rt, mz[i] - 1e-4)));
            Assert.That(here, Is.GreaterThanOrEqualTo(forward.Evaluate(Rt, mz[i] + 1e-4)));
        }
    }

    [Test]
    public void UnresolvedNeighboursMergeIntoOnePeakBetweenThem()
    {
        // 0.08 Da apart at z=8 is 0.01 Th, under one σ: an instrument reports one peak, not two.
        var a = Model(10000);
        var b = Model(10000.08);
        var (alone, _) = Centroid(a);
        var (merged, intensity) = Centroid(a, b);

        Assert.That(merged, Has.Length.EqualTo(alone.Length));
        var kernel = new IsotopeEnvelopeKernel(10000);
        double lo = kernel.NeutralMass(0).ToMz(8);
        Assert.That(merged[0], Is.GreaterThan(lo).And.LessThan(lo + 0.01));
    }

    [Test]
    public void AFaintSpeciesOnTheShoulderOfABrightOneAddsNoPeak()
    {
        // 1.5σ off a peak 100x brighter: the summed profile has no maximum there, so the union
        // sampling this replaced would have written a spurious shoulder peak.
        var bright = Model(10000, abundance: 1e7);
        var faint = Model(10000 + 1.5 * Sigma * 8, abundance: 1e5);

        var (alone, _) = Centroid(bright);
        var (both, _) = Centroid(bright, faint);

        Assert.That(both, Has.Length.EqualTo(alone.Length));
    }

    [Test]
    public void DuplicatedModelsGiveTheSamePeaksAtTwiceTheHeight()
    {
        var model = Model(10000);
        var (mz1, i1) = Centroid(model);
        var (mz2, i2) = Centroid(model, model);

        Assert.That(mz2, Has.Length.EqualTo(mz1.Length));
        for (int i = 0; i < mz1.Length; i++)
        {
            Assert.That(mz2[i], Is.EqualTo(mz1[i]).Within(1e-9));
            Assert.That(i2[i], Is.EqualTo(2 * i1[i]).Within(1e-9).Percent);
        }
    }

    [Test]
    public void ParallelCentroidingMatchesSerial()
    {
        var models = Enumerable.Range(0, 20).Select(i => Model(9000 + 137.3 * i, 1e6 * (1 + i % 3))).ToArray();
        var centroider = new ProfileCentroider(models, 6, 12, new ConstantPeakWidth(Sigma));
        double[] times = Enumerable.Range(0, 16).Select(i => Rt - 0.4 + 0.05 * i).ToArray();

        var parallel = centroider.Centroid(times);
        for (int s = 0; s < times.Length; s++)
        {
            var serial = centroider.Centroid(times[s]);
            Assert.That(parallel[s].Mz, Is.EqualTo(serial.Mz));
            Assert.That(parallel[s].Intensity, Is.EqualTo(serial.Intensity));
        }
    }
}
