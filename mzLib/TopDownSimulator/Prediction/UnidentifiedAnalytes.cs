#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using TopDownSimulator.Model;

namespace TopDownSimulator.Prediction;

/// <summary>
/// Synthetic proteoforms standing in for the analytes a run contains but its search did not
/// identify. They are rendered into the spectra and left out of the feature truth
/// (<c>Simulator.WriteMzml</c>'s <c>background</c>).
/// </summary>
/// <remarks>
/// <para>
/// Each one is a fitted model from the training run with its mass, elution time and abundance
/// moved, so the population keeps the joint shape of what was fitted: mass with charge, elution
/// width with tail. Its charge centre is rescaled with √mass, as the fitted μ_z scales, and its
/// width follows from the centre the way the fitted ones do.
/// </para>
/// <para>
/// In rep1 fract7 only about 15 % of the peaks at S/N ≥ 10 between 30 and 50 min belong to an
/// identified species. Without this component a simulation has nothing to put there except the tail
/// of its noise.
/// </para>
/// </remarks>
/// <param name="Count">How many to draw.</param>
/// <param name="LogAbundanceShift">Added to each drawn model's log10 abundance: unidentified analytes are mostly fainter.</param>
/// <param name="LogAbundanceSd">Scatter added to the log10 abundance.</param>
/// <param name="RtSd">Standard deviation of the elution-time shift, minutes.</param>
/// <param name="RelativeMassSd">Standard deviation of the relative mass change.</param>
public sealed record UnidentifiedAnalytes(
    int Count,
    double LogAbundanceShift,
    double LogAbundanceSd = 0.5,
    double RtSd = 1.0,
    double RelativeMassSd = 0.15)
{
    /// <summary>
    /// Draws <see cref="Count"/> models from <paramref name="fitted"/>, which must carry Gaussian
    /// charge distributions. Retention times are kept inside [<paramref name="minRt"/>, <paramref name="maxRt"/>].
    /// </summary>
    public ProteoformModel[] Draw(IReadOnlyList<ProteoformModel> fitted, double minRt, double maxRt, int seed = 0)
    {
        if (fitted is null) throw new ArgumentNullException(nameof(fitted));
        var pool = fitted.Where(m => m.Abundance > 0 && m.ChargeDistribution is GaussianChargeDistribution).ToArray();
        if (pool.Length == 0) throw new ArgumentException("No fitted models with a Gaussian charge distribution to draw from.", nameof(fitted));

        var rng = new Random(seed);
        double Normal() => Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());

        var drawn = new ProteoformModel[Count];
        for (int i = 0; i < Count; i++)
        {
            var source = pool[rng.Next(pool.Length)];
            var charge = (GaussianChargeDistribution)source.ChargeDistribution;

            double mass = source.MonoisotopicMass * Math.Exp(RelativeMassSd * Normal());
            double scale = Math.Sqrt(mass / source.MonoisotopicMass);
            double rt = Math.Clamp(source.RtProfile.Mu + RtSd * Normal(), minRt, maxRt);
            double logAbundance = Math.Log10(source.Abundance) + LogAbundanceShift + LogAbundanceSd * Normal();

            drawn[i] = new ProteoformModel(
                mass,
                Math.Pow(10, logAbundance),
                source.RtProfile with { Mu = rt },
                new GaussianChargeDistribution(charge.MuZ * scale, charge.SigmaZ * scale),
                $"unidentified:{i}");
        }

        return drawn;
    }
}
