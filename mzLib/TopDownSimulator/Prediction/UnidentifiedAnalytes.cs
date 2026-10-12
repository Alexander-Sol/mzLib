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
/// <param name="FollowIdentified">
/// Fraction of the drawn models whose elution time is taken from the simulated run's identified
/// load instead of their training model's. 0 keeps the training run's RT structure throughout.
/// </param>
/// <remarks>
/// Across rep2 fract2–7 (<c>NoiseAcrossRuns</c>) the ion load per minute departs from the
/// fractions' shared RT curve by about 0.45× the departure of the identified precursor intensity, in
/// log10. Half the unidentified load is therefore placed on the training run's curve and half where
/// the simulated run's identifications are.
/// </remarks>
public sealed record UnidentifiedAnalytes(
    int Count,
    double LogAbundanceShift,
    double LogAbundanceSd = 0.5,
    double RtSd = 1.0,
    double RelativeMassSd = 0.15,
    double FollowIdentified = 0)
{
    /// <summary>
    /// Fewer identifications than this say too little about a run's load to scale it.
    /// </summary>
    /// <remarks>
    /// The integrated TIC over the summed identified precursor intensity is 10^1.35 ± 0.04 in rep2
    /// fract2–7 (490–1592 identifications). In fract8, 9, 10 and 1 (241, 116, 22, 7) it is
    /// 10^2.1–3.1, with nothing to say which. Scaling fract9 by its identifications put the simulated
    /// TIC 0.63 dex low, against 0.29 dex high without scaling.
    /// </remarks>
    public const int MinimumIdentificationsForLoad = 250;

    /// <summary>
    /// The factor that scales the unidentified load from the calibration run's to the simulated
    /// run's: the ratio of their summed identified precursor intensities, which must come from the
    /// same search. Null when either run has fewer than <see cref="MinimumIdentificationsForLoad"/>
    /// identifications with a precursor intensity; the calibration run's load is then the best guess.
    /// </summary>
    public static double? LoadScale(IEnumerable<double> identifiedIntensities, IEnumerable<double> referenceIntensities)
    {
        var run = identifiedIntensities.Where(v => v > 0 && double.IsFinite(v)).ToArray();
        var reference = referenceIntensities.Where(v => v > 0 && double.IsFinite(v)).ToArray();
        if (run.Length < MinimumIdentificationsForLoad || reference.Length < MinimumIdentificationsForLoad) return null;
        return run.Sum() / reference.Sum();
    }

    /// <summary>
    /// Draws <see cref="Count"/> models from <paramref name="fitted"/>, which must carry Gaussian
    /// charge distributions. Retention times are kept inside [<paramref name="minRt"/>, <paramref name="maxRt"/>].
    /// </summary>
    /// <param name="identifiedLoad">
    /// The simulated run's identifications as (retention time, precursor intensity). A fraction
    /// <see cref="FollowIdentified"/> of the models is centred on retention times drawn from it in
    /// proportion to intensity. Ignored when null or <see cref="FollowIdentified"/> is 0.
    /// </param>
    /// <param name="loadScale">
    /// Multiplies every drawn abundance. Over rep2 fract2–7 a run's integrated TIC is
    /// 10^1.35 ± 0.04 times the summed precursor intensity of its identifications, so the ratio of
    /// that sum to the calibration run's, from the same search, scales the unidentified load.
    /// </param>
    public ProteoformModel[] Draw(
        IReadOnlyList<ProteoformModel> fitted, double minRt, double maxRt, int seed = 0,
        IReadOnlyList<(double Rt, double Intensity)>? identifiedLoad = null, double loadScale = 1)
    {
        if (fitted is null) throw new ArgumentNullException(nameof(fitted));
        if (!(loadScale > 0) || !double.IsFinite(loadScale))
            throw new ArgumentOutOfRangeException(nameof(loadScale), loadScale, "Must be finite and positive.");
        var pool = fitted.Where(m => m.Abundance > 0 && m.ChargeDistribution is GaussianChargeDistribution).ToArray();
        if (pool.Length == 0) throw new ArgumentException("No fitted models with a Gaussian charge distribution to draw from.", nameof(fitted));

        var load = FollowIdentified > 0 && identifiedLoad is not null
            ? identifiedLoad.Where(p => p.Intensity > 0 && double.IsFinite(p.Intensity) && double.IsFinite(p.Rt)).ToArray()
            : Array.Empty<(double Rt, double Intensity)>();
        var cumulative = new double[load.Length];
        for (int i = 0; i < load.Length; i++)
            cumulative[i] = (i == 0 ? 0 : cumulative[i - 1]) + load[i].Intensity;

        var rng = new Random(seed);
        // A second stream, so that FollowIdentified = 0 reproduces the earlier draws exactly.
        var followRng = new Random(seed + 7919);
        double Normal() => Math.Sqrt(-2 * Math.Log(1 - rng.NextDouble())) * Math.Cos(2 * Math.PI * rng.NextDouble());

        var drawn = new ProteoformModel[Count];
        for (int i = 0; i < Count; i++)
        {
            var source = pool[rng.Next(pool.Length)];
            var charge = (GaussianChargeDistribution)source.ChargeDistribution;

            double mass = source.MonoisotopicMass * Math.Exp(RelativeMassSd * Normal());
            double scale = Math.Sqrt(mass / source.MonoisotopicMass);
            double centre = source.RtProfile.Mu;
            if (load.Length > 0 && followRng.NextDouble() < FollowIdentified)
            {
                int at = Array.BinarySearch(cumulative, followRng.NextDouble() * cumulative[^1]);
                centre = load[Math.Min(load.Length - 1, at < 0 ? ~at : at)].Rt;
            }

            double rt = Math.Clamp(centre + RtSd * Normal(), minRt, maxRt);
            double logAbundance = Math.Log10(source.Abundance * loadScale) + LogAbundanceShift + LogAbundanceSd * Normal();

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
