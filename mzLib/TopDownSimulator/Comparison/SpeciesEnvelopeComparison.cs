#nullable enable
using System;
using TopDownSimulator.Extraction;

namespace TopDownSimulator.Comparison;

/// <summary>
/// How closely one species' simulated envelope matches the real one around its apex: isotopologue
/// intensities per charge, summed over the scans of the extraction window so that runs on different
/// scan grids compare directly.
/// </summary>
/// <param name="Cosine">Cosine over every (charge, isotopologue) intensity. Judges the envelope's shape.</param>
/// <param name="ChargeProfileCosine">Cosine over the per-charge sums. Judges the charge distribution alone.</param>
/// <param name="IntensityRatio">Simulated over real summed intensity. Judges the abundance.</param>
/// <param name="RealApexCharge">Most intense charge in the real envelope.</param>
/// <param name="SimulatedApexCharge">Most intense charge in the simulated envelope.</param>
public sealed record SpeciesEnvelopeComparison(
    double Cosine,
    double ChargeProfileCosine,
    double IntensityRatio,
    int RealApexCharge,
    int SimulatedApexCharge)
{
    /// <summary>
    /// Compares two extractions of the same species over the same charge range, one from the real
    /// run and one from the simulation. Returns null when either is empty.
    /// </summary>
    /// <remarks>
    /// Both extractions see whatever else lands in their windows: noise, and in the real run,
    /// unidentified species. The simulation's noise is matched to the real run's, so that floor is
    /// shared; unidentified signal is not, and pulls the real envelope's shape and sum up.
    /// </remarks>
    public static SpeciesEnvelopeComparison? Compare(ProteoformGroundTruth real, ProteoformGroundTruth simulated)
    {
        ArgumentNullException.ThrowIfNull(real);
        ArgumentNullException.ThrowIfNull(simulated);
        if (real.MinCharge != simulated.MinCharge || real.MaxCharge != simulated.MaxCharge)
            throw new ArgumentException("Both extractions must cover the same charges.");

        int charges = real.ChargeCount;
        double dot = 0, rNorm = 0, sNorm = 0, rSum = 0, sSum = 0;
        var rProfile = new double[charges];
        var sProfile = new double[charges];
        for (int c = 0; c < charges; c++)
        {
            int isotopologues = Math.Min(real.IsotopologueIntensities[c].Length, simulated.IsotopologueIntensities[c].Length);
            for (int k = 0; k < isotopologues; k++)
            {
                double r = Sum(real.IsotopologueIntensities[c][k]);
                double s = Sum(simulated.IsotopologueIntensities[c][k]);
                dot += r * s;
                rNorm += r * r;
                sNorm += s * s;
                rProfile[c] += r;
                sProfile[c] += s;
            }

            rSum += rProfile[c];
            sSum += sProfile[c];
        }

        if (!(rSum > 0) || !(sSum > 0)) return null;

        double profileDot = 0, rpNorm = 0, spNorm = 0;
        for (int c = 0; c < charges; c++)
        {
            profileDot += rProfile[c] * sProfile[c];
            rpNorm += rProfile[c] * rProfile[c];
            spNorm += sProfile[c] * sProfile[c];
        }

        return new SpeciesEnvelopeComparison(
            dot / Math.Sqrt(rNorm * sNorm),
            profileDot / Math.Sqrt(rpNorm * spNorm),
            sSum / rSum,
            real.MinCharge + ArgMax(rProfile),
            real.MinCharge + ArgMax(sProfile));
    }

    private static double Sum(double[] values)
    {
        double total = 0;
        foreach (double v in values) total += v;
        return total;
    }

    private static int ArgMax(double[] values)
    {
        int best = 0;
        for (int i = 1; i < values.Length; i++)
            if (values[i] > values[best]) best = i;
        return best;
    }
}
