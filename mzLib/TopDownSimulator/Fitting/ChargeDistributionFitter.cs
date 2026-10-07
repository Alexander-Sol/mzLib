using System;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace TopDownSimulator.Fitting;

public sealed record FittedChargeDistribution(
    GaussianChargeDistribution Distribution,
    double[] ApexIntensitiesByCharge,
    int ChargesUsed);

/// <summary>
/// Fits a Gaussian-on-charge distribution (μ_z, σ_z) to per-charge XIC apex intensities
/// using a weighted method-of-moments estimator. When fewer than two charges carry
/// signal, σ_z falls back to a user-supplied default (1.0 by default).
/// </summary>
/// <remarks>
/// <para>
/// Moments over every extracted charge are only right when the window is about as wide as the
/// envelope. In a wider window the outer charges carry noise and other species' peaks, which pull
/// μ_z towards the window centre and inflate σ_z. In a narrower one the envelope is truncated and
/// μ_z is pulled towards the window centre the other way.
/// </para>
/// <para>
/// A positive <c>trimFraction</c> makes the fit independent of the window: it keeps only the
/// contiguous run of charges around the most intense one whose apex is at least that fraction of
/// the maximum, ending early where the apex climbs again into another species' envelope, and fits a parabola to log apex over the run, weighted by apex² (Caruana's
/// estimator), which is exact for a Gaussian however it is truncated. It falls back to moments over
/// the run when the run is shorter than three charges or the parabola does not open downwards.
/// </para>
/// </remarks>
public sealed class ChargeDistributionFitter
{
    private readonly double _fallbackSigmaZ;
    private readonly double _trimFraction;

    /// <param name="fallbackSigmaZ">σ_z when fewer than two charges carry signal.</param>
    /// <param name="trimFraction">
    /// 0 for moments over every extracted charge; otherwise the fraction of the most intense charge's
    /// apex below which the run of charges used for the fit ends.
    /// </param>
    public ChargeDistributionFitter(double fallbackSigmaZ = 1.0, double trimFraction = 0)
    {
        if (fallbackSigmaZ <= 0) throw new ArgumentOutOfRangeException(nameof(fallbackSigmaZ));
        if (trimFraction < 0 || trimFraction >= 1) throw new ArgumentOutOfRangeException(nameof(trimFraction));
        _fallbackSigmaZ = fallbackSigmaZ;
        _trimFraction = trimFraction;
    }

    public FittedChargeDistribution Fit(ProteoformGroundTruth truth)
    {
        int nCharges = truth.ChargeCount;
        int nScans = truth.ScanCount;

        // Per-charge apex intensity from its XIC.
        double[] apex = new double[nCharges];
        int chargesUsed = 0;
        for (int c = 0; c < nCharges; c++)
        {
            double best = 0;
            for (int s = 0; s < nScans; s++)
            {
                if (truth.ChargeXics[c][s] > best) best = truth.ChargeXics[c][s];
            }
            apex[c] = best;
            if (best > 0) chargesUsed++;
        }

        if (chargesUsed == 0)
            throw new InvalidOperationException("Cannot fit charge distribution: all charge XICs are zero.");

        int first = 0, last = nCharges - 1;
        if (_trimFraction > 0)
        {
            int top = 0;
            for (int c = 1; c < nCharges; c++)
                if (apex[c] > apex[top]) top = c;
            first = RunEnd(apex, top, -1);
            last = RunEnd(apex, top, +1);

            if (TryFitLogParabola(apex, first, last, truth.MinCharge, out double mu, out double sigma))
                return new FittedChargeDistribution(new GaussianChargeDistribution(mu, sigma), apex, last - first + 1);
            chargesUsed = last - first + 1;
        }

        double sumW = 0, sumWz = 0;
        for (int c = first; c <= last; c++)
        {
            int z = truth.MinCharge + c;
            sumW += apex[c];
            sumWz += apex[c] * z;
        }
        double muZ = sumWz / sumW;

        double sigmaZ;
        if (chargesUsed < 2)
        {
            sigmaZ = _fallbackSigmaZ;
        }
        else
        {
            double sumWdd = 0;
            for (int c = first; c <= last; c++)
            {
                int z = truth.MinCharge + c;
                double d = z - muZ;
                sumWdd += apex[c] * d * d;
            }
            double var = sumWdd / sumW;
            sigmaZ = var > 0 ? Math.Sqrt(var) : _fallbackSigmaZ;
        }

        return new FittedChargeDistribution(
            new GaussianChargeDistribution(muZ, sigmaZ),
            apex,
            chargesUsed);
    }

    /// <summary>
    /// The last charge, walking from <paramref name="top"/> in direction <paramref name="step"/>, that
    /// still belongs to the envelope: the walk stops at a zero, below the trim floor, or where the
    /// apex climbs back to more than 1.5× the lowest value passed once it has fallen below half the
    /// maximum, which is another species' envelope rather than this one's flank.
    /// </summary>
    private int RunEnd(double[] apex, int top, int step)
    {
        double floor = _trimFraction * apex[top];
        double lowest = apex[top];
        int end = top;
        for (int c = top + step; c >= 0 && c < apex.Length; c += step)
        {
            double a = apex[c];
            if (a <= 0 || a < floor) break;
            if (lowest < 0.5 * apex[top] && a > 1.5 * lowest) break;
            lowest = Math.Min(lowest, a);
            end = c;
        }

        return end;
    }

    /// <summary>
    /// Caruana's estimator: ln y = a + b·z + c·z² by least squares weighted by y², which is exact for
    /// noiseless Gaussian samples and keeps the low-intensity ends from dominating through the log.
    /// </summary>
    private static bool TryFitLogParabola(double[] apex, int first, int last, int minCharge, out double mu, out double sigma)
    {
        mu = sigma = double.NaN;
        if (last - first < 2) return false;

        // Centre z on the run so the normal equations stay well conditioned.
        double z0 = minCharge + 0.5 * (first + last);
        double s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, t0 = 0, t1 = 0, t2 = 0;
        for (int c = first; c <= last; c++)
        {
            double w = apex[c] * apex[c];
            double x = minCharge + c - z0, ly = Math.Log(apex[c]);
            double x2 = x * x;
            s0 += w; s1 += w * x; s2 += w * x2; s3 += w * x2 * x; s4 += w * x2 * x2;
            t0 += w * ly; t1 += w * x * ly; t2 += w * x2 * ly;
        }

        // Solve [[s0 s1 s2] [s1 s2 s3] [s2 s3 s4]] [a b c]ᵀ = [t0 t1 t2]ᵀ by Cramer's rule.
        double det = s0 * (s2 * s4 - s3 * s3) - s1 * (s1 * s4 - s3 * s2) + s2 * (s1 * s3 - s2 * s2);
        if (!(Math.Abs(det) > 1e-12 * s0 * s0 * s0)) return false;
        double b = (s0 * (t1 * s4 - s3 * t2) - t0 * (s1 * s4 - s3 * s2) + s2 * (s1 * t2 - t1 * s2)) / det;
        double cc = (s0 * (s2 * t2 - t1 * s3) - s1 * (s1 * t2 - t1 * s2) + t0 * (s1 * s3 - s2 * s2)) / det;
        if (!(cc < 0)) return false;

        mu = z0 - b / (2 * cc);
        sigma = Math.Sqrt(-1 / (2 * cc));
        // A vertex outside the run means the run is one flank of something; the moments are safer.
        if (mu < minCharge + first - 0.5 || mu > minCharge + last + 0.5) return false;
        return double.IsFinite(mu) && double.IsFinite(sigma);
    }
}
