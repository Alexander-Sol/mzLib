#nullable enable
using System;

namespace TopDownSimulator.Noise;

/// <summary>
/// Scan-to-scan measurement error on a real peak: how far its centroid wanders in m/z and how much
/// its intensity varies between consecutive scans.
/// </summary>
/// <remarks>
/// <para>
/// Without this the simulator's XICs are exactly the analytical EMG evaluated at each scan time —
/// perfectly smooth, perfectly reproducible, and far easier than anything an instrument produces.
/// Scan-to-scan scatter is what XIC smoothing, apex detection and ppm-tolerance matching actually
/// contend with, so its absence was the largest remaining gap between simulated and real data.
/// </para>
/// <para>
/// Both laws were measured on rep2 fract7 against the top 400 MetaMorpheus identifications, using
/// estimators that cancel the elution profile rather than model it (see
/// <c>Test/TopDownSimulator/SignalJitterCharacterization.cs</c>):
/// </para>
/// <list type="bullet">
/// <item>m/z error was split into the part shared by every peak in a scan — calibration drift,
/// σ = 1.89 ppm — and the per-peak remainder, which follows σ = √(9.06²/(S/N) + 0.81²) ppm across
/// S/N from 2 to 128 (measured 5.45, 3.84, 2.77, 1.78, 1.58, 1.14 ppm; the law reproduces all six to
/// within 15 %).</item>
/// <item>Intensity error came from the scan-to-scan scatter of the ratio between adjacent
/// isotopologues of one species, which is fixed by chemistry so the elution profile divides out
/// exactly. Per-peak σ of log intensity follows √(1.026²/(S/N) + 0.165²), reproducing the six
/// measured bins to within 6 %.</item>
/// </list>
/// <para>
/// The <b>0.165 floor on log-intensity σ is the headline number</b>: even an S/N of several hundred
/// still scatters ~17 % scan to scan, where the simulator currently produces 0 %.
/// </para>
/// </remarks>
public sealed record PeakJitterModel
{
    /// <summary>The measured Orbitrap defaults; every figure is documented on the type.</summary>
    public static readonly PeakJitterModel Orbitrap = new();

    /// <summary>Jitter switched off, for isolating its effect in a comparison.</summary>
    public static readonly PeakJitterModel None = new()
    {
        MzPpmAtUnitSignalToNoise = 0,
        MzPpmFloor = 0,
        CommonModeMzPpm = 0,
        LogIntensitySigmaAtUnitSignalToNoise = 0,
        LogIntensitySigmaFloor = 0,
        DropPeaksFallingBelowFloor = false,
    };

    /// <summary>Coefficient of the 1/√(S/N) term of the per-peak m/z error, in ppm.</summary>
    public double MzPpmAtUnitSignalToNoise { get; init; } = 9.06;

    /// <summary>Irreducible per-peak m/z error at high S/N, in ppm.</summary>
    public double MzPpmFloor { get; init; } = 0.81;

    /// <summary>
    /// σ of the per-scan m/z offset shared by every peak in that scan — mass calibration drift.
    /// </summary>
    /// <remarks>
    /// Kept separate from the per-peak term because the two stress a feature finder differently: a
    /// common-mode shift moves a whole spectrum together and can be calibrated out, while
    /// independent per-peak error cannot. Measured at 1.89 ppm against a per-peak 3.98 ppm, so the
    /// per-peak term dominates but neither is negligible.
    /// </remarks>
    public double CommonModeMzPpm { get; init; } = 1.887;

    /// <summary>Coefficient of the 1/√(S/N) term of the log-intensity σ.</summary>
    public double LogIntensitySigmaAtUnitSignalToNoise { get; init; } = 1.026;

    /// <summary>Irreducible log-intensity σ at high S/N — about 17 %.</summary>
    public double LogIntensitySigmaFloor { get; init; } = 0.165;

    /// <summary>
    /// σ of a per-scan multiplicative intensity factor common to every peak, in log space.
    /// </summary>
    /// <remarks>
    /// Defaults to zero because it is <b>unmeasured</b>, not because it is believed absent. The
    /// isotopologue-ratio estimator divides out anything affecting both peaks equally, so it is
    /// blind to exactly this term; a real AGC fluctuation would produce one. Set it if you have an
    /// independent estimate.
    /// </remarks>
    public double CommonModeLogIntensitySigma { get; init; }

    /// <summary>
    /// Whether a peak whose jittered intensity falls under the detection floor is dropped.
    /// </summary>
    /// <remarks>
    /// On by default because it is what an instrument does, and because the resulting XIC dropout is
    /// realistic and worth simulating: a peak sitting just above the floor has S/N ≈ 1.6 and hence
    /// σ_log ≈ 0.84, so it genuinely disappears in a good fraction of scans. Note the consequence
    /// for <see cref="Simulation.FeatureGroundTruth"/> — a feature's first and last scan numbers
    /// become approximate under jitter, though its identity, charge, mass and apex do not.
    /// </remarks>
    public bool DropPeaksFallingBelowFloor { get; init; } = true;

    /// <summary>Per-peak m/z error in ppm at the given signal-to-noise ratio.</summary>
    public double MzPpmSigma(double signalToNoise)
    {
        double a = MzPpmAtUnitSignalToNoise;
        double b = MzPpmFloor;
        if (!(signalToNoise > 0) || !double.IsFinite(signalToNoise))
            return b;

        return Math.Sqrt(a * a / signalToNoise + b * b);
    }

    /// <summary>σ of log intensity at the given signal-to-noise ratio.</summary>
    public double LogIntensitySigma(double signalToNoise)
    {
        double c = LogIntensitySigmaAtUnitSignalToNoise;
        double d = LogIntensitySigmaFloor;
        if (!(signalToNoise > 0) || !double.IsFinite(signalToNoise))
            return d;

        return Math.Sqrt(c * c / signalToNoise + d * d);
    }

    /// <summary>
    /// A mean-preserving multiplicative intensity factor: <c>exp(N(−σ²/2, σ))</c>, whose expectation
    /// is 1 rather than <c>exp(σ²/2)</c>.
    /// </summary>
    /// <remarks>
    /// At the σ ≈ 0.84 a peak near the detection floor carries, a median-preserving draw would
    /// inflate expected intensity by 43 %, quietly brightening the faintest peaks in the file — the
    /// ones whose detectability the simulation is supposed to be testing.
    /// </remarks>
    public static double IntensityFactor(Random rng, double logSigma)
    {
        if (!(logSigma > 0))
            return 1.0;

        return Math.Exp(-0.5 * logSigma * logSigma + logSigma * NoiseFloorModel.SampleStandardNormal(rng));
    }

    /// <summary>True if any component of this model would change a peak.</summary>
    public bool IsActive =>
        MzPpmAtUnitSignalToNoise > 0 || MzPpmFloor > 0 || CommonModeMzPpm > 0 ||
        LogIntensitySigmaAtUnitSignalToNoise > 0 || LogIntensitySigmaFloor > 0 ||
        CommonModeLogIntensitySigma > 0;
}
