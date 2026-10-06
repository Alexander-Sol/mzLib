#nullable enable
using System;

namespace TopDownSimulator.Noise;

/// <summary>
/// The signal-to-noise ratio a reported noise peak carries, as a shifted exponential:
/// S/N = <see cref="Threshold"/> + Exponential(<see cref="MeanExcess"/>).
/// </summary>
/// <remarks>
/// <para>
/// The shift is the instrument's detection threshold — nothing below it is reported at all — and
/// the exponential tail is what the amplitude of a noise field above a high threshold looks like.
/// </para>
/// <para>
/// This distribution turned out to be essentially universal. Measured across four Jurkat runs, in
/// both the busy 31–35 min window and the quiet 50–55 min one, the S/N quantiles agreed to within a
/// few percent everywhere below p95 (p1 ≈ 1.37, p25 ≈ 2.09, p50 ≈ 2.75, p75 ≈ 3.90). They diverge
/// only at p99, and only in the busy window, which is real analyte rather than noise. That is why
/// <see cref="Orbitrap"/> is a fixed default while the noise <i>amplitude</i> it multiplies is a
/// per-file parameter.
/// </para>
/// </remarks>
public sealed record SignalToNoiseDistribution(double Threshold, double MeanExcess)
{
    /// <summary>
    /// Fitted to the pooled quiet-window S/N of the Jurkat runs through its median and p95, which
    /// reproduces p25 and p75 to under 2 % and p99 to about 5 %.
    /// </summary>
    /// <remarks>
    /// The hard floor at 1.565 sits slightly above the measured p1 of ~1.37, so roughly the faintest
    /// 1 % of real noise peaks are not representable. That is a deliberate trade: lowering the
    /// threshold to 1.35 and refitting degrades the fit through the body of the distribution, and
    /// the body is what a feature finder actually contends with.
    /// </remarks>
    public static readonly SignalToNoiseDistribution Orbitrap = new(Threshold: 1.565, MeanExcess: 1.707);

    public double Threshold { get; init; } = Threshold > 0 && double.IsFinite(Threshold)
        ? Threshold
        : throw new ArgumentOutOfRangeException(nameof(Threshold), Threshold, "Threshold must be finite and positive.");

    public double MeanExcess { get; init; } = MeanExcess > 0 && double.IsFinite(MeanExcess)
        ? MeanExcess
        : throw new ArgumentOutOfRangeException(nameof(MeanExcess), MeanExcess, "Mean excess must be finite and positive.");

    /// <summary>Mean S/N of a reported noise peak.</summary>
    public double Mean => Threshold + MeanExcess;

    public double Sample(Random rng)
    {
        // 1 - NextDouble() excludes 1.0 rather than 0.0, keeping the log finite.
        double u = 1.0 - rng.NextDouble();
        return Threshold - MeanExcess * Math.Log(u);
    }

    /// <summary>The S/N below which a fraction <paramref name="q"/> of noise peaks fall.</summary>
    public double Quantile(double q)
    {
        if (q < 0 || q >= 1)
            throw new ArgumentOutOfRangeException(nameof(q), q, "Quantile must be in [0, 1).");

        return Threshold - MeanExcess * Math.Log(1 - q);
    }

    public override string ToString() => $"S/N ~ {Threshold:R} + Exp(mean={MeanExcess:R})";
}
