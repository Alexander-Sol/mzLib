#nullable enable
using System;
using TopDownSimulator.Simulation;

namespace TopDownSimulator.Noise;

/// <summary>
/// A detection limit set at a multiple of the local noise amplitude — the rule a real instrument
/// actually applies.
/// </summary>
/// <remarks>
/// Using this instead of a fraction of the run's brightest peak is what keeps a noisy simulation
/// self-consistent. The default floor is 10⁻⁴ of the global maximum, which for a realistic run sits
/// far above the noise amplitude; combined with an injected noise floor it would produce a file
/// containing noise peaks fainter than any surviving signal peak, which is not something an
/// instrument can produce.
/// <para>
/// <see cref="SnThreshold"/> defaults to the same threshold the measured noise peaks respect, so
/// signal and noise are reported under one rule.
/// </para>
/// </remarks>
public sealed class NoiseRelativeFloor : IIntensityFloor
{
    private readonly NoiseFloorModel _model;

    public NoiseRelativeFloor(NoiseFloorModel model, double? snThreshold = null)
    {
        _model = model ?? throw new ArgumentNullException(nameof(model));
        double threshold = snThreshold ?? model.SignalToNoise.Threshold;
        if (!(threshold > 0) || !double.IsFinite(threshold))
            throw new ArgumentOutOfRangeException(nameof(snThreshold), threshold,
                "S/N threshold must be finite and positive.");

        SnThreshold = threshold;
    }

    /// <summary>Multiple of the local noise amplitude a peak must reach to be written.</summary>
    public double SnThreshold { get; }

    public double At(double mz) => SnThreshold * _model.NoiseLevelAt(mz);

    /// <summary>
    /// The noise level is not monotonic in m/z — it rises to a maximum near m/z 1450 and falls away
    /// again — but it is unimodal, and the minimum of a unimodal function over an interval is always
    /// at an endpoint. If the calibration table is ever replaced by one that rises and falls more
    /// than once, this has to scan the knots instead.
    /// </summary>
    public double MinOver(double minMz, double maxMz) => Math.Min(At(minMz), At(maxMz));

    public override string ToString() => $"NoiseRelative(S/N>{SnThreshold:R})";
}
