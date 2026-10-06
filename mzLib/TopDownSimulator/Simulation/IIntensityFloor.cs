#nullable enable
using System;

namespace TopDownSimulator.Simulation;

/// <summary>
/// The intensity below which a simulated peak is not written to file.
/// </summary>
/// <remarks>
/// This is an m/z-dependent quantity rather than a scalar because a real instrument's detection
/// limit is one: it reports what clears a multiple of the local noise, and the noise amplitude in an
/// Orbitrap varies about 4.6× between m/z 650 and m/z 1500. Holding the floor constant across the
/// range therefore over-reports at low m/z and under-reports at high m/z.
/// <para>
/// The same object has to serve reduction and <see cref="FeatureGroundTruth"/>, since the truth is
/// only correct if it describes exactly what survived reduction.
/// </para>
/// </remarks>
public interface IIntensityFloor
{
    /// <summary>The floor at <paramref name="mz"/>.</summary>
    double At(double mz);

    /// <summary>
    /// A lower bound on the floor over [<paramref name="minMz"/>, <paramref name="maxMz"/>]. Callers
    /// use it to discard whole envelopes cheaply; over-estimating it would drop features that are
    /// genuinely in the file.
    /// </summary>
    double MinOver(double minMz, double maxMz);
}

/// <summary>A floor that does not vary with m/z.</summary>
public sealed record ConstantIntensityFloor(double Value) : IIntensityFloor
{
    public double Value { get; init; } = Value >= 0 && double.IsFinite(Value)
        ? Value
        : throw new ArgumentOutOfRangeException(nameof(Value), Value, "Floor must be finite and non-negative.");

    public double At(double mz) => Value;

    public double MinOver(double minMz, double maxMz) => Value;

    public override string ToString() => $"Constant({Value:R})";
}
