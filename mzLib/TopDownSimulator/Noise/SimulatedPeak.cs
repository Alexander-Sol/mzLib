#nullable enable
namespace TopDownSimulator.Noise;

/// <summary>Where a peak in a simulated scan came from.</summary>
public enum PeakOrigin
{
    /// <summary>A real proteoform isotopologue produced by the forward model.</summary>
    Signal = 0,

    /// <summary>An unstructured detector/FT noise peak with no isotopic partner.</summary>
    FtNoise = 1,
}

/// <summary>
/// One peak on its way into a simulated scan, tagged with its provenance.
/// </summary>
/// <remarks>
/// The tag is what lets a benchmark distinguish "the feature finder recovered a real proteoform"
/// from "it recovered a contaminant" from "it fired on nothing", which is the whole point of
/// simulating noise rather than just adding it.
/// </remarks>
public readonly record struct SimulatedPeak(double Mz, double Intensity, PeakOrigin Origin);
