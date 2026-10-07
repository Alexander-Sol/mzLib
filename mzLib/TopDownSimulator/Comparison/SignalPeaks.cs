#nullable enable
using System;
using System.Collections.Generic;
using TopDownSimulator.Noise;

namespace TopDownSimulator.Comparison;

/// <summary>
/// The peaks of a centroid spectrum that stand clear of the noise, so that comparisons can judge
/// the signal model without the noise floor dominating them.
/// </summary>
/// <remarks>
/// The noise level comes from a <see cref="NoiseFloorModel"/>, normally one scan's
/// (<see cref="ScanNoiseConditions.FromSourceScans"/>: amplitude from injection time, shape in m/z
/// from the Jurkat measurement). To compare a simulated scan with a real one, select both with the
/// real scan's model so the threshold is the same. The median centroid intensity is not a noise
/// level: the instrument reports only local maxima above its own threshold, so the median sits well
/// above the noise and S/N 10 against it keeps a handful of peaks per scan.
/// </remarks>
public static class SignalPeaks
{
    /// <summary>
    /// The peaks in [<paramref name="minMz"/>, <paramref name="maxMz"/>] with intensity at least
    /// <paramref name="minSignalToNoise"/> times <paramref name="noise"/>'s level at their m/z, in
    /// their original order.
    /// </summary>
    public static (double[] Mz, double[] Intensity) Select(
        double[] mz, double[] intensity, NoiseFloorModel noise, double minSignalToNoise = 10,
        double minMz = double.NegativeInfinity, double maxMz = double.PositiveInfinity)
    {
        ArgumentNullException.ThrowIfNull(mz);
        ArgumentNullException.ThrowIfNull(intensity);
        ArgumentNullException.ThrowIfNull(noise);
        if (mz.Length != intensity.Length) throw new ArgumentException("The m/z and intensity arrays differ in length.");

        var keptMz = new List<double>();
        var keptIntensity = new List<double>();
        for (int i = 0; i < mz.Length; i++)
        {
            if (mz[i] < minMz || mz[i] > maxMz || intensity[i] < minSignalToNoise * noise.NoiseLevelAt(mz[i])) continue;
            keptMz.Add(mz[i]);
            keptIntensity.Add(intensity[i]);
        }

        return (keptMz.ToArray(), keptIntensity.ToArray());
    }
}
