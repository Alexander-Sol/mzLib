using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    // <summary>
    /// Used with adding low frequency noise to objects derived from the SimulatedData class. 
    /// </summary>
    /// <see cref="SimulatedData"/>
    /// <see cref="SimulatedChargeStateEnvelope"/>
    public class LowFrequencyNoiseParameters
    {
        public int PeakNumberLimitLow { get; }
        public int PeakNumberLimitHigh { get; }
        public double PeakLocationLimitLow { get; }
        public double PeakLocationLimitHigh { get; }
        public double PeakWidthLimitLow { get; }
        public double PeakWidthLimitHigh { get; }
        public double PeakIntensityLimitLow { get; }
        public double PeakIntensityLimitHigh { get; }
        public (double Min, double Max)? ExcludedZone { get; }
        public PeakShapeOptions PeakShapeOptions { get; }
        public LowFrequencyNoiseParameters(
            int peakNumberLimitLow = 100, 
            int peakNumberLimitHigh = 1000,
            double peakLocationLimitLow = 200, 
            double peakLocationLimitHigh = 2000,
            double peakWidthLimitLow = 0.01,
            double peakWidthLimitHigh = 0.02,
            double peakIntensityLimitLow = 5, 
            double peakIntensityLimitHigh = 250,
            (double Min, double Max)? excludedZone = null,
            PeakShapeOptions peakShapeOptions = PeakShapeOptions.Centroid)
        {
            PeakNumberLimitLow = peakNumberLimitLow;
            PeakNumberLimitHigh = peakNumberLimitHigh;
            PeakLocationLimitLow = peakLocationLimitLow;
            PeakLocationLimitHigh = peakLocationLimitHigh;
            PeakWidthLimitLow = peakWidthLimitLow;
            PeakWidthLimitHigh = peakWidthLimitHigh;
            PeakIntensityLimitLow = peakIntensityLimitLow;
            PeakIntensityLimitHigh = peakIntensityLimitHigh;
            ExcludedZone = excludedZone;
            PeakShapeOptions = peakShapeOptions;
        }
    }

    public enum PeakShapeOptions
    {
        Gaussian,
        Centroid,
        None
    }
}
