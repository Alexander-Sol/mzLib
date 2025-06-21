using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    /// <summary>
    /// Algorithm for modelling chromatographic peak shapes based on this paper
    /// https://pubs.acs.org/doi/10.1021/pr400727e
    /// </summary>
    public class PeakShapeAlgorithm
    {
        public double TailingFactor { get; set; } = 0.5;
        public double FrontingFactor { get; set; } = 0.5;   
        public double HeightFactor { get; set; } = 1.0;

        public PeakShapeAlgorithm(double frontingFactor = 0.5, double tailingFactor = 0.5, double heightFactor = 1.0)
        {
            TailingFactor = tailingFactor;
            FrontingFactor = frontingFactor;
            HeightFactor = heightFactor;
        }   

        public List<double> GetIntensityRange(double apexRt, double relativeStartRt, double relativeEndRt, int numPoints = 100)
        {
            List<double> intensities = new List<double>();
            double step = (relativeEndRt - relativeStartRt) / numPoints;
            for (int i = 0; i <= numPoints; i++)
            {
                double relativeRt = apexRt - (relativeStartRt + i * step);
                double intensity = GetIntensity(apexRt, relativeRt);
                intensities.Add(intensity);
            }
            return intensities;
        }

        public List<double> GetIntensityRange(double apexRt, List<double> retentionTimes)
        {
            List<double> intensities = new List<double>();
            double rtStart = retentionTimes.Min();
            apexRt = apexRt - rtStart;
            foreach (var rt in retentionTimes)
            {
                double relativeRt = rt - rtStart;
                double intensity = GetIntensity(apexRt, relativeRt);
                intensities.Add(intensity);
            }
            return intensities;
        }

        public double GetIntensity(double apexRt, double relativeRt)
        {
            // Calculate the peak shape based on the relative retention time
            double rtDiff = Math.Abs(relativeRt - apexRt);
            double elutionStandardDeviation = TailingFactor * relativeRt + FrontingFactor;
            double intensity = HeightFactor * Math.Exp(- Math.Pow(rtDiff / elutionStandardDeviation , 2) / 2);
            return intensity;
        }   

    }
}
