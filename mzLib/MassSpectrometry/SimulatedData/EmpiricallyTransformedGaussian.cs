using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MathNet.Numerics.Interpolation;
//using Vector = MathNet.Numerics.LinearAlgebra;

namespace MassSpectrometry
{
    public class EmpiricallyTransformedGaussian
    {
        public double Height { get; set; } = 100.0;

        public double FullWidthHalfMax { get; set; } = 0.25;

        //Leading parameters

        // Lambda / PreExponential Parameters. Smaller values = Thicker Peaks
        // Kexp / Exponential Parameters. Smaller values = Thicker Peaks

        public double LambdaLeading { get; set; } = 5;
        public double KexpLeading { get; set; } = 15;
        public double AlphaLeading { get; set; } = 10;

        // Tailing Parameters
        public double LambdaTailing { get; set; } = 5;
        public double KexpTailing { get; set; } = 15;
        public double AlphaTailing { get; set; } = 10;

        public EmpiricallyTransformedGaussian() { // Default constructor with default values
        }

        public EmpiricallyTransformedGaussian(Vector<double> parameters)
        {
            if (parameters.Count != 7)
                throw new ArgumentException("Expected 7 parameters: FullWidthHalfMax, 3 Leading parameters, and 3 Tailing parameters.");
            
            FullWidthHalfMax = parameters[0];
            LambdaLeading = parameters[1];
            KexpLeading = parameters[2];
            AlphaLeading = parameters[3];
            LambdaTailing = parameters[4];
            KexpTailing = parameters[5];
            AlphaTailing = parameters[6];
        }

        public double GetIntensity(double rt, double apexRt)
        {
            double halfMaxRtLeading = apexRt - FullWidthHalfMax / 2;
            double halfMaxRtTailing = apexRt + FullWidthHalfMax / 2;

            // Calculate the peak shape based on the relative retention time
            double leadingTerm = 1 + LambdaLeading * Math.Exp(KexpLeading * (halfMaxRtLeading - rt));
            double leadingTermPower = Math.Pow(halfMaxRtLeading / rt, AlphaLeading);

            double tailingTerm = 1 + LambdaTailing * Math.Exp(KexpTailing * (rt - halfMaxRtTailing));
            double tailingTermPower = Math.Pow(rt / halfMaxRtTailing, AlphaTailing);

            double denominator = Math.Pow(leadingTerm, leadingTermPower) + Math.Pow(tailingTerm, tailingTermPower) - 1;

            return Height / denominator;
        }

        public List<double> GetIntensityRange(double apexRt, List<double> retentionTimes)
        {
            List<double> intensities = new List<double>();
            foreach (var rt in retentionTimes)
            {
                double intensity = GetIntensity(rt, apexRt);
                intensities.Add(intensity);
            }
            return intensities;
        }

        public static Vector<double> DefaultSettings = Vector<double>.Build.Dense(new[]
        {
            0.25, // FullWidthHalfMax
            5.0,  // LambdaLeading
            15.0, // KexpLeading
            10.0, // AlphaLeading
            5.0,  // LambdaTailing
            15.0, // KexpTailing
            10.0  // AlphaTailing
        });

        /// <summary>
        /// This method creates an objective function for optimization based on the indexed peaks.
        /// </summary>
        /// <param name="indexedPeaks"></param>
        /// <returns></returns>
        /// <exception cref="ArgumentException"></exception>
        public static Func<Vector<double>, double> CreateObjectiveFunction(List<IIndexedPeak> indexedPeaks)
        {
            var maxIntensity = indexedPeaks.Max(p => p.Intensity);
            var apexRt = indexedPeaks.MaxBy(p => p.Intensity).RetentionTime;
            var spline = LinearSpline.Interpolate(indexedPeaks.Select(p => p.RetentionTime).ToArray(),
                indexedPeaks.Select(p => 100 * p.Intensity / maxIntensity).ToArray());
            Func<Vector<double>, double> objectiveFunction = (parameters) =>
            {
                if (parameters.Count != 7)
                    throw new ArgumentException("Expected 7 parameters: FullWidthHalfMax, 3 Leading parameters, and 3 Tailing parameters.");
                EmpiricallyTransformedGaussian etgAlgo = new EmpiricallyTransformedGaussian
                {
                    Height = 200.0,
                    FullWidthHalfMax = parameters[0],
                    LambdaLeading = parameters[1],
                    KexpLeading = parameters[2],
                    AlphaLeading = parameters[3],
                    LambdaTailing = parameters[4],
                    KexpTailing = parameters[5],
                    AlphaTailing = parameters[6]
                };

                // Example calculation: return a simple function of the parameters
                double sse = 0.0;
                foreach (var peak in indexedPeaks)
                {
                    double simulatedIntensity = etgAlgo.GetIntensity(peak.RetentionTime, apexRt);
                    sse += Math.Pow((100*peak.Intensity / maxIntensity) - simulatedIntensity, 2);
                }

                return sse;
            };
            return objectiveFunction;
        }
    }
}
