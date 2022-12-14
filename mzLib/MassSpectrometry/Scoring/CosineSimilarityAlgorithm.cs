using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry.Scoring
{
    internal class CosineSimilarityAlgorithm : ScoringAlgorithm
    {
        /// <summary>
        /// The cosine similarity returns values between 1 and -1 with 1 being closest, -1 being opposite, and 0 being orthogonal
        /// </summary>
        /// <param name="tolerance"></param>
        public CosineSimilarityAlgorithm(PpmTolerance tolerance) : base(tolerance)
        {

        }

        public override int CompareScores(double instanceScore, double argumentScore)
        {
            return DefaultCompare(instanceScore, argumentScore);
        }

        /// <summary>
        /// Returns a Cosine Similarity Score 1 and -1 with 1 being closest, -1 being opposite, and 0 being orthogonal
        /// </summary>
        /// <param name="tolerance"></param>
        public override double GetScore(double[] experimentalMz, double[] experimentalIntensity,
            double[] theoreticalMz, double[] theoreticalIntensity)
        {
            double[,] intensityPairs = new double[2, theoreticalMz.Length];
            try
            {
                intensityPairs = GetIntensityPairs(theoreticalMz, theoreticalIntensity, experimentalMz,
                    experimentalIntensity);
            }
            catch (ArgumentException)
            {
                return 0;
            }

            double numerator = 0;
            double denominatorValue1 = 0;
            double denominatorValue2 = 0;
            for (int i = 0; i < theoreticalMz.Length; i++)
            {
                numerator += intensityPairs[0, i] * intensityPairs[1, i];
                denominatorValue1 += Math.Pow(intensityPairs[0, i], 2);
                denominatorValue2 += Math.Pow(intensityPairs[1, i], 2);
            }
            double denominatorProduct = denominatorValue1 * denominatorValue2;

            //because we keep all secondary spectrum peaks, denominatorValue1 can equal zero
            if (denominatorProduct == 0)
            {
                return 0;
            }

            return numerator / Math.Sqrt(denominatorProduct);
        }
    }
}
