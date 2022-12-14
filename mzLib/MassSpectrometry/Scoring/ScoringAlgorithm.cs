using System;
using System.Collections.Generic;
using System.Dynamic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;


namespace MassSpectrometry.Scoring
{
    public abstract class ScoringAlgorithm
    {
        public PpmTolerance Tolerance { get; }

        public ScoringAlgorithm(PpmTolerance tolerance)
        {
            Tolerance = tolerance;
        }


        public abstract double GetScore(double[] theoreticalMz, double[] theoreticalIntensity, double[] experimentalMz, double[] experimentalIntensity);

        /// <summary>
        /// Compares two scores in an algorithm specific fashion.
        /// </summary>
        /// <param name="instanceScore"></param>
        /// <param name="argumentScore"></param>
        /// <returns> +1 if the second score is better, -1 if the second score is worse, 0 if the two scores are equivalent</returns>
        public abstract int CompareScores(double instanceScore, double argumentScore);

        /// <summary>
        /// Finds matching peaks within the theoretical and experimental mzArrays, then stores the corresponding intensities
        /// to a two dimensional matrix. The first dimension (rows) represents theoretical (Index 0) vs experimental (Index 1).
        /// The second dimension (columns) represents paired intensity values
        /// </summary>
        /// <param name="theoretical"></param>
        /// <param name="experimental"></param>
        /// <returns></returns>
        public double[,] GetIntensityPairs(double[] theoreticalMz, double[] theoreticalIntensity, 
            double[] experimentalMz, double[] experimentalIntensity)
        {
            double[,] intensityPairs = new double[2, theoreticalMz.Length];
            double intensitySum = 0;
            for (int i = 0; i < theoreticalMz.Length; i++)
            {
                int expIndex = experimentalMz.GetNearestIndex(theoreticalMz[i]);
                if (Tolerance.Within(theoreticalMz[i], experimentalMz[expIndex]))
                {
                    intensityPairs[0, i] = theoreticalIntensity[i];
                    intensityPairs[1, i] = experimentalIntensity[expIndex];
                    intensitySum += intensityPairs[1, i];
                }
            }

            if (!(intensitySum > 0)) throw new ArgumentException("No paired peaks were found");
            return intensityPairs;
        }

        /// <summary>
        /// The default score comparison, where higher scores are better.
        /// </summary>
        /// <param name="instanceScore"></param>
        /// <param name="argumentScore"></param>
        /// <param name="betterScore"> The higher of the two scores</param>
        /// <returns> -1 if the instance score is higher than the argument score, 1 if the argument score is higher
        /// than the instance score, 0 if the scores are equal</returns>
        internal int DefaultCompare(double instanceScore, double argumentScore)
        {
            if (instanceScore > argumentScore)
            {
                return -1;
            }
            else if (instanceScore == argumentScore)
            {
                return 0;
            }
            else
            {
                return 1;
            }
        }
    }
}
