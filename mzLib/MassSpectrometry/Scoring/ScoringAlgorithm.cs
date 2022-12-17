using System;
using System.Collections.Generic;
using System.Dynamic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using MzLibUtil;


namespace MassSpectrometry.Scoring
{
    public abstract class ScoringAlgorithm
    {
        public PpmTolerance Tolerance { get; }

        // If AllPeaks is set to true, experimental peaks not found in the
        // theoretical spectrum are used for scoring. Otherwise, they are discarded
        public bool AllPeaks { get; set; }

        public ScoringAlgorithm(PpmTolerance tolerance, bool allPeaks = true)
        {
            Tolerance = tolerance;
            AllPeaks = allPeaks;
        }


        public abstract double GetScore(double[] experimentalMz, double[] experimentalIntensity, double[] theoreticalMz, double[] theoreticalIntensity);

        /// <summary>
        /// Compares two scores in an algorithm specific fashion.
        /// </summary>
        /// <param name="instanceScore"></param>
        /// <param name="argumentScore"></param>
        /// <returns> +1 if the second score is better, -1 if the second score is worse, 0 if the two scores are equivalent</returns>
        public abstract int CompareScores(double instanceScore, double argumentScore);

        /// <summary>
        /// Finds matching peaks within the theoretical and experimental mzArrays, then stores the corresponding intensities
        /// to a list of (double, double) tuples. The first entry in the tuple corresponds to the theoretical
        /// intensity and the second entry corresponds to the experimental intensity
        /// </summary>
        /// <param name="theoretical"></param>
        /// <param name="experimental"></param>
        /// <returns></returns>
        public List<(double, double)> GetIntensityPairs(double[] experimentalMz, double[] experimentalIntensity,
            double[] theoreticalMz, double[] theoreticalIntensity)
        {
            List<(double, double)> intensityPairs = new List<(double theoretical, double experimental)>();
            SpectrumTree experimentalTree = new SpectrumTree(experimentalMz, experimentalIntensity);
            int matchesFound = 0;
            for (int i = 0; i < theoreticalMz.Length; i++)
            {
                if (experimentalTree.PopClosestPeak(theoreticalMz[i], Tolerance,
                        out var experimentalPeakMz, out var experimentalPeakIntensity))
                {
                    intensityPairs.Add((theoreticalIntensity[i], experimentalPeakIntensity));
                    matchesFound++;
                }
                else
                {
                    intensityPairs.Add((theoreticalIntensity[i], 0));
                }
            }
            if (matchesFound == 0) throw new ArgumentException("No matches exist");
            if (AllPeaks)
            {
                foreach (Node node in experimentalTree)
                {
                    intensityPairs.Add((0, node.Value));
                }
            }

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
