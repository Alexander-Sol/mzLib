using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry.Scoring;
using MzLibUtil;

namespace MassSpectrometry.Scoring;
public class Scorer
{
    public enum ScoringScheme
    {
        KullbackLeibler,
        SpectralContrastAngle
    }

    public enum NormalizationScheme
    {
        squareRootSpectrumSum,
        spectrumSum,
        mostAbundantPeak,
        unnormalized
    }

    public ScoringScheme ScoreType { get; }
    public NormalizationScheme Normalization { get; }
    public ScoringAlgorithm Algorithm { get; private set; }
    private double? _poorScore;
    public PpmTolerance Tolerance { get; }

    public double PoorScore
    {
        get
        {
            if (_poorScore.HasValue) return (double)_poorScore;
            switch (ScoreType)
            {
                case ScoringScheme.KullbackLeibler:
                    _poorScore = Double.MaxValue;
                    return (double)_poorScore;
                case ScoringScheme.SpectralContrastAngle:
                    _poorScore = 0;
                    return (double)_poorScore;
                default:
                    _poorScore = Double.MinValue;
                    return (double)_poorScore;
            }
        }
    }

    public Scorer(ScoringScheme scoringScheme, NormalizationScheme normalization, PpmTolerance tolerance)
    {
        ScoreType = scoringScheme;
        Normalization = normalization;
        ConstructScoringAlgorithm(tolerance);
    }

    /// <summary>
    /// For any peak with m/z less than the threshold m/z, m/z and corresponding intensity are
    /// set to 0 via modification in place.
    /// </summary>
    /// <param name="mzArray"></param>
    /// <param name="intensityArray"></param>
    /// <param name="thresholdMz"></param>
    public void TruncateBelow(double[] mzArray, double[] intensityArray, double thresholdMz)
    {
        for (int i = 0; i < mzArray.Length; i++)
        {
            if (mzArray[i] < thresholdMz)
            {
                mzArray[i] = 0;
                intensityArray[i] = 0;
            }
        }
    }

    /// <summary>
    /// Normalizes a given intensity array by modifying in place
    /// </summary>
    /// <param name="intensityArray"></param>
    public void Normalize(double[] intensityArray)
    {
        Func<double, double> normalizationCalculation = null;
        switch (Normalization)
        {
            case NormalizationScheme.mostAbundantPeak:
                double maxAbundance = intensityArray.Max();
                normalizationCalculation = x => x / maxAbundance;
                break;
            case NormalizationScheme.spectrumSum:
                double spectrumSum = intensityArray.Sum();
                normalizationCalculation = x => x / spectrumSum;
                break;
            case NormalizationScheme.squareRootSpectrumSum:
                double squareRootSpectrumSum = intensityArray.Select(y => Math.Sqrt(y)).Sum();
                normalizationCalculation = x => Math.Sqrt(x) / squareRootSpectrumSum;
                break;
            case NormalizationScheme.unnormalized:
                return; // No normalization required 
        }

        for (int i = 0; i < intensityArray.Length; i++)
        {
            intensityArray[i] = normalizationCalculation(intensityArray[i]);
        }
    }

    public List<(double, double)> GetIntensityPairs(double[] theoretical, double[] experimental)
    {
        List<(double theoretical, double experimental)> intensityPairs = new List<(double, double)> ();
        for (int i = 0; i < theoretical.Length; i++)
        {
            int expIndex = experimental.GetNearestIndex(theoretical[i]);
            if (Tolerance.Within(theoretical[i], experimental[expIndex]) )
            {
                intensityPairs.Add((theoretical[i], experimental[expIndex]));
            }
        }

        return intensityPairs;
    }

    //public double Score(IScoreArgs args)
    //{
    //    return GetScore(args);
    //}

    /// <summary>
    /// Here the experimental spectrum should be the longer of the two
    /// </summary>
    /// <param name="experimental"></param>
    /// <param name="theoretical"></param>
    /// <returns></returns>
    public double Score(ISpectralComparable experimental, ISpectralComparable theoretical)
    {
        throw new NotImplementedException();
    }

//public double Score(MinimalSpectrum experimentalSpectrum, MinimalSpectrum theoreticalSpectrum)
    //{
    //    IScoreArgs args = new MinimalSpectraArgs(experimentalSpectrum, theoreticalSpectrum);
    //    return ScoringScheme.GetScore(args);
    //}

    /// <summary>
    /// Compares two scores in a method specific fashion. Returns true if the instanceScore (first)
    /// is better than the argumentScore (second). Outputs the better of the two. This method is necessary
    /// because there are some metrics where lower scores are better.
    /// </summary>
    /// <param name="instanceScore"></param>
    /// <param name="argumentScore"></param>
    /// <param name="betterScore"></param>
    /// <returns></returns>
    /// <exception cref="NotImplementedException"></exception>
    public bool TestForScoreImprovement(double instanceScore, double argumentScore, out double betterScore)
    {
        switch (ScoreType)
        {
            case ScoringScheme.KullbackLeibler:
                if (instanceScore < argumentScore)
                {
                    betterScore = instanceScore;
                    return true;
                }
                else
                {
                    betterScore = argumentScore;
                    return false;
                }
            case ScoringScheme.SpectralContrastAngle:
                return DefaultCompare(instanceScore, argumentScore, out betterScore);
            default:
                return DefaultCompare(instanceScore, argumentScore, out betterScore);
        }
    }

    /// <summary>
    /// The default score comparison, where higher scores are better. Compares two scores, returns true
    /// if the instance score is higher than the argument score, returns false if instance score is lower.
    /// </summary>
    /// <param name="instanceScore"></param>
    /// <param name="argumentScore"></param>
    /// <param name="betterScore"> The higher of the two scores</param>
    /// <returns></returns>
    private bool DefaultCompare(double instanceScore, double argumentScore, out double betterScore)
    {
        if (instanceScore > argumentScore)
        {
            betterScore = instanceScore;
            return true;
        }
        else
        {
            betterScore = argumentScore;
            return false;
        }
    }

    private void ConstructScoringAlgorithm(PpmTolerance tolerance)
    {
        switch (ScoreType)
        {
            case ScoringScheme.KullbackLeibler:
                throw new NotImplementedException();
            case ScoringScheme.SpectralContrastAngle:
                Algorithm = new SpectralContrastAlgorithm(tolerance);
                break;
            default:
                throw new NotImplementedException();
        }
    }

}