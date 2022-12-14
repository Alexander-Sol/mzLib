using System;
using System.Collections.Generic;
using System.Globalization;
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
    public PpmTolerance Tolerance { get; }
    public double ThresholdMz { get; }
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
    private double? _poorScore;

    public Scorer(ScoringScheme scoringScheme, NormalizationScheme normalization, PpmTolerance tolerance, double thresholdMz = 300)
    {
        ScoreType = scoringScheme;
        Normalization = normalization;
        Tolerance = tolerance;
        ThresholdMz = thresholdMz;
        ConstructScoringAlgorithm();
    }

    /// <summary>
    /// Here the experimental spectrum should be the longer of the two
    /// </summary>
    /// <param name="experimental"></param>
    /// <param name="theoretical"></param>
    /// <returns></returns>
    public double Score(ISpectralComparable experimental, ISpectralComparable theoretical)
    {
        double[] theoreticalMz = theoretical.GetMzArrayCopy();
        double[] experimentalMz = experimental.GetMzArrayCopy();
        double[] theoreticalIntensity = theoretical.GetIntensityArrayCopy();
        double[] experimentalIntensity = experimental.GetIntensityArrayCopy();

        if (theoreticalIntensity.Length == 0 | experimentalIntensity.Length == 0)
        {
            throw new MzLibException(string.Format(CultureInfo.InvariantCulture, "Empty YArray in spectrum."));
        }
        if (theoreticalIntensity.Length != theoreticalMz.Length | experimentalIntensity.Length != experimentalMz.Length)
        {
            throw new MzLibException(string.Format(CultureInfo.InvariantCulture, "Discordance between length of mz and intensity arrays"));
        }

        TruncateBelow(theoreticalMz, theoreticalIntensity, 300);
        TruncateBelow(experimentalMz, experimentalIntensity, 300);

        Normalize(theoreticalIntensity);
        Normalize(experimentalIntensity);

        return Algorithm.GetScore(theoreticalMz, theoreticalIntensity, experimentalMz, experimentalIntensity);
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
        double cumulativeIntensity = 0;
        for (int i = 0; i < mzArray.Length; i++)
        {
            if (mzArray[i] < thresholdMz)
            {
                mzArray[i] = 0;
                intensityArray[i] = 0;
            }
            cumulativeIntensity += intensityArray[i];
        }

        if (!(cumulativeIntensity > 0)) throw new MzLibException(string.Format(CultureInfo.InvariantCulture, "Spectrum has no intensity after truncation."));
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
        if (Algorithm.CompareScores(instanceScore, argumentScore) > 0)
        {
            betterScore = argumentScore;
            return true;
        }

        betterScore = instanceScore;
        return false;
    }



    /// <summary>
    /// Assigns the specific scoring algorithm that will be used by the Scorer
    /// </summary>
    /// <exception cref="NotImplementedException"></exception>
    private void ConstructScoringAlgorithm()
    {
        switch (ScoreType)
        {
            case ScoringScheme.KullbackLeibler:
                throw new NotImplementedException();
            case ScoringScheme.SpectralContrastAngle:
                Algorithm = new SpectralContrastAlgorithm(Tolerance);
                break;
            default:
                throw new NotImplementedException();
        }
    }

}