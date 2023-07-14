using Chemistry;
using MathNet.Numerics;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    public class VerboseIsotopicEnvelope : IsotopicEnvelope
    {
        /// <summary>
        /// Dictionary where the key represents the isotopic distance from the monoisotopic species. 
        /// Key 0 represent the monoisotopic peak, 1 represent the first isotope peak, etc. 
        /// It is assumed that all isotopic peaks are C13-C12 mz apart
        /// </summary>
        public Dictionary<int, IndexedMassSpectralPeak> PeakDictionary { get; }

        public VerboseIsotopicEnvelope(
            IndexedMassSpectralPeak mostAbundantPeak, 
            List<IndexedMassSpectralPeak> allPeaks,
            int chargeState, 
            double monoisotopicMass,
            int isotopePpmTolerance = 5,
            double intensity = -1) :
            base(mostAbundantPeak, chargeState, intensity > -1 ? intensity : allPeaks.Sum(p => p.Intensity))
        {
            PeakDictionary = WritePeakDictionary(
                allPeaks.OrderBy(peak => peak.Mz).ToList(), 
                monoisotopicMass.ToMz(chargeState), 
                isotopePpmTolerance);
            RetentionTime = mostAbundantPeak.RetentionTime.Round(4);
        }

        public double RetentionTime { get; }

        public override string ToString()
        {
            return "+" + ChargeState + 
                   "|" + Intensity.ToString("F0") + 
                   "|" + IndexedPeak.RetentionTime.ToString("F3") + 
                   "|" + IndexedPeak.ZeroBasedMs1ScanIndex;
        }

        public static Dictionary<int, IndexedMassSpectralPeak> WritePeakDictionary(
            List<IndexedMassSpectralPeak> peaks, double monoisotopicMz, int isotopePpmTolerance)
        {
            Dictionary<int, IndexedMassSpectralPeak> peakDict = new();
            PpmTolerance tolerance = new PpmTolerance(isotopePpmTolerance);
            int peakListIndex = 0;
            int isotopeIndex = 0;
            while(peakDict.Count < peaks.Count)
            {
                if (tolerance.Within(peaks[peakListIndex].Mz, monoisotopicMz + Chemistry.Constants.C13MinusC12 * isotopeIndex))
                    peakDict.Add(isotopeIndex++, peaks[peakListIndex++]);
                else if (tolerance.GetMinimumValue(monoisotopicMz + Chemistry.Constants.C13MinusC12 * isotopeIndex) > peaks[peakListIndex].Mz)
                    return WritePeakDictionary(peaks, monoisotopicMz, isotopePpmTolerance + 5);
                else
                    isotopeIndex++;
            }
            return peakDict;
        }
        public static string GetIsotopePeakName((int isotopeNumber, int chargeState) key)
        {
            return "i" + key.isotopeNumber + "z" + key.chargeState;
        }
    }
}
