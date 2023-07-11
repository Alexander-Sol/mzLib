using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    public class VerboseIsotopicEnvelope : IsotopicEnvelope
    {
        public VerboseIsotopicEnvelope(
            IndexedMassSpectralPeak mostAbundantPeak, 
            List<IndexedMassSpectralPeak> allPeaks,
            int chargeState, 
            double intensity = -1) :
            base(mostAbundantPeak, chargeState, intensity > -1 ? intensity : allPeaks.Sum(p => p.Intensity))
        {
            AllPeaks = allPeaks.OrderBy(peak => peak.Mz).ToList();
        }

        public List<IndexedMassSpectralPeak> AllPeaks { get; }

        public override string ToString()
        {
            return "+" + ChargeState + 
                   "|" + Intensity.ToString("F0") + 
                   "|" + IndexedPeak.RetentionTime.ToString("F3") + 
                   "|" + IndexedPeak.ZeroBasedMs1ScanIndex;
        }
    }
}
