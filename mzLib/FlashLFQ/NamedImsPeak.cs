using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FlashLFQ
{
    [Serializable]
    public class NamedImsPeak : IndexedMassSpectralPeak
    {
        public NamedImsPeak(double mz, double intensity, int zeroBasedMs1ScanIndex, double retentionTime, int isotopeCount) :
            base(mz, intensity, zeroBasedMs1ScanIndex, retentionTime)
        {
            //TODO: Define a get name function eg. GetName(charge: 2) i1z2
        }
    }
}
