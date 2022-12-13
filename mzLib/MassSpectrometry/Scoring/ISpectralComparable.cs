using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry.Scoring
{
    public interface ISpectralComparable
    {
        public double[] GetMzArrayCopy();
        public double[] GetIntensityArrayCopy();
    }
}
