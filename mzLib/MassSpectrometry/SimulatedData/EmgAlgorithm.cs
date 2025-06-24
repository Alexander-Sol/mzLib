using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Meta.Numerics;
using Complex = Meta.Numerics.Complex;
using Meta.Numerics.Functions;


namespace MassSpectrometry;
public class EmgAlgorithm
{
    /// <summary>
    /// This class implements the EMG algorithm for peak detection in mass spectrometry data.
    /// </summary>
    public EmgAlgorithm(double stdDev, double skewFactor, double intensity)
    {
        StandardDeviation = stdDev;
        SkewFactor = skewFactor;
        Intensity = intensity;
    }

    public EmgAlgorithm()
    {
        // Default constructor with default values
    }

    public double StandardDeviation { get; set; } = 0.1;
    /// <summary>
    /// Skew factor less than zero is fronting, greater than zero is tailing.
    /// </summary>
    public double SkewFactor { get; set; } = 0.5;
    public double Intensity { get; set; } = 1000.0;

    public double GetIntensity(double rt, double apexRt)
    {
        double term1 = Intensity * apexRt * Math.Sqrt(2 * Math.PI) / (SkewFactor * 2);
        double term2 = Math.Exp( ( (StandardDeviation - rt) / SkewFactor) * (Math.Pow(apexRt, 2) / (2 * Math.Pow(SkewFactor, 2) )));

        Complex complexTermForErf = SkewFactor < 0 ?
            new Complex(re: (StandardDeviation - rt) / Math.Sqrt(2 * apexRt), im: (apexRt / Math.Sqrt(2 * Math.Abs(SkewFactor)))) :
            new Complex(re: ((StandardDeviation - rt) / Math.Sqrt(2 * apexRt)) + (apexRt / Math.Sqrt(2 * SkewFactor)), im: 0) ;
        double term3 = (double)Math.Sign(SkewFactor) - AdvancedComplexMath.Erf(complexTermForErf).Re;

        return term1 * term2 * term3;
    }

    public List<double> GetIntensityRange(double apexRt, List<double> rts)
    {
        List<double> intensities = new List<double>();
        foreach (var rt in rts)
        {
            double intensity = GetIntensity(rt, apexRt);
            intensities.Add(intensity);     
        }
        return intensities;
    }

    //static double Erf(double x)
    //{
    //    // constants
    //    double a1 = 0.254829592;
    //    double a2 = -0.284496736;
    //    double a3 = 1.421413741;
    //    double a4 = -1.453152027;
    //    double a5 = 1.061405429;
    //    double p = 0.3275911;

    //    // Save the sign of x
    //    int sign = 1;
    //    if (x < 0)
    //        sign = -1;
    //    x = Math.Abs(x);

    //    // A&S formula 7.1.26
    //    double t = 1.0 / (1.0 + p * x);
    //    double y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.Exp(-x * x);

    //    return sign * y;
    //}
}

