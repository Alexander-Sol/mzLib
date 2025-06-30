using MathNet.Numerics.Interpolation;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Optimization;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MassSpectrometry
{
    public class EtgObjectiveFunction : IObjectiveFunction
    {
        public int ParameterCount => 7; // Height, FullWidthHalfMax, and Leading/Tailing parameters

        public Vector<double> Point => throw new NotImplementedException();

        public double Value => throw new NotImplementedException();

        public bool IsGradientSupported => false;

        public Vector<double> Gradient => throw new NotImplementedException();

        public bool IsHessianSupported => false;

        public Matrix<double> Hessian => throw new NotImplementedException();

        public EtgObjectiveFunction(List<IIndexedPeak> indexedPeaks)
        {
            var maxIntensity = indexedPeaks.Max(p => p.Intensity);
            LinearSpline = LinearSpline.Interpolate(indexedPeaks.Select(p => p.RetentionTime).ToArray(), 
                indexedPeaks.Select(p => 100* p.Intensity / maxIntensity).ToArray());
        }

        LinearSpline LinearSpline { get; }

        public void EvaluateAt(Vector<double> points)
        {
            double sse = 0.0;
            foreach(var point in points)
            {
                var expIntensity = LinearSpline.Interpolate(point);
            }
            
        }


        public static Func<Vector<double>, double> CreateObjectiveFunction(List<IIndexedPeak> indexedPeaks)
        {
            var maxIntensity = indexedPeaks.Max(p => p.Intensity);
            var apexRt = indexedPeaks.MaxBy(p => p.Intensity).RetentionTime;
            var spline = LinearSpline.Interpolate(indexedPeaks.Select(p => p.RetentionTime).ToArray(), 
                indexedPeaks.Select(p => 100 * p.Intensity / maxIntensity).ToArray());
            Func<Vector<double>, double> objectiveFunction = (parameters) =>
            {
                if (parameters.Count != 2)
                    throw new ArgumentException("Expected 7 parameters: FullWidthHalfMax, 3 Leading parameters, and 3 Tailing parameters.");
                EmpiricallyTransformedGaussian etgAlgo = new EmpiricallyTransformedGaussian
                {
                    Height = 200.0,
                    FullWidthHalfMax = parameters[0],
                    LambdaLeading = parameters[1],
                    KexpLeading = parameters[2],
                    AlphaLeading = parameters[3],
                    LambdaTailing = parameters[4],
                    KexpTailing = parameters[5],
                    AlphaTailing = parameters[6]
                };

                // Example calculation: return a simple function of the parameters
                double sse = 0.0;
                foreach (var peak in indexedPeaks)
                {
                    double expectedIntensity = etgAlgo.GetIntensity(peak.RetentionTime, apexRt);
                    sse += Math.Pow( ( 100 * peak.Intensity / maxIntensity ) - expectedIntensity, 2);
                }
                
                return sse;
            };
            return objectiveFunction;
        }

        public double ValueAt(double[] parameters)
        {
            if (parameters.Length != ParameterCount)
                throw new ArgumentException($"Expected {ParameterCount} parameters, but got {parameters.Length}.");
            double height = parameters[0];
            double fullWidthHalfMax = parameters[1];
            double leading = parameters[2]; // This could be a combination of leading and tailing parameters
            // Example calculation: return a simple function of the parameters
            return height * Math.Exp(-Math.Pow(fullWidthHalfMax, 2) / (2 * leading * leading));
        }
        public IObjectiveFunction DerivativeAt(double[] parameters)
        {
            // Implement derivative logic if needed
            throw new NotImplementedException("Derivative calculation is not implemented.");
        }


        public IObjectiveFunction Fork()
        {
            throw new NotImplementedException();
        }

        public IObjectiveFunction CreateNew()
        {
            throw new NotImplementedException();
        }
    }
}
