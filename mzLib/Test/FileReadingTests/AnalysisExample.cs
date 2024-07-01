using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using System.Xml.Serialization;
using Chemistry;
using MassSpectrometry;
using MzIdentML;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using Readers;
using Stopwatch = System.Diagnostics.Stopwatch;
using Plotly.NET.CSharp;

namespace Readers
{
    [TestFixture]
    internal class AnalysisExample
    {

        [Test]
        public static void ControlXIC()
        {
            string mzmlPath = @"D:\Human_Ecoli_TwoProteome_60minGradient\CalibrateSearch_4_19_24\Human_Calibrated_Files\04-12-24_Human_C18_3mm_50msec_stnd-60min_1-calib.mzML";
            var reader = MsDataFileReader.GetDataFile(mzmlPath);
            reader.LoadAllStaticData();
            var ms2Scans = reader.GetAllScansList().Where(scan => scan.MsnOrder > 1).ToList();
            Tolerance tolerance = new PpmTolerance(20);
            double[] rtArray = new double[ms2Scans.Count];
            double[] intensityArray = new double[ms2Scans.Count];

            for (int i = 0; i < ms2Scans.Count; i++)
            {
                rtArray[i] = ms2Scans[i].RetentionTime;
                intensityArray[i] = 0;
                int idx240 = ms2Scans[i].MassSpectrum.GetClosestPeakIndex(240.17);
                if (!tolerance.Within(ms2Scans[i].MassSpectrum.XArray[idx240], 240.17)) continue;

                int idx509 = ms2Scans[i].MassSpectrum.GetClosestPeakIndex(509.31);
                if (!tolerance.Within(ms2Scans[i].MassSpectrum.XArray[idx509], 509.31)) continue;

                intensityArray[i] = ms2Scans[i].MassSpectrum.YArray[idx240] + ms2Scans[i].MassSpectrum.YArray[idx509];
            }

            Chart.Line<double, double, string>(rtArray, intensityArray)
                .WithTraceInfo("Diagnostic Ions")
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Retention Time (min)"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity of Diagnostic Ions"))
                .WithSize(Width: 1000, Height: 500)
                .Show();
        }
    }
}
