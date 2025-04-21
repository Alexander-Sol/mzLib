using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzIdentML;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using Readers;
using Stopwatch = System.Diagnostics.Stopwatch;
using Plotly.NET.CSharp;
using MathNet.Numerics.Distributions;
using static Plotly.NET.StyleParam.DrawingStyle;

namespace Test.FileReadingTests
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

        [Test]
        public static void AnalyzeFlashData()
        {
            string dataPath = @"D:\JurkatTopdown\MM108_DecoyPTMs_MOxVariable\FlashLFQ_Tweaked\QuantifiedPeaks.tsv";
            var peaks = FileReader.ReadFile<QuantifiedPeakFile>(dataPath);
            int mixedDecoyRealModPeaks = 0;
            int successfulDisambiguations = 0;
            int totalDecoyMods = 0;
            foreach (var peak in peaks)
            {
                if (!peak.FullSequence.Contains(":Decoy")) continue;
                if (peak.BestMatchingSequence.Contains(":Decoy"))
                {
                    totalDecoyMods++;
                }
                var x = peak.FullSequence.Split('|');
                if (x.Length < 2) continue;
                if(x.Any(s => s.Contains(":Decoy")) && !x.All(s => s.Contains(":Decoy")))
                {
                    mixedDecoyRealModPeaks++;
                    if (!peak.BestMatchingSequence.Contains(":Decoy")) successfulDisambiguations++;
                }
            }

            int placeholder = 0;
        }

        [Test]
        public static void GetIsoEnv()
        {
            string histoneSeq = "MPEPAKSAPAPKKGSKKAVTKAQKKDGKKRKRSRKESYSVYVYKVLKQVHPDTGISSKAMGIMNSFVNDIFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSAK";

            ChemicalFormula cf = new Peptide(histoneSeq).GetChemicalFormula();
            IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
            double[] mz = dist.Masses.Select(v => v.ToMz(20)).ToArray();
            double[] intensities = dist.Intensities.Select(v => v * 100).ToArray();
            double rt = 1;


            double[] plotMz = new double[mz.Length * 3];
            double[] plotIntensity = new double[mz.Length * 3];

            for (int i = 0; i < mz.Length; i++)
            {
                int first = (3 * i);
                int second = (3 * i + 1);
                int third = (3 * i + 2);

                plotMz[first] = mz[i] - 0.00001;
                plotIntensity[first] = 0;

                plotMz[second] = mz[i];
                plotIntensity[second] = intensities[i];

                plotMz[third] = mz[i] + 0.00001;
                plotIntensity[third] = 0;
            }
            

            // add the scan
            MsDataScan scan = new MsDataScan(massSpectrum: new MzSpectrum(plotMz, plotIntensity, false), oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: rt, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (1));

            Chart.Line<double, double, string>(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray)
                .WithTraceInfo("Theoretical Envelope")
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
                .WithSize(Width: 1000, Height: 500)
                .Show();

        }


        [Test]
        public static void GetCombinedIsoEnv()
        {
            string histoneSeq = "PEPAKSAPAPKKGSKKAVTKAQKKDGKKRKRSRKESYSIYVYKVLKQVHPDTGISSKAMGIMNSFVNDIFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSSK";

            int baseScalingFactor = 75000;
            int modScalingFactor = 25000;
            //int z = 30;

            ChemicalFormula cf = new Peptide(histoneSeq).GetChemicalFormula();
            IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);

            ChemicalFormula methyl = ChemicalFormula.ParseFormula("CH2");
            ChemicalFormula acetyl = ChemicalFormula.ParseFormula("C2H2O");
            ChemicalFormula oxidation = ChemicalFormula.ParseFormula("O");

            ChemicalFormula methylH = ChemicalFormula.Combine(new List<ChemicalFormula> { methyl, cf });
            ChemicalFormula acetylH = ChemicalFormula.Combine(new List<ChemicalFormula> { acetyl, cf });
            ChemicalFormula acetylMethylH = ChemicalFormula.Combine(new List<ChemicalFormula> { acetyl, methyl, cf });
            ChemicalFormula phospho = ChemicalFormula.Combine(new List<ChemicalFormula> { ChemicalFormula.ParseFormula("PO4H3"), cf });
            //ChemicalFormula methionineLoss = ChemicalFormula.Combine(new List<ChemicalFormula> { ChemicalFormula.ParseFormula("PO4H3"), cf })

            List<ChemicalFormula> chemicalFormulas = new List<ChemicalFormula> { cf, methylH, acetylH, phospho, acetylMethylH };
            chemicalFormulas = chemicalFormulas.SelectMany(f => 
                new List<ChemicalFormula>() { f, ChemicalFormula.Combine( new List<ChemicalFormula> { f, oxidation }) })
                .ToList();
            List<IsotopicDistribution> isoDistributions = chemicalFormulas.Select(cf => IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8)).ToList();

            List<double[]> mzArrayList = new();
            List<int[]> intensityArrayList = new();

            Random randomScaling = new();

            for(int z = 14; z < 23; z++)
            {
                int chargeStateScalingFactor = 6 - Math.Abs(z - 18);
                double randomScalingFactor = randomScaling.NextDouble() * 0.2 + 0.9;

                // unmodified
                mzArrayList.Add(dist.Masses.Select(v => v.ToMz(z)).ToArray());
                intensityArrayList.Add(dist.Intensities.Select(v => (int)(v * baseScalingFactor * chargeStateScalingFactor * randomScalingFactor)).ToArray());

                foreach(var distribution in isoDistributions)
                {
                    randomScalingFactor = randomScaling.NextDouble() * 0.2 + 0.9;
                    mzArrayList.Add(distribution.Masses.Select(v => v.ToMz(z)).ToArray());
                    intensityArrayList.Add(distribution.Intensities.Select(v => (int)(v * modScalingFactor * chargeStateScalingFactor * randomScalingFactor)).ToArray());
                }
            }

            var noiseParams = new LowFrequencyNoiseParameters(peakNumberLimitLow:5999, peakNumberLimitHigh: 6000, 
                peakLocationLimitLow: 625, peakLocationLimitHigh: 1000, peakIntensityLimitHigh: 1500);

            var lowFreqNoiseSpectrum = SimulatedData.BuildLowFrequencyNoiseSpectrum(noiseParams);
            mzArrayList.Add(lowFreqNoiseSpectrum.XArray);
            intensityArrayList.Add(lowFreqNoiseSpectrum.YArray.Select(d => (int)d).ToArray());

            var displaySpectrum = TofSpectraMerger.MergeArraysToMs2Spectrum(mzArrayList, intensityArrayList, ppmTolerance: 5);

            PlotSpectrum(displaySpectrum.XArray, displaySpectrum.YArray);

            SimulatedData.AddHighFrequencyNoiseToYArray(displaySpectrum.YArray, new Normal(250, 1500));

            MsDataScan[] scans = new MsDataScan[5];

            // add the scan
            scans[0] = new MsDataScan(massSpectrum: displaySpectrum, oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: 1, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: displaySpectrum.SumOfAllY, injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (1),
                isolationMZ: 692, isolationWidth: 4);

            scans[1] = new MsDataScan(massSpectrum: displaySpectrum, oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: 1, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: displaySpectrum.SumOfAllY, injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (1),
                isolationMZ: 729, isolationWidth: 4);

            scans[2] = new MsDataScan(massSpectrum: displaySpectrum, oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: 1, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: displaySpectrum.SumOfAllY, injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (1),
                isolationMZ: 770, isolationWidth: 4);

            scans[3] = new MsDataScan(massSpectrum: displaySpectrum, oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: 1, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: displaySpectrum.SumOfAllY, injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (1),
                isolationMZ: 815, isolationWidth: 4);

            scans[4] = new MsDataScan(massSpectrum: displaySpectrum, oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: 1, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: displaySpectrum.SumOfAllY, injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (1),
                isolationMZ: 866, isolationWidth: 4);

            var myMsDataFile = new FakeMsDataFile(scans);
            //MsDataScan scan = myMsDataFile.GetAllScansList()[0];

            DeconvolutionParameters deconParams = new ClassicDeconvolutionParameters(1, 60, 4, 3);

            var precursors = new List<IsotopicEnvelope>();
            foreach(var scan in myMsDataFile.GetAllScansList())
            {
                precursors.AddRange(scan.GetIsolatedMassesAndCharges(scan.MassSpectrum, deconParams));
            }

            var deconvolutedMonoisotopicMasses = precursors.Select(p => p.MonoisotopicMass).OrderBy(d => d).Distinct().ToList();
            var actualMonoisotopicMasses = chemicalFormulas.Select(cf => cf.MonoisotopicMass).OrderBy(d => d).ToList();

            Console.WriteLine("Ground truth monoisotopic masses:\n" + String.Join("\n", actualMonoisotopicMasses));
            Console.WriteLine("Deconvoluted monoisotopic masses:\n" + String.Join("\n", deconvolutedMonoisotopicMasses));

            using(var sw = new StreamWriter(@"C:\Users\Alex\Documents\TopDownVision\ActualVsDecon.tsv"))
            {
                sw.Write("Actual\t");
                sw.WriteLine(string.Join('\t', actualMonoisotopicMasses));
                sw.WriteLine(string.Join('\t', deconvolutedMonoisotopicMasses.OrderBy(m => m).ToList()));
            }

            PlotSpectrum(displaySpectrum.XArray, displaySpectrum.YArray);
        }


        public static void PlotSpectrum(double[] xArray, double[] yArray)
        {
            double[] plotMz = new double[xArray.Length * 3];
            double[] plotIntensity = new double[yArray.Length * 3];

            for (int i = 0; i < xArray.Length; i++)
            {
                int first = (3 * i);
                int second = (3 * i + 1);
                int third = (3 * i + 2);

                plotMz[first] = xArray[i] - 0.00001;
                plotIntensity[first] = 0;

                plotMz[second] = xArray[i];
                plotIntensity[second] = Math.Max(0, yArray[i]);

                plotMz[third] = xArray[i] + 0.00001;
                plotIntensity[third] = 0;
            }

            Chart.Line<double, double, string>(plotMz, plotIntensity)
                .WithTraceInfo("Theoretical Envelope")
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
                .WithSize(Width: 1000, Height: 500)
                .Show();
        }

        [Test]
        public static void Ms1Example()
        {
            string mzmlPath = @"D:\Human_Ecoli_TwoProteome_60minGradient\CalibrateSearch_4_19_24\Human_Calibrated_Files\04-12-24_Human_C18_3mm_50msec_stnd-60min_1-calib.mzML";
            var reader = MsDataFileReader.GetDataFile(mzmlPath);
            reader.LoadAllStaticData();
            var ms2Scans = reader.GetAllScansList().Where(scan => scan.MsnOrder == 1).ToList();
            Tolerance tolerance = new PpmTolerance(20);
            double[] rtArray = new double[ms2Scans.Count];
            double[] intensityArray = new double[ms2Scans.Count];


            var scan = ms2Scans[30 + ms2Scans.Count / 2];

            //for (int i = 0; i < ms2Scans.Count; i++)
            //{
            //    rtArray[i] = ms2Scans[i].RetentionTime;
            //    intensityArray[i] = 0;
            //    int idx240 = ms2Scans[i].MassSpectrum.GetClosestPeakIndex(240.17);
            //    if (!tolerance.Within(ms2Scans[i].MassSpectrum.XArray[idx240], 240.17)) continue;

            //    int idx509 = ms2Scans[i].MassSpectrum.GetClosestPeakIndex(509.31);
            //    if (!tolerance.Within(ms2Scans[i].MassSpectrum.XArray[idx509], 509.31)) continue;

            //    intensityArray[i] = ms2Scans[i].MassSpectrum.YArray[idx240] + ms2Scans[i].MassSpectrum.YArray[idx509];
            //}

            Chart.Point<double, double, string>(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray)
                .WithTraceInfo("Diagnostic Ions")
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
                .WithSize(Width: 1000, Height: 500)
                .Show();
        }

        [Test]
        public static void Ms2LogExample()
        {
            string mzmlPath = @"D:\Human_Ecoli_TwoProteome_60minGradient\CalibrateSearch_4_19_24\Human_Calibrated_Files\04-12-24_Human_C18_3mm_50msec_stnd-60min_1-calib.mzML";
            var reader = MsDataFileReader.GetDataFile(mzmlPath);
            reader.LoadAllStaticData();
            var ms2Scans = reader.GetAllScansList().Where(scan => scan.MsnOrder == 1).ToList();
            Tolerance tolerance = new PpmTolerance(20);
            double[] rtArray = new double[ms2Scans.Count];
            double[] intensityArray = new double[ms2Scans.Count];


            var scan = ms2Scans[30 + ms2Scans.Count / 2];

            //for (int i = 0; i < ms2Scans.Count; i++)
            //{
            //    rtArray[i] = ms2Scans[i].RetentionTime;
            //    intensityArray[i] = 0;
            //    int idx240 = ms2Scans[i].MassSpectrum.GetClosestPeakIndex(240.17);
            //    if (!tolerance.Within(ms2Scans[i].MassSpectrum.XArray[idx240], 240.17)) continue;

            //    int idx509 = ms2Scans[i].MassSpectrum.GetClosestPeakIndex(509.31);
            //    if (!tolerance.Within(ms2Scans[i].MassSpectrum.XArray[idx509], 509.31)) continue;

            //    intensityArray[i] = ms2Scans[i].MassSpectrum.YArray[idx240] + ms2Scans[i].MassSpectrum.YArray[idx509];
            //}

            Chart.Bar<double, double, string>(scan.MassSpectrum.XArray.Select(x => Math.Log(x)), scan.MassSpectrum.YArray)
                .WithTraceInfo("Diagnostic Ions")
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
                .WithSize(Width: 1000, Height: 500)
                .Show();
        }
    }
}
