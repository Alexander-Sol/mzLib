using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Readers;
using Readers.QuantificationResults;
using MzLibUtil;
using FlashLFQ;
using Easy.Common.Extensions;
using Plotly.NET.CSharp;
using Plotly.NET;
using Chart = Plotly.NET.CSharp.Chart;
using FSharp.Data;
using System.Net;
using Plotly.NET.TraceObjects;
using Plotly.NET.LayoutObjects;
using static Deedle.Vectors.VectorConstruction;
using Chemistry;
using System.Reflection;
using MathNet.Numerics;



namespace Test.FileReadingTests
{
    [ExcludeFromCodeCoverage]
    internal class TestQuantifiedPeak
    {
        internal static string TestDirectory;
        internal static string TestFilePath;

        [OneTimeSetUp]
        public void SetUp()
        {
            TestDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ReadingWritingTests");
            TestFilePath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                @"FileReadingTests\ExternalFileTypes\FlashLFQ_MzLib1.0.549_QuantifiedPeaks.tsv");
            Directory.CreateDirectory(TestDirectory);
        }

        [OneTimeTearDown]
        public void TearDown()
        {
            Directory.Delete(TestDirectory, true);
        }

        [Test]
        public static void TestFileLoadsAndCountCorrect()
        {
            QuantifiedPeakFile file = new QuantifiedPeakFile(TestFilePath);
            Assert.That(file.Count(), Is.EqualTo(4));
            Assert.That(file.CanRead(TestFilePath));

            file = FileReader.ReadFile<QuantifiedPeakFile>(TestFilePath);
            Assert.That(file.Count(), Is.EqualTo(4));
            Assert.That(file.CanRead(TestFilePath));
        }

        [Test]
        public static void LocalDataAnalysis()
        {
            string inhousePath = @"D:\Human_Ecoli_TwoProteome_60minGradient\Human_FlashLFQ_326_DonorPepQ_0pt1\QuantifiedPeaks.tsv";
            QuantifiedPeakFile file = new QuantifiedPeakFile(inhousePath);

            var msmsPeaks = file.Where(peak => peak.FileName.Contains("Human_C18"))
                .GroupBy(peak => peak.FileName)
                .MaxBy(group => group.Count(peak => peak.PeakDetectionType == "MSMS"))
                .Where(peak => peak.PeakMz != null && peak.PeakRTApex != null && peak.PeakDetectionType == "MSMS")
                .OrderBy(peak => peak.PeakMz)
                .ToList();

            Tolerance ppmTolerance = new PpmTolerance(10);
            List<List<QuantifiedPeak>> coelutingPeaks = new List<List<QuantifiedPeak>>();

            HashSet<string> peptideSeqs = msmsPeaks.Select(msmsPeaks => msmsPeaks.FullSequence).ToHashSet();
            
            for (int i = 0; i < msmsPeaks.Count; i++)
            {
                int j = i+1;
                List<QuantifiedPeak> peakGroup = new List<QuantifiedPeak>();
                while (j < msmsPeaks.Count)
                {
                    var peak = msmsPeaks[i];
                    var nextPeak = msmsPeaks[j];
                    if (peak.BaseSequence == nextPeak.BaseSequence)
                    {
                        j++;
                        continue;
                    }
                    
                    if (ppmTolerance.Within((double)peak.PeakMz, (double)nextPeak.PeakMz))
                    {
                        j++;
                        if(Math.Abs((double)peak.PeakRTApex - (double)nextPeak.PeakRTApex) < 0.397
                            && peak.PeakCharge == nextPeak.PeakCharge)
                        {
                            peakGroup.Add(peak);
                            peakGroup.Add(nextPeak);
                        }
                    }
                    else
                    {
                        if(peakGroup.IsNotNullOrEmpty())
                            coelutingPeaks.Add(peakGroup);
                        break;
                    }
                }
            }

            int totalCoelutingPeaks = coelutingPeaks.Sum(group => group.Count);

            List<double> rtArray = coelutingPeaks.SelectMany(group => group.Select(peak => (double)peak.PeakRTApex)).ToList();
            List<double> mzArray = coelutingPeaks.SelectMany(coelutingPeaks => coelutingPeaks.Select(peak => (double)peak.PeakMz)).ToList();

            var pairedScatter = Chart.Point<double, double, string>(x: rtArray, y: mzArray,
                MarkerColor: Color.fromKeyword(ColorKeyword.Red), Name: "Co-Eluting Peaks",
                Opacity: 0.75);

            List<double> rtArray2 = msmsPeaks.Select(peak => (double)peak.PeakRTApex).ToList();
            List<double> mzArray2 = msmsPeaks.Select(peak => (double)peak.PeakMz).ToList();

            var allScatter = Chart.Point<double, double, string>(rtArray2, mzArray2, Name: "All Peaks",
                MarkerColor: Color.fromKeyword(ColorKeyword.Blue),
                Opacity: 0.5);

            var combined = Chart.Combine([allScatter, pairedScatter])
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Retention Time (min)"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
                .WithSize(Width: 1000, Height: 500);

            Plotly.NET.CSharp.GenericChartExtensions.Show(combined);

            Assert.That(file.Count(), Is.EqualTo(4));
            Assert.That(file.CanRead(TestFilePath));

            file = FileReader.ReadFile<QuantifiedPeakFile>(TestFilePath);
            Assert.That(file.Count(), Is.EqualTo(4));
            Assert.That(file.CanRead(TestFilePath));
        }

        [Test]
        public static void NperVizInhouse()
        {
            string inhousePath = @"D:\Human_Ecoli_TwoProteome_60minGradient\Human_FlashLFQ_326_DonorPepQ_0pt1\QuantifiedPeaks.tsv";
            string censoredPath = @"D:\Human_Ecoli_TwoProteome_60minGradient\CensoredHuman_FlashLFQ_323_DonorPepQ_0pt1\QuantifiedPeaks.tsv";

            QuantifiedPeakFile file = new QuantifiedPeakFile(inhousePath);
            QuantifiedPeakFile censoredFile = new QuantifiedPeakFile(censoredPath);



            var originalFileMsmsPeaks = file.Where(peak => peak.FileName.Contains("Human_C18"))
                .GroupBy(peak => peak.FileName)
                .OrderByDescending(group => group.Count(peak => peak.PeakDetectionType == "MSMS"))
                .Skip(1).Take(1)
                .SelectMany(group => group.Where(peak => peak.PeakMz != null && peak.PeakRTApex != null))
                .OrderBy(peak => peak.PeakMz)
                .ToList();

            string dataFile = originalFileMsmsPeaks.First().FileName;

            int distinctPeptides = originalFileMsmsPeaks.Select(peak => peak.FullSequence).Distinct().Count();

            var multiPeaks = originalFileMsmsPeaks.GroupBy(peak => peak.FullSequence)
                .Where(group => group.Count() > 1);
                

            string censoredPsms = @"D:\Human_Ecoli_TwoProteome_60minGradient\CensoredHumanData\CensoredPsms.psmtsv";
            HashSet<string> censoredPeakSeqs = new();
            using (StreamReader reader = new(censoredPsms))
            {
                while (!reader.EndOfStream)
                {
                    var line = reader.ReadLine().Split('\t');
                    if (line[0] == dataFile)
                    {
                        censoredPeakSeqs.Add(line[9]);
                    }
                }
            }

            var originalCensoredPeaks = originalFileMsmsPeaks.Where(peak => censoredPeakSeqs.Contains(peak.FullSequence));

            var censoredFileMbrPeaks = censoredFile
                .Where(peak => peak.PeakMz != null && peak.PeakRTApex != null && peak.FileName == originalFileMsmsPeaks.First().FileName + "-censored" 
                    && censoredPeakSeqs.Contains(peak.FullSequence) && !peak.RandomRt)
                .OrderBy(peak => peak.PeakMz)
                .ToList();
                

            var peakDict = originalCensoredPeaks
                .GroupBy(peak => peak.FullSequence)
                .ToDictionary(group => group.Key, group => group);


            List<(QuantifiedPeak mbrPeak, QuantifiedPeak msmsPeak)> peakPairs = new();
            foreach(var peak in censoredFileMbrPeaks)
            {
                if(peakDict.ContainsKey(peak.FullSequence))
                {
                    QuantifiedPeak msmsPeak = peakDict[peak.FullSequence]
                        .MinBy(msPeak => Math.Abs((double)msPeak.PeakRTApex - (double)peak.PeakRTApex));
                    if (msmsPeak.PeakCharge != peak.PeakCharge) continue;
                    peak.PeakRtDiff = (double)peak.PeakRTApex - (double)msmsPeak.PeakRTApex;
                    if( (peak.PeakRtDiff <= 0.39 || peak.PeakRtDiff >= -0.39)
                        && Math.Abs(peak.PeakRtDiff) > 0.05)
                    {                         
                        peakPairs.Add((peak, msmsPeak));
                    }
                }
                    
            }

            string mzmlDirectory = @"D:\Human_Ecoli_TwoProteome_60minGradient\CalibrateSearch_4_19_24\Human_Calibrated_Files";
            string mzmlPath = Path.Combine(mzmlDirectory, dataFile + ".mzML");

            peakPairs.RemoveAt(0);

            PeakIndexingEngine indexingEngine = new();
            SpectraFileInfo fileInfo = new(mzmlPath, "A", 1, 1, 1);
            Dictionary<SpectraFileInfo, Ms1ScanInfo[]> ms1Scans = new();

            var test = indexingEngine.IndexMassSpectralPeaks(fileInfo, true, ms1Scans);

            foreach(var pair in peakPairs)
            {
                var combinedPlot = GetXicPlot(pair.mbrPeak, pair.msmsPeak, indexingEngine, ms1Scans[fileInfo]);

                Plotly.NET.CSharp.GenericChartExtensions.Show(combinedPlot);
            }

        }

        [Test]
        public static void NperVizGygi()
        {
            string gygiPath = @"D:\GygiTwoProteome_PXD014415\FlashLFQ_326_DonorPepQ_0pt1\QuantifiedPeaks.tsv";
            string censoredPath = @"D:\GygiTwoProteome_PXD014415\CensoredData_FlashLFQ_326_DonorPepQ_0pt1\QuantifiedPeaks.tsv";

            QuantifiedPeakFile file = new QuantifiedPeakFile(gygiPath);
            QuantifiedPeakFile censoredFile = new QuantifiedPeakFile(censoredPath);



            var originalFileMsmsPeaks = file.Where(peak => peak.FileName.Contains("uman_90"))
                .GroupBy(peak => peak.FileName)
                .OrderByDescending(group => group.Count(peak => peak.PeakDetectionType == "MSMS"))
                .Skip(4)
                .Take(1)
                .SelectMany(group => group.Where(peak => peak.PeakMz != null && peak.PeakRTApex != null))
                .OrderBy(peak => peak.PeakMz)
                .ToList();

            string dataFile = originalFileMsmsPeaks.First().FileName;

            int distinctPeptides = originalFileMsmsPeaks.Select(peak => peak.FullSequence).Distinct().Count();

            var multiPeaks = originalFileMsmsPeaks.GroupBy(peak => peak.FullSequence)
                .Where(group => group.Count() > 1);


            string censoredPsms = @"D:\GygiTwoProteome_PXD014415\MsConvertmzMLs\MM105_7_17_24\MetaMorpheusCensoredmzMLs\CensoredPsms.psmtsv";
            HashSet<string> censoredPeakSeqs = new();
            using (StreamReader reader = new(censoredPsms))
            {
                while (!reader.EndOfStream)
                {
                    var line = reader.ReadLine().Split('\t');
                    if (line[0] == dataFile)
                    {
                        censoredPeakSeqs.Add(line[9]);
                    }
                }
            }

            var originalCensoredPeaks = originalFileMsmsPeaks.Where(peak => censoredPeakSeqs.Contains(peak.FullSequence));

            var censoredFileMbrPeaks = censoredFile
                .Where(peak => peak.PeakMz != null && peak.PeakRTApex != null && peak.FileName == originalFileMsmsPeaks.First().FileName + "-censored"
                    && censoredPeakSeqs.Contains(peak.FullSequence) && !peak.RandomRt && peak.PeakDetectionType == "MBR")
                .OrderBy(peak => peak.PeakMz)
                .ToList();


            var peakDict = originalCensoredPeaks
                .GroupBy(peak => peak.FullSequence)
                .ToDictionary(group => group.Key, group => group);


            List<(QuantifiedPeak mbrPeak, QuantifiedPeak msmsPeak)> peakPairs = new();
            foreach (var peak in censoredFileMbrPeaks)
            {
                if (peakDict.ContainsKey(peak.FullSequence))
                {
                    QuantifiedPeak msmsPeak = peakDict[peak.FullSequence]
                        .MinBy(msPeak => Math.Abs((double)msPeak.PeakRTApex - (double)peak.PeakRTApex));
                    if (msmsPeak.PeakCharge != peak.PeakCharge) continue;
                    peak.PeakRtDiff = (double)peak.PeakRTApex - (double)msmsPeak.PeakRTApex;
                    if (Math.Abs(peak.PeakRtDiff) > 0.05)
                    {
                        peakPairs.Add((peak, msmsPeak));
                    }
                }

            }

            string mzmlDirectory = @"D:\GygiTwoProteome_PXD014415\MsConvertmzMLs\MM105_7_17_24\Task1-CalibrateTask";
            string mzmlPath = Path.Combine(mzmlDirectory, dataFile + ".mzML");

            peakPairs.RemoveAt(0);

            PeakIndexingEngine indexingEngine = new();
            SpectraFileInfo fileInfo = new(mzmlPath, "A", 1, 1, 1);
            Dictionary<SpectraFileInfo, Ms1ScanInfo[]> ms1Scans = new();

            var test = indexingEngine.IndexMassSpectralPeaks(fileInfo, true, ms1Scans);

            foreach (var pair in peakPairs)
            {
                var combinedPlot = GetXicPlot(pair.mbrPeak, pair.msmsPeak, indexingEngine, ms1Scans[fileInfo]);

                Plotly.NET.CSharp.GenericChartExtensions.Show(combinedPlot);
            }

        }

        public static GenericChart GetXicPlot(QuantifiedPeak mbrPeak, QuantifiedPeak msmsPeak, PeakIndexingEngine indexingEngine, Ms1ScanInfo[] scans)
        {
            int precursorScanIndex = -1;
            int lastScanIndex = -1;
            double startRetentionTime = Math.Min((double)mbrPeak.PeakRTStart, (double)msmsPeak.PeakRTStart) - 1;
            double endRetentionTime = Math.Max((double)mbrPeak.PeakRTEnd, (double)msmsPeak.PeakRTEnd) + 1;
            foreach (Ms1ScanInfo ms1Scan in scans)
            {
                if (precursorScanIndex < 0 && ms1Scan.RetentionTime > startRetentionTime)
                {
                    precursorScanIndex = ms1Scan.ZeroBasedMs1ScanIndex;
                }
                else if (ms1Scan.RetentionTime > endRetentionTime)
                {
                    lastScanIndex = ms1Scan.ZeroBasedMs1ScanIndex;
                    break;
                }
            }

            PpmTolerance tolerance = new(10);
            double peakFindingMz = (double)mbrPeak.PeakMz;
            int peakCharge = (int)mbrPeak.PeakCharge;

            List<IndexedMassSpectralPeak> indexedPeaks = new();
            List<IndexedMassSpectralPeak> indexedPeaksPlus1 = new();
            List<IndexedMassSpectralPeak> indexedPeaksPlus2 = new();
            double peakFindingMass = peakFindingMz.ToMass(peakCharge);
            double peakFindingMassPlus1 = peakFindingMz.ToMass(peakCharge) + Chemistry.Constants.C13MinusC12;
            double peakFindingMassPlus2 = peakFindingMz.ToMass(peakCharge) + 2 * Chemistry.Constants.C13MinusC12;
            for (int i = precursorScanIndex; i < lastScanIndex; i++)
            {
                var imsPeak = indexingEngine.GetIndexedPeak(peakFindingMass, i, tolerance, peakCharge);
                var imsPeakPlus1 = indexingEngine.GetIndexedPeak(peakFindingMassPlus1, i, tolerance, peakCharge);
                indexedPeaks.Add(imsPeak);
                indexedPeaksPlus1.Add(imsPeakPlus1);
                indexedPeaksPlus2.Add(indexingEngine.GetIndexedPeak(peakFindingMassPlus2, i, tolerance, peakCharge));
            }

            var maxIntensity = indexedPeaks.Where(p => p != null).Concat( indexedPeaksPlus1.Where(p => p != null)).Concat(indexedPeaksPlus2.Where(p => p != null)).Max(p => p.Intensity);

            var xicPlot = GetXicTrace(indexedPeaks);
            var xicPlotPlus1 = GetXicTrace(indexedPeaksPlus1);
            var xicPlotPlus2 = GetXicTrace(indexedPeaksPlus2);
            var mbrApex = GetApexPlot(mbrPeak, indexedPeaks, indexedPeaksPlus1, indexedPeaksPlus2);
            var msmsApex = GetApexPlot(msmsPeak, indexedPeaks, indexedPeaksPlus1, indexedPeaksPlus2);

            var msmsStart = GetVerticalTrace((double)msmsPeak.PeakRTStart, maxIntensity / 1.5, "MSMS Start/End");
            var msmsEnd = GetVerticalTrace((double)msmsPeak.PeakRTEnd, maxIntensity / 1.5, "MSMS Start/End");
            var msmsScan = GetVerticalTrace((double)msmsPeak.MS2RetentionTime, maxIntensity / 1.25, "MSMS Scan");

            var mbrStart = GetVerticalTrace((double)mbrPeak.PeakRTStart, maxIntensity / 2.0, "MBR Start/End");
            var mbrEnd = GetVerticalTrace((double)mbrPeak.PeakRTEnd, maxIntensity / 2.0, "MBR Start/End");

            string title = "XIC for " + mbrPeak.BaseSequence + "\n" + mbrPeak.PeakRtDiff.Round(2).ToString() + " minute difference in retention time";

            var combined = Chart.Combine([xicPlot, xicPlotPlus1, xicPlotPlus2, mbrApex, mbrStart, mbrEnd, msmsApex, msmsStart, msmsScan, msmsEnd])
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Retention Time (min)"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
                .WithTitle(title)
                .WithSize(Width: 1000, Height: 500);


            return combined;
        }

        public static GenericChart GetXicTrace(IEnumerable<IndexedMassSpectralPeak> indexedPeaks)
        {
            var intensityArray = indexedPeaks.Where(peak => peak != null).Select(peak => peak.Intensity).ToArray();
            var rtArray = indexedPeaks.Where(peak => peak != null).Select(peak => peak.RetentionTime).ToArray();

            string mzString = indexedPeaks.First(peak => peak != null).Mz.Round(2).ToString();

            var xicPlot = Chart.Line<double, double, string>(rtArray, intensityArray, Name: mzString + " m/z")
            .WithSize(Width: 1000, Height: 500);

            return xicPlot;
        }

        public static GenericChart GetVerticalTrace(double retentionTime, double intensity, string name)
        {
            var color = name.Contains("MBR") ? Color.fromKeyword(ColorKeyword.Purple) : Color.fromKeyword(ColorKeyword.Orange);
            if(name.Contains("Scan"))
            {
                color = Color.fromKeyword(ColorKeyword.Black);
            }

            var vlinePlot = Chart.Line<double, double, string>(new double[] { retentionTime, retentionTime }, new double[] { 0, intensity }, Name: name,
                MarkerColor: color, Opacity: 0.5)
                .WithSize(Width: 1000, Height: 500);
            return vlinePlot;
        }

        public static GenericChart GetApexPlot(QuantifiedPeak qPeak, 
            IEnumerable<IndexedMassSpectralPeak> indexedPeaks,
            IEnumerable<IndexedMassSpectralPeak> indexedPeaks1,
            IEnumerable<IndexedMassSpectralPeak> indexedPeaks2)
        {
            double apexRt = (double)qPeak.PeakRTApex;
            var peak = indexedPeaks.Where(p => p!= null).FirstOrDefault(p => Math.Abs(p.RetentionTime - apexRt) < 0.0001);
            var peak1 = indexedPeaks1.Where(p => p != null).FirstOrDefault(p => Math.Abs(p.RetentionTime - apexRt) < 0.0001);
            var peak2 = indexedPeaks2.Where(p => p != null).FirstOrDefault(p => Math.Abs(p.RetentionTime - apexRt) < 0.0001);
            
            List<double> intensities = new();
            if (peak != null)
            {
                intensities.Add(peak.Intensity);
            }
            if (peak1 != null)
            {
                intensities.Add(peak1.Intensity);
            }
            if (peak2 != null)
            {
                intensities.Add(peak2.Intensity);
            }

            return getApexPlot(apexRt, intensities, qPeak.PeakDetectionType + " Apex");
        }

        private static GenericChart getApexPlot(double retentionTime, List<double> intensities, string name)
        {
            var color = name.Contains("MBR") ? Color.fromKeyword(ColorKeyword.Purple) : Color.fromKeyword(ColorKeyword.Orange);

            double[] rts = new double[intensities.Count];
            for(int i = 0; i < rts.Length; i++)
            {
                rts[i] = retentionTime;
            }
            var apex = Chart.Point<double, double, string>(rts, intensities , Name: name, MarkerColor: color)
                .WithSize(Width: 1000, Height: 500);
            return apex;
        }

        [Test]
        public static void TestFileFirstAndLastAreCorrect()
        {
            QuantifiedPeakFile file = new QuantifiedPeakFile(TestFilePath);
            var first = file.First();

            Assert.That(first.FileName, Is.EqualTo("20100721_Velos1_TaGe_SA_A549_06-calib-averaged"));
            Assert.That(first.BaseSequence, Is.EqualTo("DICNDVLSLLEK"));
            Assert.That(first.FullSequence, Is.EqualTo("DIC[Common Fixed:Carbamidomethyl on C]NDVLSLLEK"));
            Assert.That(first.ProteinGroup, Is.EqualTo("P63104"));
            Assert.That(first.PeptideMonoisotopicMass, Is.EqualTo(1417.712281191));
            Assert.That(first.MS2RetentionTime, Is.EqualTo(188.04623));
            Assert.That(first.PrecursorCharge, Is.EqualTo(3));
            Assert.That(first.TheoreticalMZ, Is.EqualTo(473.578036863879));
            Assert.That(first.PeakIntensity, Is.EqualTo(61339740.974156484));
            Assert.That(first.PeakRTStart, Is.EqualTo(187.71813));
            Assert.That(first.PeakRTApex, Is.EqualTo(188.27129333333335));
            Assert.That(first.PeakRTEnd, Is.EqualTo(195.134625));
            Assert.That(first.PeakMz, Is.EqualTo(709.8649279277056));
            Assert.That(first.PeakCharge, Is.EqualTo(2));
            Assert.That(first.NumChargeStatesObserved, Is.EqualTo(3));
            Assert.That(first.PeakDetectionType, Is.EqualTo("MSMS"));
            Assert.That(first.MBRScore, Is.EqualTo(0));
            Assert.That(first.PSMsMapped, Is.EqualTo(2));
            Assert.That(first.BaseSequencesMapped, Is.EqualTo(1));
            Assert.That(first.FullSequencesMapped, Is.EqualTo(1));
            Assert.That(first.PeakSplitValleyRT, Is.EqualTo(0));
            Assert.That(first.PeakApexMassError, Is.EqualTo(2.1314131880687888));

            var last = file.Last();
            Assert.That(last.FileName, Is.EqualTo("20101230_Velos1_TaGe_SA_Jurkat6-calib-averaged"));
            Assert.That(last.BaseSequence, Is.EqualTo("QDLEAQIRGLREEVEK"));
            Assert.That(last.FullSequence, Is.EqualTo("QDLEAQIRGLREEVEK"));
            Assert.That(last.ProteinGroup, Is.EqualTo("A1A5D9"));
            Assert.That(last.PeptideMonoisotopicMass, Is.EqualTo(1912.001403494));
            Assert.That(last.MS2RetentionTime, Is.EqualTo(71.64922));
            Assert.That(last.PrecursorCharge, Is.EqualTo(4));
            Assert.That(last.TheoreticalMZ, Is.EqualTo(479.007627340379));
            Assert.That(last.PeakIntensity, Is.EqualTo(0));
            Assert.That(last.PeakRTStart, Is.Null);
            Assert.That(last.PeakRTApex, Is.Null);
            Assert.That(last.PeakRTEnd, Is.Null);
            Assert.That(last.PeakMz, Is.Null);
            Assert.That(last.PeakCharge, Is.Null);
            Assert.That(last.NumChargeStatesObserved, Is.EqualTo(0));
            Assert.That(last.PeakDetectionType, Is.EqualTo("MSMS"));
            Assert.That(last.MBRScore, Is.EqualTo(0));
            Assert.That(last.PSMsMapped, Is.EqualTo(1));
            Assert.That(last.BaseSequencesMapped, Is.EqualTo(1));
            Assert.That(last.FullSequencesMapped, Is.EqualTo(1));
            Assert.That(last.PeakSplitValleyRT, Is.EqualTo(0));
            Assert.That(last.PeakApexMassError, Is.EqualTo(double.NaN));
        }

        [Test]
        public static void TestFileReadWrite_WithoutExtensionInPath()
        {
            var file = FileReader.ReadFile<QuantifiedPeakFile>(TestFilePath);
            var testOutputPath = Path.Combine(TestDirectory, "TestOutput");

            file.WriteResults(testOutputPath);
            var newPath = testOutputPath + file.FileType.GetFileExtension();
            Assert.That(File.Exists(newPath));

            var writtenFile = new QuantifiedPeakFile(newPath);
            Assert.That(file.Count(), Is.EqualTo(writtenFile.Count()));

            for (int i = 0; i < file.Count(); i++)
            {
                var originalPeak = file.Results[i];
                var writtenPeak = writtenFile.Results[i];
                Assert.That(originalPeak.FileName, Is.EqualTo(writtenPeak.FileName));
                Assert.That(originalPeak.BaseSequence, Is.EqualTo(writtenPeak.BaseSequence));
                Assert.That(originalPeak.FullSequence, Is.EqualTo(writtenPeak.FullSequence));
                Assert.That(originalPeak.ProteinGroup, Is.EqualTo(writtenPeak.ProteinGroup));
                Assert.That(originalPeak.PeptideMonoisotopicMass, Is.EqualTo(writtenPeak.PeptideMonoisotopicMass).Within(0.0000001));
                Assert.That(originalPeak.MS2RetentionTime, Is.EqualTo(writtenPeak.MS2RetentionTime).Within(0.0000001));
                Assert.That(originalPeak.PrecursorCharge, Is.EqualTo(writtenPeak.PrecursorCharge));
                Assert.That(originalPeak.TheoreticalMZ, Is.EqualTo(writtenPeak.TheoreticalMZ).Within(0.0000001));
                Assert.That(originalPeak.PeakIntensity, Is.EqualTo(writtenPeak.PeakIntensity).Within(0.0000001));
                Assert.That(originalPeak.PeakRTStart, Is.EqualTo(writtenPeak.PeakRTStart).Within(0.0000001));
                Assert.That(originalPeak.PeakRTApex, Is.EqualTo(writtenPeak.PeakRTApex).Within(0.0000001));
                Assert.That(originalPeak.PeakRTEnd, Is.EqualTo(writtenPeak.PeakRTEnd).Within(0.0000001));
                Assert.That(originalPeak.PeakMz, Is.EqualTo(writtenPeak.PeakMz).Within(0.0000001));
                Assert.That(originalPeak.PeakCharge, Is.EqualTo(writtenPeak.PeakCharge));
                Assert.That(originalPeak.NumChargeStatesObserved, Is.EqualTo(writtenPeak.NumChargeStatesObserved));
                Assert.That(originalPeak.PeakDetectionType, Is.EqualTo(writtenPeak.PeakDetectionType));
                Assert.That(originalPeak.MBRScore, Is.EqualTo(writtenPeak.MBRScore).Within(0.0000001));
                Assert.That(originalPeak.PSMsMapped, Is.EqualTo(writtenPeak.PSMsMapped));
                Assert.That(originalPeak.BaseSequencesMapped, Is.EqualTo(writtenPeak.BaseSequencesMapped));
                Assert.That(originalPeak.FullSequencesMapped, Is.EqualTo(writtenPeak.FullSequencesMapped));
                Assert.That(originalPeak.PeakSplitValleyRT, Is.EqualTo(writtenPeak.PeakSplitValleyRT));
                Assert.That(originalPeak.PeakApexMassError, Is.EqualTo(writtenPeak.PeakApexMassError).Within(0.0000001));
            }
        }

        [Test]
        public static void TestFileReadWrite_WithExtensionInPath()
        {
            var file = FileReader.ReadFile<QuantifiedPeakFile>(TestFilePath);
            var testOutputPath = Path.Combine(TestDirectory, "TestOutput_QuantifiedPeaks.tsv");

            file.WriteResults(testOutputPath);
            Assert.That(File.Exists(testOutputPath));

            var writtenFile = new QuantifiedPeakFile(testOutputPath);
            Assert.That(file.Count(), Is.EqualTo(writtenFile.Count()));

            for (int i = 0; i < file.Count(); i++)
            {
                var originalPeak = file.Results[i];
                var writtenPeak = writtenFile.Results[i];
                Assert.That(originalPeak.FileName, Is.EqualTo(writtenPeak.FileName));
                Assert.That(originalPeak.BaseSequence, Is.EqualTo(writtenPeak.BaseSequence));
                Assert.That(originalPeak.FullSequence, Is.EqualTo(writtenPeak.FullSequence));
                Assert.That(originalPeak.ProteinGroup, Is.EqualTo(writtenPeak.ProteinGroup));
                Assert.That(originalPeak.PeptideMonoisotopicMass, Is.EqualTo(writtenPeak.PeptideMonoisotopicMass).Within(0.0000001));
                Assert.That(originalPeak.MS2RetentionTime, Is.EqualTo(writtenPeak.MS2RetentionTime).Within(0.0000001));
                Assert.That(originalPeak.PrecursorCharge, Is.EqualTo(writtenPeak.PrecursorCharge));
                Assert.That(originalPeak.TheoreticalMZ, Is.EqualTo(writtenPeak.TheoreticalMZ).Within(0.0000001));
                Assert.That(originalPeak.PeakIntensity, Is.EqualTo(writtenPeak.PeakIntensity).Within(0.0000001));
                Assert.That(originalPeak.PeakRTStart, Is.EqualTo(writtenPeak.PeakRTStart).Within(0.0000001));
                Assert.That(originalPeak.PeakRTApex, Is.EqualTo(writtenPeak.PeakRTApex).Within(0.0000001));
                Assert.That(originalPeak.PeakRTEnd, Is.EqualTo(writtenPeak.PeakRTEnd).Within(0.0000001));
                Assert.That(originalPeak.PeakMz, Is.EqualTo(writtenPeak.PeakMz).Within(0.0000001));
                Assert.That(originalPeak.PeakCharge, Is.EqualTo(writtenPeak.PeakCharge));
                Assert.That(originalPeak.NumChargeStatesObserved, Is.EqualTo(writtenPeak.NumChargeStatesObserved));
                Assert.That(originalPeak.PeakDetectionType, Is.EqualTo(writtenPeak.PeakDetectionType));
                Assert.That(originalPeak.MBRScore, Is.EqualTo(writtenPeak.MBRScore).Within(0.0000001));
                Assert.That(originalPeak.PSMsMapped, Is.EqualTo(writtenPeak.PSMsMapped));
                Assert.That(originalPeak.BaseSequencesMapped, Is.EqualTo(writtenPeak.BaseSequencesMapped));
                Assert.That(originalPeak.FullSequencesMapped, Is.EqualTo(writtenPeak.FullSequencesMapped));
                Assert.That(originalPeak.PeakSplitValleyRT, Is.EqualTo(writtenPeak.PeakSplitValleyRT));
                Assert.That(originalPeak.PeakApexMassError, Is.EqualTo(writtenPeak.PeakApexMassError).Within(0.0000001));
            }
        }
    }
}
