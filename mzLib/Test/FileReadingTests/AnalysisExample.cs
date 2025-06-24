using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MzIdentML;
using MzLibUtil;
using NUnit.Framework;
using Plotly.NET;
using Plotly.NET.CSharp;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using System.Printing;
using System.Windows.Documents;
using System.Xml.Serialization;
using static Plotly.NET.StyleParam.DrawingStyle;
using Chart = Plotly.NET.CSharp.Chart;
using GenericChartExtensions = Plotly.NET.CSharp.GenericChartExtensions;
using Stopwatch = System.Diagnostics.Stopwatch;
using Peptide = Proteomics.AminoAcidPolymer.Peptide;
using IsotopicEnvelope = MassSpectrometry.IsotopicEnvelope;
using System.Printing.Interop;
using Microsoft.ML.Trainers;

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

            var x = Chart.Line<double, double, string>(rtArray, intensityArray)
                .WithTitle("Diagnostic Ions")
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Retention Time (min)"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity of Diagnostic Ions"))
                .WithSize(Width: 1000, Height: 500);

            GenericChartExtensions.Show(x);
        }

        [Test]
        public static void AnalyzeFlashData()
        {
            string dataPath = @"D:\JurkatTopdown\MM108_DecoyPTMs_MOxVariable\FlashLFQ_Tweaked\QuantifiedPeaks.tsv";
            var peaks = FileReader.ReadFile<QuantifiedPeakFile>(dataPath);
            int mixedDecoyRealModPeaks = 0;
            int successfulDisambiguations = 0;
            int totalDecoyMods = 0;
            List<QuantifiedPeak> realPeaks = new List<QuantifiedPeak>();
            foreach (var peak in peaks)
            {
                if (!peak.FullSequence.Contains(":Decoy") && peak.FileName.Contains("jurkat_td_rep2_fract5") )
                {
                    realPeaks.Add(peak);
                }

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

            

            var mostIntensePeak = realPeaks.First(p => p.FullSequence == "PELAKSAPAPKKGSKKAVTKAQKKDGKKRKRSRKESYSVYVYKVLKQVHPDTGISSKAMGIMNSFVNDIFERIASEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSSK");
            List<double> mzChannels = new List<double>();
            for (int i = 0; i < 5; i++)
            {
                double baseMz = (double)mostIntensePeak.PeakMz;
                double isotopeSpacing = 1.0 / (double)mostIntensePeak.PeakCharge;
                double isotopeDeltaMz = isotopeSpacing * (i - 2);
                // Add the mz for the current isotope channel
                mzChannels.Add(baseMz + isotopeDeltaMz);
            }

            var datafile = MsDataFileReader.GetDataFile(@"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract5.raw");
            var indexingEngine = PeakIndexingEngine.InitializeIndexingEngine(datafile);

            var tol = new PpmTolerance(10);
            List<GenericChart> xicPlots = new List<GenericChart>();
            foreach (var mz in mzChannels)
            {
                var xic = indexingEngine.GetXic(mz, (double)mostIntensePeak.PeakRTApex, tol, missedScansAllowed: 10, maxPeakHalfWidth: 5);
                xicPlots.Add(GetXicChart(xic));
            }

            GenericChartExtensions.Show(Chart.Combine(xicPlots));
        }

        [Test]
        public static void SimulatedXIC()
        {
            var peakShapeAlgorithm = new PeakShapeAlgorithm(frontingFactor: 0.252, tailingFactor: 0.002, heightFactor: 400);
            List<double> retentionTimes = Enumerable.Range(0, 20).Select(i => i * 0.1 + 39.5).ToList(); // Simulated retention times in minutes
            List<double> intensities = peakShapeAlgorithm.GetIntensityRange(40.10, retentionTimes);

            var chart = Chart.Line<double, double, string>(retentionTimes.ToArray(), intensities.ToArray())
                .WithTitle("XIC")
                //.WithLayout(Layout.init<IConvertible>(PlotBGColor: Plotly.NET.Color.fromString("white")))
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("RT (min)"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Relative Abundance"))
                .WithLineStyle(Width: 3)
                //.WithSize(Width: 3000, Height: 1600);
                .WithSize(Width: 750, Height: 400);
            GenericChartExtensions.Show(chart);
        }

        [Test]
        public static void SimulatedXICEmg()
        {
            var peakShapeAlgorithm = new EmgAlgorithm()
            {
                SkewFactor = -0.5, // Negative skew for fronting peak shape 
                Intensity = 400,
                StandardDeviation = 0.5
            };
            List<double> retentionTimes = Enumerable.Range(0, 20).Select(i => i * 0.1 + 39.5).ToList(); // Simulated retention times in minutes
            List<double> intensities = peakShapeAlgorithm.GetIntensityRange(40.10, retentionTimes);

            var chart = Chart.Line<double, double, string>(retentionTimes.ToArray(), intensities.ToArray())
                .WithTitle("XIC")
                //.WithLayout(Layout.init<IConvertible>(PlotBGColor: Plotly.NET.Color.fromString("white")))
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("RT (min)"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Relative Abundance"))
                .WithLineStyle(Width: 3)
                //.WithSize(Width: 3000, Height: 1600);
                .WithSize(Width: 750, Height: 400);
            GenericChartExtensions.Show(chart);
        }

        [Test]
        public static void SimulatedEtg()
        {
            var peakShapeAlgorithm = new EmpiricallyTransformedGaussian()
            {
                FullWidthHalfMax = 0.3,

                LambdaLeading = 3,
                KexpLeading = 5,
                AlphaLeading = 0.5,

                LambdaTailing = 5,
                KexpTailing = 15,
                AlphaTailing = 10
            };
            List<double> retentionTimes = Enumerable.Range(0, 1200).Select(i => i * 0.01 + 35).ToList(); // Simulated retention times in minutes
            List<double> intensities = peakShapeAlgorithm.GetIntensityRange(40.10, retentionTimes);

            var chart = Chart.Line<double, double, string>(retentionTimes.ToArray(), intensities.ToArray())
                .WithTitle("XIC")
                //.WithLayout(Layout.init<IConvertible>(PlotBGColor: Plotly.NET.Color.fromString("white")))
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("RT (min)"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Relative Abundance"))
                .WithLineStyle(Width: 3)
                //.WithSize(Width: 3000, Height: 1600);
                .WithSize(Width: 750, Height: 400);
            GenericChartExtensions.Show(chart);
        }


        public static GenericChart GetXicChart(List<IIndexedPeak> peaks)
        {
            var xarray = peaks.Select(peak => peak.RetentionTime).ToArray();
            var yarray = peaks.Select(peak => peak.Intensity).ToArray();

            return Chart.Line<double, double, string>(xarray, yarray)
                .WithTitle("XIC")
                //.WithLayout(Layout.init<IConvertible>(PlotBGColor: Plotly.NET.Color.fromString("white")))
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("RT (min)"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Relative Abundance"))
                .WithLineStyle(Width: 3)
                //.WithSize(Width: 3000, Height: 1600);
                .WithSize(Width: 750, Height: 400);
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

            //Chart.Line<double, double, string>(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray)
            //    .WithTitle("Theoretical Envelope")
            //    .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
            //    .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
            //    .WithSize(Width: 1000, Height: 500)
            //    .Show();

        }

        [Test]
        public static void ProFormaReader()
        {
            string prosightPath = @"C:\Users\Alex\Downloads\Jurkat_ProsightPdChimeras_Rep2_15_10ppm_PSMs_ProformaFile.tsv";
            string mmPath = @"C:\Users\Alex\Downloads\Jurkat_MetaMorpheus_Rep2_WithLibrary_NewPEP_NoNorm_PSMs_ProformaFile.tsv";
            string mspftPath = @"C:\Users\Alex\Downloads\Jurkat_MsPathFinderTWithMods_15Rep2_Final_PSMs_ProformaFile.tsv";

            var prosightFile = new ProformaFile(prosightPath)
                .Where(id => id.FileName == "2_05" && id.ScanNumber > 2200 && id.ScanNumber < 2320)
                .DistinctBy(id => id.FullSequence)
                .ToList();

            var mmFile = new ProformaFile(mmPath)
                .Where(id => id.FileName == "2_05" && id.ScanNumber > 2200 && id.ScanNumber < 2320)
                .DistinctBy(id => id.FullSequence)
                .ToList();

            //var mspftFile = new ProformaFile(mspftPath)
            //    .Where(id => id.FileName == "2_05" && id.ScanNumber > 2100 && id.ScanNumber < 2420)
            //    .DistinctBy(id => id.FullSequence)
            //    .ToList();

            var sequenceIntersect = prosightFile.Select(p => p.FullSequence).Intersect(mmFile.Select(m => m.FullSequence)).ToHashSet();


            var normal = new Normal(2264, 40);

            List<PeptideWithSetModifications> sharedForms = new();
            List<double> intensities = new();
            List<double> manualScaling = new();

            using(StreamWriter sw = new StreamWriter(@"C:\Users\Alex\Documents\TopDownVision\ConsensusSimulatedProteoforms.tsv"))
            {
                sw.WriteLine("ScanNumber\tAccession\tMost Abundant Mass\tMonoisotopic Mass\tModLocations\tFullSequence\tElutionIntensity");
                foreach (var record in prosightFile)
                {
                    if (!sequenceIntersect.Contains(record.FullSequence)) continue;
                    var pep = new PeptideWithSetModifications(record.FullSequence, ModificationConverter.AllKnownModsDictionary);
                    sharedForms.Add(pep);
                    var dist = IsotopicDistribution.GetDistribution(sharedForms.Last().FullChemicalFormula);
                    double mostAbundant = dist.MostAbundantMass.ToMz(18);
                    intensities.Add(normal.Density(record.ScanNumber) * 100);
                    sw.WriteLine(record.ScanNumber + "\t" +
                                 record.ProteinAccession + "\t" +
                                mostAbundant + "\t" + 
                                pep.MonoisotopicMass.ToMz(18) + "\t" +
                                string.Join(";", pep.AllModsOneIsNterminus.Select(kvp => kvp.Key + "-" + kvp.Value.IdWithMotif)) + "\t" + 
                                record.FullSequence + "\t" + 
                                normal.Density(record.ScanNumber) * 100);
                }
            }
            

            double baseScalingFactor = 75000;

            HashSet<int> alreadyObservedMasses = new();

            List<double[]> mzArrayList = new();
            List<int[]> intensityArrayList = new();
            PopulateArrayLists(mzArrayList, intensityArrayList, sharedForms, intensities, alreadyObservedMasses, baseScalingFactor, customScalingFactor: 0.1);
            var displaySpectrum = TofSpectraMerger.MergeArraysToMs2Spectrum(mzArrayList, intensityArrayList, ppmTolerance: 5);
            SimulatedData.AddHighFrequencyNoiseToYArray(displaySpectrum.YArray, new Normal(250, 250));
            Plotly.NET.GenericChart sharedPlot = GetSpectrumPlot(displaySpectrum.XArray, displaySpectrum.YArray, "green", opacity: 0.95);


            var metaSequences = mmFile.Select(m => m.FullSequence).Where(s => !sequenceIntersect.Contains(s)).ToHashSet();
            sharedForms.Clear();
            intensities.Clear();
            using (StreamWriter sw = new StreamWriter(@"C:\Users\Alex\Documents\TopDownVision\MetaMorpheusSimulatedProteoforms2.tsv"))
            {
                sw.WriteLine("ScanNumber\tAccession\tMost Abundant Mass z=18\tMost Abundant Mass z = 14\tModLocations\tFullSequence\tElutionIntensity");
                foreach (var record in mmFile)
                {
                    if (!metaSequences.Contains(record.FullSequence)) continue;
                    var pep = new PeptideWithSetModifications(record.FullSequence, ModificationConverter.AllKnownModsDictionary);
                    sharedForms.Add(pep);
                    var dist = IsotopicDistribution.GetDistribution(sharedForms.Last().FullChemicalFormula);
                    double mostAbundant = dist.MostAbundantMass.ToMz(18);
                    intensities.Add(normal.Density(record.ScanNumber) * 100);
                    sw.WriteLine(record.ScanNumber + "\t" +
                                 record.ProteinAccession + "\t" +
                                 mostAbundant + "\t" +
                                 dist.MostAbundantMass.ToMz(14) + "\t" +
                                 string.Join(";", pep.AllModsOneIsNterminus.Select(kvp => kvp.Key + "-" + kvp.Value.IdWithMotif)) + "\t" +
                                 record.FullSequence + "\t" +
                                 normal.Density(record.ScanNumber) * 100);
                }
            }

            mzArrayList.Clear();
            intensityArrayList.Clear();
            PopulateArrayLists(mzArrayList, intensityArrayList, sharedForms, intensities, alreadyObservedMasses, baseScalingFactor: 25000, meta: true);
            displaySpectrum = TofSpectraMerger.MergeArraysToMs2Spectrum(mzArrayList, intensityArrayList, ppmTolerance: 5);  
            SimulatedData.AddHighFrequencyNoiseToYArray(displaySpectrum.YArray, new Normal(250, 250));
            Plotly.NET.GenericChart mmPlot = GetSpectrumPlot(displaySpectrum.XArray, displaySpectrum.YArray, "blue", scalingFactor: 0.4);


            var prosightSequences = prosightFile.Select(p => p.FullSequence).Where(s => !sequenceIntersect.Contains(s)).ToHashSet();
            sharedForms.Clear();
            intensities.Clear();

            using (StreamWriter sw = new StreamWriter(@"C:\Users\Alex\Documents\TopDownVision\ProSightSimulatedProteoforms.tsv"))
            {
                sw.WriteLine("ScanNumber\tAccession\tMost Abundant Mass z=18\tMost Abundant Mass z = 15\tModLocations\tFullSequence\tElutionIntensity");
                foreach (var record in prosightFile)
                {
                    if (!prosightSequences.Contains(record.FullSequence)) continue;
                    var pep = new PeptideWithSetModifications(record.FullSequence, ModificationConverter.AllKnownModsDictionary);
                    sharedForms.Add(pep);
                    var dist = IsotopicDistribution.GetDistribution(sharedForms.Last().FullChemicalFormula);
                    double mostAbundant = dist.MostAbundantMass.ToMz(18);
                    intensities.Add(normal.Density(record.ScanNumber) * 100);
                    sw.WriteLine(record.ScanNumber + "\t" +
                                 record.ProteinAccession + "\t" +
                                 mostAbundant + "\t" +
                                 dist.MostAbundantMass.ToMz(15) + "\t" +
                    string.Join(";", pep.AllModsOneIsNterminus.Select(kvp => kvp.Key + "-" + kvp.Value.IdWithMotif)) + "\t" +
                                 record.FullSequence + "\t" +
                                 normal.Density(record.ScanNumber) * 100);
                }
            }

            mzArrayList.Clear();
            intensityArrayList.Clear();
            PopulateArrayLists(mzArrayList, intensityArrayList, sharedForms, intensities, alreadyObservedMasses, baseScalingFactor, prosight: true);
            displaySpectrum = TofSpectraMerger.MergeArraysToMs2Spectrum(mzArrayList, intensityArrayList, ppmTolerance: 5);
            SimulatedData.AddHighFrequencyNoiseToYArray(displaySpectrum.YArray, new Normal(250, 250));
            Plotly.NET.GenericChart prosightPlot = GetSpectrumPlot(displaySpectrum.XArray, displaySpectrum.YArray, "yellow", scalingFactor: 0.13);

            var experimentalSpectraFile = MsDataFileReader.GetDataFile(@"C:\Users\Alex\Documents\TopDownVision\02-18-20_jurkat_td_rep2_fract5-AveragedHistone.raw");
            experimentalSpectraFile.LoadAllStaticData();
            var expSpectrum = experimentalSpectraFile.Scans.First().MassSpectrum;
            var expPlot = GetSpectrumPlot(expSpectrum.XArray, expSpectrum.YArray, color: "black", opacity: 0.8, mirror: true);

            var envelopePlots = GetEnvelopePlots();

            //var combinedChart = Chart.Combine(new[] {   prosightPlot, mmPlot, sharedPlot, expPlot, envelopePlots[0], envelopePlots[1], envelopePlots[2] });
            var combinedChart = Chart.Combine(new[] { sharedPlot, expPlot, envelopePlots[0], envelopePlots[1], envelopePlots[2] });

            GenericChartExtensions.Show(combinedChart);

            var testPep = new PeptideWithSetModifications(prosightFile.First().FullSequence, ModificationConverter.AllKnownModsDictionary);
            var mass = testPep.FullChemicalFormula.MonoisotopicMass;

            /// Histone target scan = 2264 +- 120

            int x = 0;
            //var noiseParams = new LowFrequencyNoiseParameters(peakNumberLimitLow: 5999, peakNumberLimitHigh: 6000,
            //    peakLocationLimitLow: 625, peakLocationLimitHigh: 1000, peakIntensityLimitHigh: 1500);

            //var lowFreqNoiseSpectrum = SimulatedData.BuildLowFrequencyNoiseSpectrum(noiseParams);
            //mzArrayList.Add(lowFreqNoiseSpectrum.XArray);
            //intensityArrayList.Add(lowFreqNoiseSpectrum.YArray.Select(d => (int)d).ToArray());
        }

        [Test]
        public static void GetEnvelopePlotsTest()
        {
            //colors ??= new[] { "#a51c30", "#74121d", "#580c1f"};
            string[] colors = new[] { "#fb747d", "#e30613", "#77030b" };

            var experimentalSpectraFile = MsDataFileReader.GetDataFile(@"C:\Users\Alex\Documents\TopDownVision\02-18-20_jurkat_td_rep2_fract5-AveragedHistone.raw");
            experimentalSpectraFile.LoadAllStaticData();
            var expSpectrum = experimentalSpectraFile.Scans.First().MassSpectrum;

            var scan = new MsDataScan(massSpectrum: expSpectrum, oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: 1, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: expSpectrum.SumOfAllY, injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (1),
                isolationMZ: 763.75, isolationWidth: 2.4);

            DeconvolutionParameters isoDecParams = new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false);

            var precursors = new List<IsotopicEnvelope>();
            precursors.AddRange(scan.GetIsolatedMassesAndCharges(scan.MassSpectrum, isoDecParams).OrderBy(e => e.MostAbundantObservedIsotopicMass));
            var deconvolutedMonoisotopicMasses = precursors.Select(p => p.MonoisotopicMass).OrderBy(d => d).Distinct().ToList();
            var deconvolutedMzs = deconvolutedMonoisotopicMasses.Select(m => m.ToMz(18)).ToList();
            var estimatedMostAbundantMzs = deconvolutedMzs.Select(m => m + 0.445768).ToList();


            List<GenericChart> plots = new List<GenericChart>();
            for (int i = 2; i < 5; i++)
            {
                List<double> xs = new();
                List<double> ys = new();
                for (int j = 0; j < precursors[i].Peaks.Count; j++)
                {
                    if (precursors[i].Peaks[j].mz < 764.9)
                    {
                        xs.Add(precursors[i].Peaks[j].mz);
                        ys.Add(precursors[i].Peaks[j].intensity);
                    }
                }
                //double[] xArray = precursors[i].Peaks.Select(p => p.mz).ToArray();
                plots.Add(GetSpectrumPlot(
                    xs.ToArray(),
                    ys.ToArray(),
                    color: colors[i - 2],
                    opacity: 1,
                    mirror: true,
                    maxIntensity: scan.MassSpectrum.YofPeakWithHighestY));
            }

            //return plots;
        }

        public static List<GenericChart> GetEnvelopePlots(string[] colors = null)
        {
            //colors ??= new[] { "#a51c30", "#74121d", "#580c1f"};
            colors ??= new[] { "#fb747d", "#e30613", "#77030b" };

            var experimentalSpectraFile = MsDataFileReader.GetDataFile(@"C:\Users\Alex\Documents\TopDownVision\02-18-20_jurkat_td_rep2_fract5-AveragedHistone.raw");
            experimentalSpectraFile.LoadAllStaticData();
            var expSpectrum = experimentalSpectraFile.Scans.First().MassSpectrum;

            var scan = new MsDataScan(massSpectrum: expSpectrum, oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: 1, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: expSpectrum.SumOfAllY, injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (1),
                isolationMZ: 763.75, isolationWidth: 2.4);

            DeconvolutionParameters isoDecParams = new IsoDecDeconvolutionParameters(reportMultipleMonoisos: false);

            var precursors = new List<IsotopicEnvelope>();
            precursors.AddRange(scan.GetIsolatedMassesAndCharges(scan.MassSpectrum, isoDecParams).OrderBy(e => e.MostAbundantObservedIsotopicMass));
            var deconvolutedMonoisotopicMasses = precursors.Select(p => p.MonoisotopicMass).OrderBy(d => d).Distinct().ToList();
            var deconvolutedMzs = deconvolutedMonoisotopicMasses.Select(m => m.ToMz(18)).ToList();
            var estimatedMostAbundantMzs = deconvolutedMzs.Select(m => m + 0.445768).ToList();


            List<GenericChart> plots = new List<GenericChart>();
            for (int i = 2; i < 5; i++)
            {
                List<double> xs = new();
                List<double> ys = new();
                for (int j = 0; j < precursors[i].Peaks.Count; j++)
                {
                    if (precursors[i].Peaks[j].mz < 764.9)
                    {
                        xs.Add(precursors[i].Peaks[j].mz);
                        ys.Add(precursors[i].Peaks[j].intensity);
                    }
                } 
                //double[] xArray = precursors[i].Peaks.Select(p => p.mz).ToArray();
                plots.Add(GetSpectrumPlot(
                    xs.ToArray(),
                    ys.ToArray(),
                    color: colors[i-2],
                    opacity: 1,
                    mirror: true,
                    maxIntensity: scan.MassSpectrum.YofPeakWithHighestY));
            }

            return plots;
        }


        public static GenericChart GetSpectrumPlot(double[] xArray, double[] yArray, string color = "red", double opacity = 0.8, bool mirror = false,
            double scalingFactor = 1, bool experimental = false, double? maxIntensity = null)
        {
            //filtering step
            List<double> xList = new();
            List<double> yList =new();
            for (int i = 0; i < yArray.Length; i++)
            {
                if (xArray[i] > 754 & xArray[i] < 776)
                {
                    xList.Add(xArray[i]);
                    yList.Add(yArray[i]);
                }
                
            }
            xArray = xList.ToArray();
            yArray = yList.ToArray();


            double[] plotMz = new double[xArray.Length * 3];
            double[] plotIntensity = new double[yArray.Length * 3];
            maxIntensity ??= yArray.Max();

            for (int i = 0; i < xArray.Length; i++)
            {
                int first = (3 * i);
                int second = (3 * i + 1);
                int third = (3 * i + 2);

                plotMz[first] = xArray[i] - 0.00001;
                plotIntensity[first] = 0;

                plotMz[second] = xArray[i];
                plotIntensity[second] = Math.Max(0, yArray[i]) * 100 * scalingFactor / (double)maxIntensity;
                if (mirror) plotIntensity[second] *= -1;

                plotMz[third] = xArray[i] + 0.00001;
                plotIntensity[third] = 0;
            }

            return Chart.Line<double, double, string>(plotMz, plotIntensity, LineColor: Plotly.NET.Color.fromString(color), Opacity: opacity)
                .WithTitle("Theoretical Envelope")
                .WithLayout(Layout.init<IConvertible>(PlotBGColor: Plotly.NET.Color.fromString("white")))
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Relative Abundance") )
                .WithLineStyle(Width: 3)
                //.WithSize(Width: 3000, Height: 1600);
                .WithSize(Width: 750, Height: 400);
        }


        public static void PopulateArrayLists(List<double[]> mzArrayList, List<int[]> intensityArrayList, List<PeptideWithSetModifications> forms, List<double> intensities, 
            HashSet<int> alreadyObservedMasses, double baseScalingFactor = 75000, double customScalingFactor = 0.05, bool prosight = false, bool meta = false)
        {
            for (int i = 0; i < forms.Count; i++)
            {
                for (int z = 14; z < 23; z++)
                {
                    int chargeStateScalingFactor = 6 - Math.Abs(z - 18);
                    double randomScalingFactor = intensities[i];

                    var dist = IsotopicDistribution.GetDistribution(forms[i].FullChemicalFormula, 0.125, 1e-8);
                    if (z == 14)
                    {
                        if (alreadyObservedMasses.Contains((int)(dist.MostAbundantMass * 1000)))
                        {
                            break;
                        }
                        else alreadyObservedMasses.Add((int)(dist.MostAbundantMass * 1000));
                    }

                    if (z == 18)
                    {
                        double mostAbundant = dist.MostAbundantMass.ToMz(18);
                        if (Math.Abs(mostAbundant - 766.2) < 0.1)
                        {
                            randomScalingFactor *= 2;

                        }
                        else if (Math.Abs(mostAbundant - 765.43) < 0.1)
                        {
                            randomScalingFactor *= 0.6;
                        }
                        else if (Math.Abs(mostAbundant - 766.31) < 0.1)
                        {
                            randomScalingFactor *= 1;
                        }
                        else if (Math.Abs(mostAbundant - 767.04) < 0.1)
                        {
                            randomScalingFactor *= 1;
                        }
                        else if (Math.Abs(mostAbundant - 767.87) < 0.1)
                        {
                            randomScalingFactor *= 0.5;
                        }
                        else if (Math.Abs(mostAbundant - 767.59) < 0.1)
                        {
                            randomScalingFactor *= 0.2;
                        }
                        else if (meta || prosight)
                        {
                            
                            randomScalingFactor *= customScalingFactor/2;
                        }
                        else
                        {
                            randomScalingFactor *= customScalingFactor;
                        }
                    }

                    if (z == 15)
                    {
                        double mostAbundant = dist.MostAbundantMass.ToMz(15);
                        if (Math.Abs(mostAbundant - 754.7) < 0.1)
                        {
                            randomScalingFactor *= 0.5;

                        }
                        else if (Math.Abs(mostAbundant - 755.63) < 0.1)
                        {
                            randomScalingFactor *= 0.1;
                        }
                        else if (Math.Abs(mostAbundant - 756.56) < 0.1)
                        {
                            randomScalingFactor *= 0.2;
                        }
                        else if (Math.Abs(mostAbundant - 756.36) < 0.1)
                        {
                            randomScalingFactor *= 0.2;
                        }
                        else if (prosight && Math.Abs(mostAbundant - 757.52) < 0.1)
                        {
                            randomScalingFactor *= 0.1;
                        }
                        else if (Math.Abs(mostAbundant - 760.00) < 0.1)
                        {
                            randomScalingFactor *= 0.5;
                        }
                    }

                    if (z == 14)
                    {
                        randomScalingFactor *= 0.3;
                    }

                    // unmodified
                    mzArrayList.Add(dist.Masses.Select(v => v.ToMz(z)).ToArray());
                    intensityArrayList.Add(dist.Intensities.Select(v => (int)(v * baseScalingFactor * chargeStateScalingFactor * randomScalingFactor)).ToArray());
                }
            }
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

            //Chart.Line<double, double, string>(plotMz, plotIntensity)
            //    .WithTitle("Theoretical Envelope")
            //    .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
            //    .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
            //    .WithSize(Width: 1000, Height: 500)
            //    .Show();
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

            //Chart.Point<double, double, string>(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray)
            //    .WithTitle("Diagnostic Ions")
            //    .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
            //    .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
            //    .WithSize(Width: 1000, Height: 500)
            //    .Show();
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

            //Chart.Bar<double, double, string>(scan.MassSpectrum.XArray.Select(x => Math.Log(x)), scan.MassSpectrum.YArray)
            //    .WithTitle("Diagnostic Ions")
            //    .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
            //    .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
            //    .WithSize(Width: 1000, Height: 500)
            //    .Show();
        }
    }
}
