using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MzIdentML;
using MzLibUtil;
using NUnit.Framework;
using Plotly.NET.CSharp;
using Proteomics.AminoAcidPolymer;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using System.Xml.Serialization;
using Test.FileReadingTests;
using TopDownSimulator.Extraction;
using TopDownSimulator.Fitting;
using TopDownSimulator.Model;
using TopDownSimulator.Noise;
using TopDownSimulator.Simulation;
using UsefulProteomicsDatabases;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test.FileReadingTests
{
    [TestFixture]
    internal class AnalysisExample
    {
        private enum SimulationRunMode
        {
            QuickDev,
            FullFidelity,
        }

        private sealed record SimulationRunProfile(
            SimulationRunMode Mode,
            int MaxRecords,
            double RtHalfWidth,
            double IntensityThresholdFraction,
            double MinIntensityThreshold);

        private static SimulationRunProfile GetSimulationRunProfile()
        {
            string? mode = Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_MODE");
            if (string.Equals(mode, "full", StringComparison.OrdinalIgnoreCase)
                || string.Equals(mode, "fullfidelity", StringComparison.OrdinalIgnoreCase))
            {
                return new SimulationRunProfile(
                    Mode: SimulationRunMode.FullFidelity,
                    MaxRecords: 200,
                    RtHalfWidth: 0.40,
                    IntensityThresholdFraction: 1e-5,
                    MinIntensityThreshold: 0.1);
            }

            return new SimulationRunProfile(
                Mode: SimulationRunMode.QuickDev,
                MaxRecords: 25,
                RtHalfWidth: 0.25,
                IntensityThresholdFraction: 1e-4,
                MinIntensityThreshold: 1.0);
        }

        private static bool GetDeduplicateProteoforms()
        {
            if (IsTrueEnvironmentVariable("MZLIB_TOPDOWN_SIM_NO_DEDUP"))
                return false;

            return true;
        }

        private static bool GetGlobalAbundanceRefitEnabled()
        {
            if (IsTrueEnvironmentVariable("MZLIB_TOPDOWN_SIM_NO_GLOBAL_ABUNDANCE_REFIT"))
                return false;

            return true;
        }

        private static int GetGlobalAbundanceRefitMaxModels()
        {
            // Was 200 while building the refit basis evaluated every model at every sample. The
            // basis now visits only models whose envelope reaches each sample, so a full run's
            // thousand-odd models fit comfortably; the cap is kept as a guard.
            const int defaultMaxModels = 10000;
            var raw = Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_GLOBAL_REFIT_MAX_MODELS");
            if (string.IsNullOrWhiteSpace(raw))
                return defaultMaxModels;

            return int.TryParse(raw, out int parsed) && parsed > 0
                ? parsed
                : defaultMaxModels;
        }

        /// <summary>
        /// Whether to fit σ_m as a function of m/z rather than collapsing every record to one
        /// scalar. On by default; a constant width makes high-m/z envelopes far easier to resolve
        /// than they are on an instrument, which is precisely the regime a feature-finder benchmark
        /// needs to be honest about.
        /// </summary>
        private static bool GetPeakWidthModelEnabled() =>
            !IsTrueEnvironmentVariable("MZLIB_TOPDOWN_SIM_CONSTANT_PEAK_WIDTH");

        /// <summary>
        /// How far the free-slope diagnostic may sit from 1.5 before the fitted width law is
        /// rejected outright rather than shipped with a warning.
        /// </summary>
        private const double MaxSlopeDeviationInStandardErrors = 3.0;

        /// <summary>
        /// Overrides how finely a window must be sampled to count as a peak shape. Set to 0 to
        /// disable the guard, which is how the centroid-scatter artifact is reproduced: without it
        /// the pooled fit happily returns σ ∝ (m/z)^0.94.
        /// </summary>
        private static double GetMinimumSamplesPerSigma()
        {
            var raw = Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_MIN_SAMPLES_PER_SIGMA");
            if (string.IsNullOrWhiteSpace(raw))
                return EnvelopeWidthFitter.DefaultMinimumSamplesPerSigma;

            return double.TryParse(raw, System.Globalization.NumberStyles.Float,
                       System.Globalization.CultureInfo.InvariantCulture, out double parsed) && parsed >= 0
                ? parsed
                : EnvelopeWidthFitter.DefaultMinimumSamplesPerSigma;
        }

        /// <summary>
        /// An explicitly supplied k for σ_m = k·(m/z)^1.5, or null to fit it from the data.
        /// </summary>
        /// <remarks>
        /// σ_m is not measurable from a centroided source — the peaks have already been reduced to
        /// positions — so on a centroided run the pooled fit is refused and the simulation falls
        /// back to a constant width, which is the defect this was meant to remove. Supplying k lets
        /// the merged-envelope regime be reached anyway. For an Orbitrap, k ≈ FWHM/(2.355·(m/z)^1.5)
        /// at any m/z where the resolving power is known: 60k at m/z 400 gives k ≈ 1.4e-7.
        /// </remarks>
        private static double? GetExplicitPeakWidthK()
        {
            var raw = Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_PEAK_WIDTH_K");
            if (string.IsNullOrWhiteSpace(raw))
                return null;

            return double.TryParse(raw, System.Globalization.NumberStyles.Float,
                       System.Globalization.CultureInfo.InvariantCulture, out double parsed) && parsed > 0
                ? parsed
                : null;
        }

        private static bool IsTrueEnvironmentVariable(string variableName)
        {
            var value = Environment.GetEnvironmentVariable(variableName);
            return string.Equals(value, "1", StringComparison.OrdinalIgnoreCase)
                   || string.Equals(value, "true", StringComparison.OrdinalIgnoreCase)
                   || string.Equals(value, "yes", StringComparison.OrdinalIgnoreCase)
                   || string.Equals(value, "on", StringComparison.OrdinalIgnoreCase);
        }


        [Test]
        [Explicit("Interactive Plotly demo for a synthetic MS1 spectrum")]
        public void PlotSyntheticSpectrum()
        {
            var simulation = BuildSimulation();
            var scan = simulation.Scans[simulation.Scans.Length / 2];

            Chart.Line<double, double, string>(scan.MassSpectrum.XArray, scan.MassSpectrum.YArray)
                .WithTraceInfo("Synthetic MS1 Spectrum")
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("m/z"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
                .WithSize(Width: 1000, Height: 500)
                .Show();
        }

        [Test]
        [Explicit("Interactive Plotly demo for a synthetic charge XIC")]
        public void PlotSyntheticChargeXic()
        {
            const double mass = 10000.0;
            var simulation = BuildSimulation();
            var extractor = new GroundTruthExtractor(simulation.Scans, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
            var truth = extractor.Extract(mass, rtCenter: 20.0, rtHalfWidth: 1.0, minCharge: 6, maxCharge: 11);

            int chargeOffset = 8 - truth.MinCharge;
            Chart.Line<double, double, string>(truth.ScanTimes, truth.ChargeXics[chargeOffset])
                .WithTraceInfo("Charge 8 XIC")
                .WithXAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Retention Time (min)"))
                .WithYAxisStyle<double, double, string>(Title: Plotly.NET.Title.init("Intensity"))
                .WithSize(Width: 1000, Height: 500)
                .Show();
        }

        private static SimulationResult BuildSimulation()
        {
            var model = new ProteoformModel(
                MonoisotopicMass: 10000.0,
                Abundance: 1.5e6,
                RtProfile: new EmgProfile(Mu: 20.0, Sigma: 0.22, Tau: 0.08),
                ChargeDistribution: new GaussianChargeDistribution(MuZ: 8.3, SigmaZ: 1.15));

            double[] scanTimes = Enumerable.Range(0, 31).Select(i => 18.5 + i * 0.1).ToArray();
            return new Simulator().Simulate(new[] { model }, minCharge: 6, maxCharge: 11, sigmaMz: 0.012, scanTimes: scanTimes);
        }

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
        public static void RawDataLoadingTimer()
        {
            var path = @"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw";
            long size = new FileInfo(path).Length;
            var sizeMb = size / (1024.0 * 1024.0);
            var sw = Stopwatch.StartNew();
            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData();
            sw.Stop();
            Console.WriteLine($"Loaded {sizeMb:F2} MB .raw file in {sw.Elapsed.ToString(@"hh\:mm\:ss\.fff")}");

            sw = Stopwatch.StartNew();
            var ms1Scans = reader.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToArray();

            sw.Stop();
            Console.WriteLine($"Collected {ms1Scans.Length} MS1 scans in {sw.Elapsed.ToString(@"hh\:mm\:ss\.fff")}");
        }

        [Test]
        public static void MzmlDataLoadingTimer()
        {
            var path = @"D:\JurkatTopdown\02-17-20_jurkat_td_rep2_fract2-calib-averaged.mzML";
            long size = new FileInfo(path).Length;
            var sizeMb = size / (1024.0 * 1024.0);
            var sw = Stopwatch.StartNew();
            var reader = MsDataFileReader.GetDataFile(path);
            reader.LoadAllStaticData();
            sw.Stop();
            Console.WriteLine($"Loaded {sizeMb:F2} MB .mzML file in {sw.Elapsed}");
        }

        [Test]
        public static void GetIsoEnv()
        {
            string histoneSeq = "MPEPAKSAPAPKKGSKKAVTKAQKKDGKKRKRSRKESYSVYVYKVLKQVHPDTGISSKAMGIMNSFVNDIFERIAGEASRLAHYNKRSTITSREIQTAVRLLLPGELAKHAVSEGTKAVTKYTSAK";

            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(histoneSeq).GetChemicalFormula();
            IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
            double[] mz = dist.Masses.Select(v => v.ToMz(20)).ToArray();
            double[] intensities = dist.Intensities.Select(v => v * 100).ToArray();
            double rt = 1;

            ChemicalFormula methyl = ChemicalFormula.Combine(new List<ChemicalFormula> { ChemicalFormula.ParseFormula("CH2"), cf } );
            ChemicalFormula acetyl = ChemicalFormula.Combine(new List<ChemicalFormula> { ChemicalFormula.ParseFormula("C2H2O"), cf });
            ChemicalFormula phospho = ChemicalFormula.Combine(new List<ChemicalFormula> { ChemicalFormula.ParseFormula("PO4H3"), cf });




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
        public static void ConvertJurkatTopdownRawToMzml()
        {
            string rawPath = @"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw";
            string outPath = Path.ChangeExtension(rawPath, ".mzML");

            var reader = MsDataFileReader.GetDataFile(rawPath);
            reader.LoadAllStaticData();

            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(reader, outPath, false);
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

        [Test]
        [Explicit("Fits TopDownSimulator models from rep2 and writes centroided simulated mzML")]
        public static void SimulateJurkatRep2AndWriteMzml()
        {
            var totalSw = Stopwatch.StartNew();
            var profile = GetSimulationRunProfile();
            bool deduplicate = GetDeduplicateProteoforms();

            Console.WriteLine($"Simulation mode: {profile.Mode} (set MZLIB_TOPDOWN_SIM_MODE=full for full-fidelity)");
            Console.WriteLine($"Deduplicate proteoforms: {deduplicate} (set MZLIB_TOPDOWN_SIM_NO_DEDUP=1 to disable)");

            string rawPath = ResolveLocalPath(@"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw");
            string mmPath = ResolveLocalPath(@"D:\JurkatTopdown\Frac7_GPTMD_Search\Task2-TopDownSearch\AllProteoforms.psmtsv");

            string stem = Path.GetFileNameWithoutExtension(rawPath);
            string outDir = Path.GetDirectoryName(rawPath)!;
            string mzmlOutPath = Path.Combine(outDir, stem + ".simulated.mzML");

            var stageSw = Stopwatch.StartNew();
            var reader = MsDataFileReader.GetDataFile(rawPath);
            reader.LoadAllStaticData();
            Console.WriteLine($"Loaded source file in {stageSw.Elapsed}");

            stageSw.Restart();
            var ms1Scans = reader.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToArray();
            Console.WriteLine($"Collected {ms1Scans.Length} MS1 scans in {stageSw.Elapsed}");

            stageSw.Restart();
            var extractor = new GroundTruthExtractor(ms1Scans, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
            Console.WriteLine($"Built ground-truth extractor in {stageSw.Elapsed}");

            stageSw.Restart();
            var loader = new MmResultLoader();
            var records = loader.Load(mmPath)
                .Where(r => string.Equals(r.FileNameWithoutExtension, stem, StringComparison.OrdinalIgnoreCase))
                .OrderByDescending(r => r.Score)
                .Take(profile.MaxRecords)
                .ToArray();

            var simulationRecords = deduplicate ? DeduplicateBySpecies(records) : AsSpecies(records);
            Console.WriteLine($"Loaded and filtered {records.Length} MM records in {stageSw.Elapsed}");
            if (simulationRecords.Length != records.Length)
                Console.WriteLine($"Deduplicated records: {simulationRecords.Length}");

            Assert.That(simulationRecords, Is.Not.Empty, "No MetaMorpheus proteoform rows matched the rep2 raw filename.");

            var fitter = new ParameterFitter(widthFitter: new EnvelopeWidthFitter(
                fallbackSigmaMz: 0.012,
                minSamplesPerSigma: GetMinimumSamplesPerSigma()));
            var fitted = new List<FittedProteoform>(simulationRecords.Length);
            var truths = new List<ProteoformGroundTruth>(simulationRecords.Length);

            int globalMinCharge = int.MaxValue;
            int globalMaxCharge = int.MinValue;

            int counter = 0;
            foreach (var (record, _, minCharge, maxCharge) in simulationRecords)
            {
                if (minCharge > maxCharge)
                    continue;

                Console.WriteLine($"Starting fit {counter}/{simulationRecords.Length}: {record.Identifier}, charge range {minCharge}-{maxCharge}");

                var truth = extractor.Extract(record.MonoisotopicMass, record.RetentionTime, rtHalfWidth: profile.RtHalfWidth, minCharge: minCharge, maxCharge: maxCharge);
                FittedProteoform fit;
                try
                {
                    fit = fitter.Fit(truth, record.Identifier);
                }
                catch (InvalidOperationException)
                {
                    continue;
                }

                if (double.IsNaN(fit.Model.Abundance) || fit.Model.Abundance <= 0)
                    continue;

                fitted.Add(fit);
                truths.Add(truth);
                globalMinCharge = Math.Min(globalMinCharge, minCharge);
                globalMaxCharge = Math.Max(globalMaxCharge, maxCharge);
                counter++;
                Console.WriteLine($"Completed fit {counter}/{simulationRecords.Length}: {record.Identifier}, charge range {minCharge}-{maxCharge}, fitted abundance {fit.Model.Abundance:F2}, sigmaMz {fit.SigmaMz:F6}");
            }

            Assert.That(fitted, Is.Not.Empty, "No fit records were produced from rep2 MM results.");

            var sigmaCandidates = fitted
                .Select(f => f.SigmaMz)
                .Where(s => !double.IsNaN(s) && !double.IsInfinity(s) && s > 0)
                .OrderBy(s => s)
                .ToArray();
            double sigmaMz = sigmaCandidates.Length == 0 ? 0.012 : sigmaCandidates[sigmaCandidates.Length / 2];

            if (globalMinCharge == int.MaxValue)
                globalMinCharge = 2;
            if (globalMaxCharge == int.MinValue)
                globalMaxCharge = 80;

            var (widthModel, fittedUnderSharedWidth) = FitSharedPeakWidth(
                fitted.ToArray(), truths.ToArray(), sigmaMz, GetPeakWidthModelEnabled());

            var scanTimes = ms1Scans.Select(s => s.RetentionTime).ToArray();
            var models = fittedUnderSharedWidth.Select(f => f.Model).ToArray();

            stageSw.Restart();
            var export = new Simulator().WriteMzml(
                models,
                globalMinCharge,
                globalMaxCharge,
                widthModel,
                scanTimes,
                mzmlOutPath,
                new ScanReductionOptions
                {
                    RelativeIntensityThreshold = profile.IntensityThresholdFraction,
                    MinimumIntensity = profile.MinIntensityThreshold,
                });
            Console.WriteLine($"Simulated and wrote mzML in {stageSw.Elapsed}");

            Console.WriteLine($"Simulated models: {models.Length}");
            Console.WriteLine($"Median per-record sigmaMz: {sigmaMz:F6}");
            Console.WriteLine($"Charge range: {globalMinCharge}-{globalMaxCharge}");
            Console.WriteLine($"Wrote simulated mzML: {export.MzmlPath}");
            Console.WriteLine($"Parameter truth: {export.GroundTruthPath}");
            Console.WriteLine($"Feature truth: {export.FeatureGroundTruthPath} ({export.FeatureCount} features)");
            Console.WriteLine($"Simulated scans: {export.ScanCount}, peak count: {export.PeakCount}");
            Console.WriteLine($"Total elapsed: {totalSw.Elapsed}");
        }

        [Test]
        [Explicit("Writes centroided simulated mzML for the 31-35 min slice and the full q<=0.01 run of rep2 fract7")]
        public static void ExportRep2SliceAndQValueSimulations()
        {
            const double rtStart = 31.0;
            const double rtEnd = 35.0;
            const double qValueThreshold = 0.01;
            const double rtHalfWidth = 0.25;
            bool deduplicate = GetDeduplicateProteoforms();
            bool useGlobalAbundanceRefit = GetGlobalAbundanceRefitEnabled();
            int globalRefitMaxModels = GetGlobalAbundanceRefitMaxModels();

            string rawPath = ResolveLocalPath(@"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw");
            string resultPath = ResolveLocalPath(@"D:\JurkatTopdown\Frac7_GPTMD_Search\Task2-TopDownSearch\Individual File Results\02-18-20_jurkat_td_rep2_fract7_Proteoforms.psmtsv");
            string stem = Path.GetFileNameWithoutExtension(rawPath);
            string outDir = Path.GetDirectoryName(rawPath)!;

            string simSliceMzmlPath = Path.Combine(outDir, stem + ".rt31-35.simulated.q001.mzML");
            string simFullMzmlPath = Path.Combine(outDir, stem + ".full.simulated.q001.mzML");

            var reduction = new ScanReductionOptions
            {
                RelativeIntensityThreshold = 1e-4,
                MinimumIntensity = 1.0,
            };

            var totalSw = Stopwatch.StartNew();
            Console.WriteLine($"Global abundance refit: {useGlobalAbundanceRefit} (set MZLIB_TOPDOWN_SIM_NO_GLOBAL_ABUNDANCE_REFIT=1 to disable)");
            Console.WriteLine($"Global abundance refit max models: {globalRefitMaxModels} (override with MZLIB_TOPDOWN_SIM_GLOBAL_REFIT_MAX_MODELS)");

            var reader = MsDataFileReader.GetDataFile(rawPath);
            reader.LoadAllStaticData();
            var allMs1Scans = reader.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToArray();
            Assert.That(allMs1Scans, Is.Not.Empty, "No MS1 scans were found in the raw file.");

            var rtSliceMs1Scans = allMs1Scans
                .Where(s => s.RetentionTime >= rtStart && s.RetentionTime <= rtEnd)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToArray();
            Assert.That(rtSliceMs1Scans, Is.Not.Empty, "No MS1 scans were found in the 31-35 min window.");

            // Real data is not re-exported here; ConvertJurkatTopdownRawToMzml handles raw -> mzML.
            // Its centroid status is logged because mzLib cannot read a profile-mode mzML back.
            Console.WriteLine($"Real slice scans: {rtSliceMs1Scans.Length} " +
                              $"(centroided: {rtSliceMs1Scans.Count(s => s.IsCentroid)}/{rtSliceMs1Scans.Length})");

            var loadedRecords = LoadQualifiedMmRecords(resultPath, stem, qValueThreshold, rtStart: null, rtEnd: null);
            var qFilteredRecords = deduplicate ? DeduplicateBySpecies(loadedRecords) : AsSpecies(loadedRecords);
            if (deduplicate)
                Console.WriteLine($"Deduplicated q<=0.01 records: {loadedRecords.Length} -> {qFilteredRecords.Length} species");
            var qFilteredSliceRecords = qFilteredRecords
                .Where(s => s.Anchor.RetentionTime >= rtStart && s.Anchor.RetentionTime <= rtEnd)
                .ToArray();

            Assert.That(qFilteredSliceRecords, Is.Not.Empty, "No proteoforms passed q<=0.01 inside 31-35 min.");
            Assert.That(qFilteredRecords, Is.Not.Empty, "No proteoforms passed q<=0.01 for full run.");

            var extractor = new GroundTruthExtractor(allMs1Scans, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);

            // The slice records are a subset of the full set, so everything is extracted and fitted
            // once and the slice is then selected out of it. The global refit is a joint fit, so it
            // is still applied separately to each output set — the slice refit is conditioned on the
            // slice, exactly as before.
            var stageSw = Stopwatch.StartNew();
            var allFits = FitProteoforms(qFilteredRecords, extractor, rtHalfWidth, GetPeakWidthModelEnabled());
            Console.WriteLine($"Fitted {allFits.Fits.Length}/{qFilteredRecords.Length} records in {stageSw.Elapsed}");
            Assert.That(allFits.Fits, Is.Not.Empty, "No simulated models were fitted for full q<=0.01 run.");

            var sliceSelection = Enumerable.Range(0, allFits.Fits.Length)
                .Where(i => allFits.Records[i].RetentionTime >= rtStart && allFits.Records[i].RetentionTime <= rtEnd)
                .ToArray();
            Assert.That(sliceSelection, Is.Not.Empty, "No fitted proteoforms fell inside 31-35 min.");

            var sliceModels = ApplyGlobalRefit(
                sliceSelection.Select(i => allFits.Fits[i]).ToArray(),
                sliceSelection.Select(i => allFits.Truths[i]).ToArray(),
                allFits.MinCharge, allFits.MaxCharge, allFits.WidthModel,
                useGlobalAbundanceRefit, globalRefitMaxModels, "31-35 min slice");

            var sliceScanTimes = rtSliceMs1Scans.Select(s => s.RetentionTime).ToArray();
            var sliceExport = new Simulator().WriteMzml(
                sliceModels,
                allFits.MinCharge,
                allFits.MaxCharge,
                allFits.WidthModel,
                sliceScanTimes,
                simSliceMzmlPath,
                reduction);

            Console.WriteLine($"Simulated slice mzML written: {sliceExport.MzmlPath}");
            Console.WriteLine($"Parameter truth: {sliceExport.GroundTruthPath}");
            Console.WriteLine($"Feature truth: {sliceExport.FeatureGroundTruthPath} ({sliceExport.FeatureCount} features)");
            Console.WriteLine($"Simulated slice models: {sliceModels.Length}, scans: {sliceExport.ScanCount}, peaks: {sliceExport.PeakCount}");

            var fullModels = ApplyGlobalRefit(
                allFits.Fits, allFits.Truths,
                allFits.MinCharge, allFits.MaxCharge, allFits.WidthModel,
                useGlobalAbundanceRefit, globalRefitMaxModels, "full q<=0.01 run");

            var fullScanTimes = allMs1Scans.Select(s => s.RetentionTime).ToArray();
            var fullExport = new Simulator().WriteMzml(
                fullModels,
                allFits.MinCharge,
                allFits.MaxCharge,
                allFits.WidthModel,
                fullScanTimes,
                simFullMzmlPath,
                reduction);

            Console.WriteLine($"Simulated full mzML written: {fullExport.MzmlPath}");
            Console.WriteLine($"Parameter truth: {fullExport.GroundTruthPath}");
            Console.WriteLine($"Feature truth: {fullExport.FeatureGroundTruthPath} ({fullExport.FeatureCount} features)");
            Console.WriteLine($"Simulated full models: {fullModels.Length}, scans: {fullExport.ScanCount}, peaks: {fullExport.PeakCount}");
            Console.WriteLine($"Total elapsed: {totalSw.Elapsed}");
        }

        /// <summary>
        /// Writes a matched pair of simulated mzMLs — one clean, one with an injected noise floor —
        /// over an RT window that spans both the tail of elution and the post-elution region where
        /// the real file's noise is unobscured, so the two can be opened next to the real data.
        /// </summary>
        /// <remarks>
        /// The noise amplitude is the one measured in this file's own 50-55 min window
        /// (<c>NoiseCharacterization</c>), so the simulated floor should sit at the same absolute
        /// intensity as the real one rather than merely having the right shape.
        /// </remarks>
        [Test]
        [Explicit("Writes clean and noisy simulated mzML for rep2 fract7 over the post-elution window")]
        public static void ExportRep2NoisySimulationForInspection() =>
            ExportRep2NoisySimulation(rtStart: 45.0, rtEnd: 56.0, label: "rt45-56", writeClean: true);

        /// <summary>
        /// The whole run rather than a window, which is what a feature-finding or deconvolution
        /// benchmark actually needs — a slice cannot exercise anything that depends on elution
        /// order, on the full dynamic range, or on the post-elution region being reached from a busy
        /// one.
        /// </summary>
        /// <remarks>
        /// Expect roughly 2900 MS1 scans, ~50 M peaks and an mzML around a gigabyte. A clean
        /// full-run simulation already exists as <c>*.full.simulated.q001.mzML</c>, so only the noisy
        /// file is written here.
        /// </remarks>
        [Test]
        [Explicit("Writes a full-run noisy simulated mzML for rep2 fract7 (~1 GB, several minutes)")]
        public static void ExportRep2FullNoisySimulation() =>
            ExportRep2NoisySimulation(rtStart: null, rtEnd: null, label: "full", writeClean: false);

        /// <summary>
        /// The same full-run export for rep1 fract7, fitted to its own raw file. It is the upper bound
        /// an ID-only simulation of rep1 is compared against; see <c>IdOnlySimulation</c>.
        /// </summary>
        [Test]
        [Explicit("Writes a full-run noisy simulated mzML for rep1 fract7, fitted to its own raw file")]
        public static void ExportRep1FullNoisySimulation() =>
            ExportNoisySimulation(
                @"D:\JurkatTopdown\02-18-20_jurkat_td_rep1_fract7.raw",
                @"D:\JurkatTopdown\Frac7_GPTMD_Search\Task2-TopDownSearch\Individual File Results\02-18-20_jurkat_td_rep1_fract7_Proteoforms.psmtsv",
                rtStart: null, rtEnd: null, label: "full", writeClean: false);

        private static void ExportRep2NoisySimulation(double? rtStart, double? rtEnd, string label, bool writeClean) =>
            ExportNoisySimulation(
                @"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw",
                @"D:\JurkatTopdown\Frac7_GPTMD_Search\Task2-TopDownSearch\Individual File Results\02-18-20_jurkat_td_rep2_fract7_Proteoforms.psmtsv",
                rtStart, rtEnd, label, writeClean);

        private static void ExportNoisySimulation(
            string rawFile, string resultFile, double? rtStart, double? rtEnd, string label, bool writeClean)
        {
            const double qValueThreshold = 0.01;
            const double rtHalfWidth = 0.25;

            // Median noise amplitude at m/z 650 measured in rep2 fract7's 50-55 min window. Only the
            // unconditioned path uses it; per-scan conditioning takes the amplitude from each scan's
            // injection time.
            const double measuredNoiseLevel = 663.0;

            double densityScale = GetNoiseDensityScale();

            string rawPath = ResolveLocalPath(rawFile);
            string resultPath = ResolveLocalPath(resultFile);
            string stem = Path.GetFileNameWithoutExtension(rawPath);
            string outDir = Path.GetDirectoryName(rawPath)!;

            label += GetOutputTag();
            string cleanPath = Path.Combine(outDir, $"{stem}.{label}.clean.simulated.mzML");
            string noisyPath = Path.Combine(outDir, $"{stem}.{label}.noisy.simulated.mzML");

            var totalSw = Stopwatch.StartNew();

            var reader = MsDataFileReader.GetDataFile(rawPath);
            reader.LoadAllStaticData();
            var allMs1Scans = reader.GetAllScansList()
                .Where(s => s.MsnOrder == 1)
                .OrderBy(s => s.OneBasedScanNumber)
                .ToArray();

            var windowScans = allMs1Scans
                .Where(s => !rtStart.HasValue || s.RetentionTime >= rtStart.Value)
                .Where(s => !rtEnd.HasValue || s.RetentionTime <= rtEnd.Value)
                .ToArray();
            Assert.That(windowScans, Is.Not.Empty, $"No MS1 scans between {rtStart} and {rtEnd} min.");

            double realPeaksPerScan = windowScans.Average(s => s.MassSpectrum.XArray.Length);
            Console.WriteLine($"Real scans ({label}, RT {windowScans[0].RetentionTime:F2}-" +
                              $"{windowScans[^1].RetentionTime:F2} min): {windowScans.Length}, " +
                              $"mean {realPeaksPerScan:F0} peaks/scan");

            // Fitted against the whole run so that proteoforms eluting just before the window still
            // contribute their tails, then simulated only on the window's scan grid.
            var records = DeduplicateBySpecies(
                LoadQualifiedMmRecords(resultPath, stem, qValueThreshold, rtStart: null, rtEnd: null));
            Console.WriteLine($"Deduplicated q<={qValueThreshold} records: {records.Length}");

            var extractor = new GroundTruthExtractor(allMs1Scans, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
            var stageSw = Stopwatch.StartNew();
            var allFits = FitProteoforms(records, extractor, rtHalfWidth, GetPeakWidthModelEnabled());
            Console.WriteLine($"Fitted {allFits.Fits.Length}/{records.Length} records in {stageSw.Elapsed}");
            Assert.That(allFits.Fits, Is.Not.Empty, "No proteoforms were fitted.");

            var models = ApplyGlobalRefit(
                allFits.Fits, allFits.Truths, allFits.MinCharge, allFits.MaxCharge, allFits.WidthModel,
                GetGlobalAbundanceRefitEnabled(), GetGlobalAbundanceRefitMaxModels(), label);

            double[] scanTimes = windowScans.Select(s => s.RetentionTime).ToArray();
            var simulator = new Simulator();

            if (writeClean)
            {
                stageSw.Restart();
                var clean = simulator.WriteMzml(
                    models, allFits.MinCharge, allFits.MaxCharge, allFits.WidthModel, scanTimes, cleanPath,
                    new ScanReductionOptions { RelativeIntensityThreshold = 1e-4, MinimumIntensity = 1.0 });
                Console.WriteLine($"Clean mzML: {cleanPath}");
                Console.WriteLine($"  scans {clean.ScanCount}, peaks {clean.PeakCount} " +
                                  $"({clean.PeakCount / (double)clean.ScanCount:F0}/scan), " +
                                  $"features {clean.FeatureCount}, {stageSw.Elapsed}");
            }

            var noise = new NoiseFloorModel(
                noiseLevelAtReferenceMz: measuredNoiseLevel,
                densityScale: densityScale);
            var conditioning = GetNoiseConditioning();
            var scanNoise = conditioning == NoiseConditioning.None
                ? null
                : ScanNoiseConditions.FromSourceScans(
                    windowScans, noise, conditionDensity: conditioning == NoiseConditioning.Full);
            Console.WriteLine(scanNoise is null
                ? $"Noise model: level {measuredNoiseLevel} at m/z {NoiseFloorModel.ReferenceMz}, " +
                  $"density scale {densityScale}, {noise.ExpectedPeaksPerScan:F0} peaks/scan expected"
                : $"Noise model: conditioned per scan ({conditioning}) on the source scans, density scale {densityScale}, " +
                  $"{scanNoise.Average(m => m.ExpectedPeaksPerScan):F0} peaks/scan expected on average");

            stageSw.Restart();
            var noisy = simulator.WriteMzml(
                models, allFits.MinCharge, allFits.MaxCharge, allFits.WidthModel, scanTimes, noisyPath,
                noise: noise, scanNoise: scanNoise);
            Console.WriteLine($"Noisy mzML: {noisyPath}");
            Console.WriteLine($"  scans {noisy.ScanCount}, peaks {noisy.PeakCount} " +
                              $"({noisy.PeakCount / (double)noisy.ScanCount:F0}/scan), " +
                              $"features {noisy.FeatureCount}, {stageSw.Elapsed}");
            Console.WriteLine($"  signal {noisy.Noise!.SignalPeaks}, noise {noisy.Noise.NoisePeaks}, " +
                              $"merged {noisy.Noise.MergedPeaks}, " +
                              $"dropped by jitter {noisy.Noise.DroppedByJitter} " +
                              $"({noisy.Noise.DroppedByJitter / (double)(noisy.Noise.SignalPeaks + noisy.Noise.DroppedByJitter):P1} " +
                              "of signal)");
            Console.WriteLine($"  real {realPeaksPerScan:F0} peaks/scan vs simulated " +
                              $"{noisy.PeakCount / (double)noisy.ScanCount:F0}");

            foreach (string path in new[] { cleanPath, noisyPath }.Where(File.Exists))
                Console.WriteLine($"{Path.GetFileName(path)}: {new FileInfo(path).Length / 1024.0 / 1024.0:F1} MB");

            Console.WriteLine($"Peak working set: {Environment.WorkingSet / 1024.0 / 1024.0 / 1024.0:F1} GB");
            Console.WriteLine($"Total elapsed: {totalSw.Elapsed}");
        }

        /// <summary>
        /// Writes a run in which every identified proteoform's monoisotopic mass is offset by a
        /// constant X daltons and every MS2 isolation window is offset by X/z to follow it.
        /// </summary>
        /// <remarks>
        /// MS1 is simulated at the shifted masses; the MS2 scans are the real ones, carried across
        /// with their peak lists untouched. The fragments therefore still describe the unshifted
        /// proteoform, so the file reads as an unlocalized modification of mass X on every
        /// identification — which is what makes it usable as a decoy or open-search benchmark.
        /// </remarks>
        [Test]
        [Explicit("Writes a precursor-mass-shifted simulated mzML for rep2 fract7 (MS1 simulated, MS2 carried over)")]
        public static void ExportRep2MassShiftedSimulation() =>
            ExportRep2ShiftedSimulation(injectNoise: false);

        [Test]
        [Explicit("Writes a precursor-mass-shifted simulated mzML for rep2 fract7 with an injected noise floor")]
        public static void ExportRep2MassShiftedNoisySimulation() =>
            ExportRep2ShiftedSimulation(injectNoise: true);

        private static void ExportRep2ShiftedSimulation(bool injectNoise)
        {
            const double qValueThreshold = 0.01;
            const double rtHalfWidth = 0.25;

            // Median noise amplitude at m/z 650 measured in this file's 50-55 min window.
            const double measuredNoiseLevel = 663.0;

            double massShiftDa = GetMassShiftDaltons();

            string rawPath = ResolveLocalPath(@"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw");
            string resultPath = ResolveLocalPath(@"D:\JurkatTopdown\Frac7_GPTMD_Search\Task2-TopDownSearch\Individual File Results\02-18-20_jurkat_td_rep2_fract7_Proteoforms.psmtsv");
            string stem = Path.GetFileNameWithoutExtension(rawPath);
            string outDir = Path.GetDirectoryName(rawPath)!;

            string shiftLabel = massShiftDa.ToString("+0.###;-0.###", System.Globalization.CultureInfo.InvariantCulture);
            string outPath = Path.Combine(outDir,
                $"{stem}.shift{shiftLabel}Da{(injectNoise ? ".noisy" : "")}.simulated.mzML");

            var totalSw = Stopwatch.StartNew();
            Console.WriteLine($"Precursor mass shift: {massShiftDa} Da " +
                              "(override with MZLIB_TOPDOWN_SIM_MASS_SHIFT_DA)");

            var reader = MsDataFileReader.GetDataFile(rawPath);
            reader.LoadAllStaticData();
            var allScans = reader.GetAllScansList().OrderBy(s => s.OneBasedScanNumber).ToArray();

            var allMs1Scans = allScans.Where(s => s.MsnOrder == 1).ToArray();
            var allMs2Scans = allScans.Where(s => s.MsnOrder > 1).ToArray();
            Assert.That(allMs1Scans, Is.Not.Empty, "No MS1 scans were found in the raw file.");
            Assert.That(allMs2Scans, Is.Not.Empty, "No MS2 scans were found in the raw file.");
            Console.WriteLine($"Source scans: {allMs1Scans.Length} MS1, {allMs2Scans.Length} MS2");

            var records = DeduplicateBySpecies(
                LoadQualifiedMmRecords(resultPath, stem, qValueThreshold, rtStart: null, rtEnd: null));
            Console.WriteLine($"Deduplicated q<={qValueThreshold} records: {records.Length}");

            var extractor = new GroundTruthExtractor(allMs1Scans, ppmTolerance: 20.0, mzWindowHalfWidth: 0.05);
            var stageSw = Stopwatch.StartNew();
            var allFits = FitProteoforms(records, extractor, rtHalfWidth, GetPeakWidthModelEnabled());
            Console.WriteLine($"Fitted {allFits.Fits.Length}/{records.Length} records in {stageSw.Elapsed}");
            Assert.That(allFits.Fits, Is.Not.Empty, "No proteoforms were fitted.");

            var models = ApplyGlobalRefit(
                allFits.Fits, allFits.Truths, allFits.MinCharge, allFits.MaxCharge, allFits.WidthModel,
                GetGlobalAbundanceRefitEnabled(), GetGlobalAbundanceRefitMaxModels(), "mass-shifted run");

            var noise = injectNoise
                ? new NoiseFloorModel(
                    noiseLevelAtReferenceMz: measuredNoiseLevel,
                    densityScale: GetNoiseDensityScale())
                : null;

            double[] scanTimes = allMs1Scans.Select(s => s.RetentionTime).ToArray();

            stageSw.Restart();
            var export = new Simulator().WriteShiftedMzml(
                models,
                allFits.MinCharge,
                allFits.MaxCharge,
                allFits.WidthModel,
                scanTimes,
                allMs2Scans,
                massShiftDa,
                outPath,
                noise is null
                    ? new ScanReductionOptions { RelativeIntensityThreshold = 1e-4, MinimumIntensity = 1.0 }
                    : null,
                noise: noise);
            Console.WriteLine($"Simulated and wrote shifted mzML in {stageSw.Elapsed}");

            var shift = export.Shift!;
            Console.WriteLine($"Shifted mzML: {export.MzmlPath}");
            Console.WriteLine($"  scans {export.ScanCount} ({export.ScanCount - shift.Ms2ScansWritten} MS1, " +
                              $"{shift.Ms2ScansWritten} MS2), peaks {export.PeakCount}");
            Console.WriteLine($"  shifted proteoforms: {models.Length}, features {export.FeatureCount}");
            Console.WriteLine($"  MS2 dropped: {shift.Ms2ScansDroppedWithoutChargeState} with no charge state, " +
                              $"{shift.Ms2ScansDroppedWithoutPrecursorMz} with no precursor m/z, " +
                              $"{shift.Ms2ScansDroppedProfileMode} profile mode, " +
                              $"{shift.Ms2ScansDroppedOutsideMs1TimeRange} outside the MS1 time range, " +
                              $"{shift.Ms2ScansDroppedWithEmptyPrecursor} with an empty precursor scan");
            Console.WriteLine($"  parameter truth: {export.GroundTruthPath}");
            Console.WriteLine($"  feature truth:   {export.FeatureGroundTruthPath}");
            Console.WriteLine($"  mass shift map:  {shift.ShiftSidecarPath}");

            if (File.Exists(outPath))
                Console.WriteLine($"{Path.GetFileName(outPath)}: {new FileInfo(outPath).Length / 1024.0 / 1024.0:F1} MB");

            Console.WriteLine($"Total elapsed: {totalSw.Elapsed}");
        }

        /// <summary>
        /// Checks a written mass-shifted run against the source raw file and its own sidecars.
        /// </summary>
        /// <remarks>
        /// The unit tests prove the shift on synthetic scans; this proves it on the file that was
        /// actually produced, which is where a scan-numbering or precursor-linking mistake would
        /// show up. Reading the mzML back at all is itself a check: mzLib refuses a file containing
        /// any profile spectrum.
        /// </remarks>
        [Test]
        [Explicit("Verifies a previously written rep2 mass-shifted mzML against the source raw file")]
        public static void VerifyRep2MassShiftedSimulation()
        {
            double massShiftDa = GetMassShiftDaltons();

            string rawPath = ResolveLocalPath(@"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw");
            string stem = Path.GetFileNameWithoutExtension(rawPath);
            string outDir = Path.GetDirectoryName(rawPath)!;
            string shiftLabel = massShiftDa.ToString("+0.###;-0.###", System.Globalization.CultureInfo.InvariantCulture);
            string shiftedPath = Path.Combine(outDir, $"{stem}.shift{shiftLabel}Da.simulated.mzML");

            Assert.That(File.Exists(shiftedPath), Is.True,
                $"{shiftedPath} does not exist; run ExportRep2MassShiftedSimulation first.");

            var shiftedReader = MsDataFileReader.GetDataFile(shiftedPath);
            shiftedReader.LoadAllStaticData();
            var shiftedScans = shiftedReader.GetAllScansList().ToArray();

            var rawReader = MsDataFileReader.GetDataFile(rawPath);
            rawReader.LoadAllStaticData();
            var rawMs2 = rawReader.GetAllScansList()
                .Where(s => s.MsnOrder > 1)
                .OrderBy(s => s.RetentionTime)
                .ToArray();

            var shiftedMs1 = shiftedScans.Where(s => s.MsnOrder == 1).ToArray();
            var shiftedMs2 = shiftedScans.Where(s => s.MsnOrder > 1).OrderBy(s => s.RetentionTime).ToArray();
            Console.WriteLine($"Shifted file: {shiftedScans.Length} scans ({shiftedMs1.Length} MS1, {shiftedMs2.Length} MS2)");

            // Numbering has to equal position, or every precursor reference in the file is wrong.
            for (int i = 0; i < shiftedScans.Length; i++)
                Assert.That(shiftedScans[i].OneBasedScanNumber, Is.EqualTo(i + 1));

            Assert.That(shiftedScans.Select(s => s.RetentionTime), Is.Ordered.Ascending);

            foreach (var scan in shiftedMs2)
            {
                Assert.That(scan.OneBasedPrecursorScanNumber, Is.Not.Null,
                    $"MS2 scan {scan.OneBasedScanNumber} has no precursor reference.");

                var precursor = shiftedScans[scan.OneBasedPrecursorScanNumber.Value - 1];
                Assert.That(precursor.MsnOrder, Is.EqualTo(1),
                    $"MS2 scan {scan.OneBasedScanNumber} points at an MS{precursor.MsnOrder}.");
                Assert.That(precursor.RetentionTime, Is.LessThanOrEqualTo(scan.RetentionTime));
            }

            // The written MS2 scans are a retention-time-ordered subsequence of the source ones:
            // scans sitting behind an empty survey scan are deliberately left out, so this is a
            // two-pointer match rather than a 1:1 zip.
            Assert.That(shiftedMs2.Length, Is.LessThanOrEqualTo(rawMs2.Length));
            Console.WriteLine($"MS2 carried over: {shiftedMs2.Length}/{rawMs2.Length} " +
                              $"({rawMs2.Length - shiftedMs2.Length} dropped)");

            double worstError = 0;
            int r = 0;
            foreach (var scan in shiftedMs2)
            {
                while (r < rawMs2.Length && Math.Abs(rawMs2[r].RetentionTime - scan.RetentionTime) > 1e-9)
                    r++;

                Assert.That(r, Is.LessThan(rawMs2.Length),
                    $"No source MS2 scan at retention time {scan.RetentionTime}.");
                var source = rawMs2[r];
                r++;

                double expectedShift = massShiftDa / Math.Abs(source.SelectedIonChargeStateGuess!.Value);
                double observedShift = scan.IsolationMz!.Value - source.IsolationMz!.Value;

                Assert.That(observedShift, Is.EqualTo(expectedShift).Within(1e-4),
                    $"Scan {scan.OneBasedScanNumber} isolation window moved by {observedShift}, expected {expectedShift}.");
                Assert.That(scan.IsolationWidth!.Value,
                    Is.EqualTo(source.IsolationWidth!.Value).Within(1e-6), "The window widened.");
                Assert.That(scan.SelectedIonChargeStateGuess, Is.EqualTo(source.SelectedIonChargeStateGuess));

                // Fragments must not have moved: this is what makes the shift read as an
                // unlocalized modification rather than a uniformly heavier molecule. Comparing the
                // positions, not just the count, also confirms the two-pointer match is aligned.
                Assert.That(scan.MassSpectrum.XArray,
                    Is.EqualTo(source.MassSpectrum.XArray).Within(1e-6).AsCollection);

                worstError = Math.Max(worstError, Math.Abs(observedShift - expectedShift));
            }

            Console.WriteLine($"Isolation windows: all {shiftedMs2.Length} moved by shift/z, worst error {worstError:E2} m/z");

            // A precursor scan with no peaks is fatal downstream, not merely unhelpful:
            // MsDataScan.RefineSelectedMzAndIntensity throws on an empty precursor spectrum, and
            // MetaMorpheus's GetMs2Scans catches that with a `continue` that skips the assignment of
            // its scansWithPrecursors slot — the resulting null then throws ArgumentNullException out
            // of a SelectMany and takes down the whole search.
            int emptyMs1 = shiftedMs1.Count(s => s.MassSpectrum.XArray.Length == 0);
            var orphanedMs2 = shiftedMs2
                .Where(s => shiftedScans[s.OneBasedPrecursorScanNumber!.Value - 1].MassSpectrum.XArray.Length == 0)
                .ToArray();

            Console.WriteLine($"Empty MS1 scans: {emptyMs1}/{shiftedMs1.Length}; " +
                              $"MS2 scans whose precursor scan is empty: {orphanedMs2.Length}");

            Assert.That(orphanedMs2, Is.Empty,
                $"{orphanedMs2.Length} MS2 scans reference an MS1 scan with no peaks, which makes " +
                "MetaMorpheus throw ArgumentNullException out of GetMs2Scans for the whole file.");

            // The feature truth must describe the shifted file, and every scan it quotes must be MS1.
            string featurePath = Path.ChangeExtension(shiftedPath, ".features.tsv");
            var featureLines = File.ReadAllLines(featurePath);
            var header = featureLines[0].Split('\t');
            int massColumn = Array.IndexOf(header, "MonoisotopicMass");
            int apexScanColumn = Array.IndexOf(header, "ApexScanNumber");

            string shiftMapPath = Path.ChangeExtension(shiftedPath, ".massshift.tsv");
            var shiftMapLines = File.ReadAllLines(shiftMapPath);
            var shiftMapHeader = shiftMapLines[0].Split('\t');
            int shiftedMassColumn = Array.IndexOf(shiftMapHeader, "ShiftedMonoisotopicMass");
            int originalMassColumn = Array.IndexOf(shiftMapHeader, "OriginalMonoisotopicMass");

            var shiftedMasses = shiftMapLines.Skip(1)
                .Select(l => l.Split('\t'))
                .Select(r => double.Parse(r[shiftedMassColumn], System.Globalization.CultureInfo.InvariantCulture))
                .ToHashSet();

            foreach (var row in shiftMapLines.Skip(1).Select(l => l.Split('\t')))
            {
                double original = double.Parse(row[originalMassColumn], System.Globalization.CultureInfo.InvariantCulture);
                double shifted = double.Parse(row[shiftedMassColumn], System.Globalization.CultureInfo.InvariantCulture);
                Assert.That(shifted - original, Is.EqualTo(massShiftDa).Within(1e-6));
            }

            foreach (var row in featureLines.Skip(1).Select(l => l.Split('\t')))
            {
                double mass = double.Parse(row[massColumn], System.Globalization.CultureInfo.InvariantCulture);
                Assert.That(shiftedMasses.Any(m => Math.Abs(m - mass) < 1e-6), Is.True,
                    $"Feature mass {mass} is not one of the shifted proteoform masses.");

                int apexScan = int.Parse(row[apexScanColumn]);
                Assert.That(shiftedScans[apexScan - 1].MsnOrder, Is.EqualTo(1),
                    $"Feature apex scan {apexScan} is an MS{shiftedScans[apexScan - 1].MsnOrder}.");
            }

            Console.WriteLine($"Feature truth: {featureLines.Length - 1} features, all at shifted masses, " +
                              "all apex scans are MS1");
            Console.WriteLine($"Mass shift map: {shiftMapLines.Length - 1} proteoforms, all offset by {massShiftDa} Da");
        }

        /// <summary>
        /// The neutral offset applied to every identified proteoform. The default is deliberately
        /// not a common modification mass and not a multiple of the 1.00335 Da isotope spacing, so a
        /// shifted feature cannot be mistaken for a real mod or for a misassigned monoisotopic peak.
        /// </summary>
        private static double GetMassShiftDaltons()
        {
            const double defaultShift = 10.0;
            var raw = Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_MASS_SHIFT_DA");
            if (string.IsNullOrWhiteSpace(raw))
                return defaultShift;

            return double.TryParse(raw, System.Globalization.NumberStyles.Float,
                       System.Globalization.CultureInfo.InvariantCulture, out double parsed)
                   && double.IsFinite(parsed) && parsed != 0
                ? parsed
                : defaultShift;
        }

        /// <summary>
        /// Fraction of the measured noise density to emit. The full density is ~16 900 peaks per
        /// scan, which is what the instrument reports but which makes for a large file; lower it to
        /// iterate.
        /// </summary>
        private static double GetNoiseDensityScale()
        {
            var raw = Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_NOISE_DENSITY_SCALE");
            if (string.IsNullOrWhiteSpace(raw))
                return 1.0;

            return double.TryParse(raw, out double parsed) && parsed >= 0 ? parsed : 1.0;
        }

        private enum NoiseConditioning { None, Amplitude, Full }

        /// <summary>
        /// How the injected noise follows the source run. <c>full</c> (the default) takes each
        /// scan's amplitude from its injection time and its density from its own low-S/N peaks;
        /// <c>amplitude</c> conditions only the amplitude; <c>none</c> uses one model for every
        /// scan, as before. Set with MZLIB_TOPDOWN_SIM_NOISE_CONDITIONING.
        /// </summary>
        private static NoiseConditioning GetNoiseConditioning()
        {
            var raw = Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_NOISE_CONDITIONING");
            return raw?.Trim().ToLowerInvariant() switch
            {
                "none" => NoiseConditioning.None,
                "amplitude" => NoiseConditioning.Amplitude,
                _ => NoiseConditioning.Full,
            };
        }

        /// <summary>
        /// Appended to the output label, e.g. ".v2", so a new export does not overwrite one that
        /// downstream results were computed from. Set with MZLIB_TOPDOWN_SIM_OUTPUT_TAG.
        /// </summary>
        private static string GetOutputTag() =>
            Environment.GetEnvironmentVariable("MZLIB_TOPDOWN_SIM_OUTPUT_TAG")?.Trim() ?? string.Empty;

        private static MmResultRecord[] LoadQualifiedMmRecords(
            string psmTsvPath,
            string expectedFileStem,
            double qValueThreshold,
            double? rtStart,
            double? rtEnd)
        {
            var file = new PsmFromTsvFile(psmTsvPath, new SpectrumMatchParsingParameters
            {
                ParseMatchedFragmentIons = false,
            });
            file.LoadResults();

            var records = file.Results
                .Where(p => p is not null)
                .Where(p => string.Equals(p.FileNameWithoutExtension, expectedFileStem, StringComparison.OrdinalIgnoreCase))
                .Where(p => p.MonoisotopicMass > 0 && p.RetentionTime >= 0)
                .Where(p => !double.IsNaN(p.QValue) && p.QValue <= qValueThreshold)
                .Where(p => !rtStart.HasValue || p.RetentionTime >= rtStart.Value)
                .Where(p => !rtEnd.HasValue || p.RetentionTime <= rtEnd.Value)
                .Select(p => new MmResultRecord(
                    FileNameWithoutExtension: p.FileNameWithoutExtension,
                    PrecursorScanNumber: p.PrecursorScanNum,
                    Ms2ScanNumber: p.Ms2ScanNumber,
                    PrecursorCharge: p.PrecursorCharge,
                    MonoisotopicMass: p.MonoisotopicMass,
                    RetentionTime: p.RetentionTime,
                    Score: p.Score,
                    FullSequence: p.FullSequence,
                    Accession: p.Accession,
                    Identifier: BuildIdentifier(p),
                    PrecursorIntensity: p.PrecursorIntensity))
                .OrderByDescending(r => r.Score)
                .ThenBy(r => r.RetentionTime)
                .ToArray();

            Console.WriteLine($"Loaded qualified records from {psmTsvPath}: {records.Length} (q<={qValueThreshold})");
            return records;
        }

        private sealed record FitBatch(
            FittedProteoform[] Fits,
            ProteoformGroundTruth[] Truths,
            MmResultRecord[] Records,
            int MinCharge,
            int MaxCharge,
            double SigmaMz,
            IPeakWidthModel WidthModel);

        /// <summary>
        /// Fits σ_m = k·(m/z)^1.5 from the pooled per-isotopologue width measurements of every
        /// record, and re-fits every abundance under it so the fitted numbers describe the file
        /// that will actually be rendered.
        /// </summary>
        /// <remarks>
        /// Returns the median-σ constant model when the pooled fit has too little to work with. The
        /// free-slope diagnostic is printed rather than assumed: if it does not come out near 1.5
        /// the data is saying something about the instrument or about the estimator, and that is
        /// worth seeing before trusting the shipped fit.
        /// </remarks>
        private static (IPeakWidthModel Model, FittedProteoform[] Fits) FitSharedPeakWidth(
            FittedProteoform[] fits,
            ProteoformGroundTruth[] truths,
            double medianSigmaMz,
            bool enabled)
        {
            var constant = new ConstantPeakWidth(medianSigmaMz);
            if (!enabled)
            {
                Console.WriteLine($"Peak width model: {constant} (set by MZLIB_TOPDOWN_SIM_CONSTANT_PEAK_WIDTH)");
                return (constant, fits);
            }

            IPeakWidthModel widthModel;
            double? explicitK = GetExplicitPeakWidthK();
            if (explicitK.HasValue)
            {
                widthModel = new OrbitrapPeakWidth(explicitK.Value);
                Console.WriteLine($"Peak width model: {widthModel} (supplied via MZLIB_TOPDOWN_SIM_PEAK_WIDTH_K)");
            }
            else
            {
                var measurements = fits.SelectMany(f => f.WidthMeasurements).ToArray();
                var widthFit = new PeakWidthModelFitter().Fit(measurements);
                if (widthFit is null)
                {
                    int coarse = fits.Sum(f => f.WidthWindowsTooCoarselySampled);
                    Console.WriteLine($"Peak width model: {constant}");
                    Console.WriteLine($"  Refused to fit a width law: {measurements.Length} usable width measurements" +
                                      (coarse > 0 ? $", with {coarse} windows rejected as too coarsely sampled to be peak shapes." : "."));
                    Console.WriteLine("  This is what a centroided source looks like — the peaks were reduced to positions");
                    Console.WriteLine("  before mzLib saw them, so sigma_m is not recoverable from this file. Supply k via");
                    Console.WriteLine("  MZLIB_TOPDOWN_SIM_PEAK_WIDTH_K to simulate m/z-dependent widths anyway.");
                    return (constant, fits);
                }

                Console.WriteLine($"Peak width fit: {widthFit.Model}");
                Console.WriteLine($"  measurements: {widthFit.MeasurementsUsed}/{widthFit.MeasurementCount} used, " +
                                  $"median sigma {widthFit.MedianSigmaMz:F6} at median m/z {widthFit.MedianMz:F2}");
                Console.WriteLine($"  m/z span covered: {widthFit.MinMz:F1} - {widthFit.MaxMz:F1}");
                Console.WriteLine($"  free-slope diagnostic: {widthFit.FreeSlope:F3} +/- {widthFit.FreeSlopeStandardError:F3} " +
                                  $"(expected ~{OrbitrapPeakWidth.Exponent}; " +
                                  $"{widthFit.SlopeDeviationInStandardErrors:F1} standard errors away)");

                // The slope is the whole reason it is fitted free rather than assumed. If the data
                // disagrees with the width law by more than a few standard errors, the measurements
                // are describing something other than instrument peak width, and shipping k from
                // them would put a confident wrong number into every simulated spectrum.
                if (widthFit.ContradictsWidthLaw(MaxSlopeDeviationInStandardErrors))
                {
                    Console.WriteLine($"Peak width model: {constant}");
                    Console.WriteLine("  Refused the fitted width law: the free slope disagrees with (m/z)^1.5 by more than");
                    Console.WriteLine($"  {MaxSlopeDeviationInStandardErrors} standard errors. A slope near 1 is the signature of width measured from");
                    Console.WriteLine("  centroided input: the extraction window is half the isotopologue spacing, 0.5/z, and");
                    Console.WriteLine("  m/z ~ M/z, so centroid scatter fills a window proportional to m/z. Supply k via");
                    Console.WriteLine("  MZLIB_TOPDOWN_SIM_PEAK_WIDTH_K to simulate m/z-dependent widths anyway.");
                    return (constant, fits);
                }

                if (widthFit.Clamped)
                    Console.WriteLine("  WARNING: the fitted width law hit its plausibility clamp.");

                widthModel = widthFit.Model;
                Console.WriteLine($"Peak width model: {widthModel}");
            }

            Console.WriteLine($"  sigma at m/z 600 / 1000 / 1500: " +
                              $"{widthModel.SigmaAt(600):F5} / {widthModel.SigmaAt(1000):F5} / {widthModel.SigmaAt(1500):F5}");

            // Peak height goes as 1/sigma, so an abundance fitted under this record's own sigma and
            // then rendered under a shared model is off by their ratio. Re-fitting closes that gap.
            var abundanceFitter = new AbundanceFitter();
            var refitted = new FittedProteoform[fits.Length];
            int notRefitted = 0;
            for (int i = 0; i < fits.Length; i++)
            {
                try
                {
                    var abundance = abundanceFitter.Fit(
                        truths[i], widthModel, fits[i].Model.RtProfile, fits[i].Model.ChargeDistribution);
                    refitted[i] = fits[i] with
                    {
                        Model = fits[i].Model with { Abundance = abundance.Abundance },
                        Residual = abundance.Residual,
                    };
                }
                catch (InvalidOperationException)
                {
                    // The model predicts zero everywhere for this record, so there is nothing to fit
                    // against. Its abundance stays as fitted under its own per-record sigma, which
                    // means it is the one case where the fitting/rendering mismatch survives.
                    refitted[i] = fits[i];
                    notRefitted++;
                }
            }

            if (notRefitted > 0)
                Console.WriteLine($"  {notRefitted}/{fits.Length} records could not be re-fitted under the shared width " +
                                  "model and keep an abundance fitted under their own sigma.");

            return (widthModel, refitted);
        }

        /// <summary>
        /// Extracts and fits every record. No global refit is applied here — that is a joint fit
        /// over a chosen set of models, so it belongs to whoever assembles that set.
        /// </summary>
        /// <remarks>
        /// The record loop runs in parallel. Each iteration builds its own IsotopeEnvelopeKernel
        /// and its own fitters, and GroundTruthExtractor is read-only after construction, so
        /// nothing is shared across threads. Results are written into preallocated slots and then
        /// compacted in record order, so the surviving set, the model ordering and the median sigma
        /// are all independent of how the work happened to be scheduled.
        /// </remarks>
        private static FitBatch FitProteoforms(
            IReadOnlyList<IdentifiedSpecies> species,
            GroundTruthExtractor extractor,
            double rtHalfWidth,
            bool fitPeakWidthModel = true)
        {
            var records = species.Select(s => s.Anchor).ToArray();
            var fits = new FittedProteoform[species.Count];
            var truths = new ProteoformGroundTruth[species.Count];
            var chargeRanges = new (int Min, int Max)[species.Count];
            int completed = 0;
            double minSamplesPerSigma = GetMinimumSamplesPerSigma();

            Parallel.For(0, species.Count, i =>
            {
                var record = species[i].Anchor;
                int minCharge = species[i].MinCharge;
                int maxCharge = species[i].MaxCharge;
                if (minCharge > maxCharge)
                    return;

                var truth = extractor.Extract(record.MonoisotopicMass, record.RetentionTime, rtHalfWidth, minCharge, maxCharge);

                FittedProteoform fit;
                try
                {
                    fit = new ParameterFitter(widthFitter: new EnvelopeWidthFitter(
                            fallbackSigmaMz: 0.012,
                            minSamplesPerSigma: minSamplesPerSigma))
                        .Fit(truth, record.Identifier);
                }
                catch (InvalidOperationException)
                {
                    return;
                }

                if (double.IsNaN(fit.Model.Abundance) || fit.Model.Abundance <= 0)
                    return;

                fits[i] = fit;
                truths[i] = truth;
                chargeRanges[i] = (minCharge, maxCharge);

                int done = Interlocked.Increment(ref completed);
                if (done % 100 == 0)
                    Console.WriteLine($"Fitted {done}/{species.Count} records");
            });

            var keptFits = new List<FittedProteoform>(species.Count);
            var keptTruths = new List<ProteoformGroundTruth>(species.Count);
            var keptRecords = new List<MmResultRecord>(species.Count);

            int minChargeGlobal = int.MaxValue;
            int maxChargeGlobal = int.MinValue;

            for (int i = 0; i < species.Count; i++)
            {
                if (fits[i] is null)
                    continue;

                keptFits.Add(fits[i]);
                keptTruths.Add(truths[i]);
                keptRecords.Add(records[i]);
                minChargeGlobal = Math.Min(minChargeGlobal, chargeRanges[i].Min);
                maxChargeGlobal = Math.Max(maxChargeGlobal, chargeRanges[i].Max);
            }

            if (minChargeGlobal == int.MaxValue)
                minChargeGlobal = 2;
            if (maxChargeGlobal == int.MinValue)
                maxChargeGlobal = 80;

            var sigmaCandidates = keptFits
                .Select(f => f.SigmaMz)
                .Where(s => !double.IsNaN(s) && !double.IsInfinity(s) && s > 0)
                .OrderBy(s => s)
                .ToArray();
            double sigmaMz = sigmaCandidates.Length == 0 ? 0.012 : sigmaCandidates[sigmaCandidates.Length / 2];

            var (widthModel, fitsUnderSharedWidth) = FitSharedPeakWidth(
                keptFits.ToArray(), keptTruths.ToArray(), sigmaMz, fitPeakWidthModel);

            return new FitBatch(
                fitsUnderSharedWidth,
                keptTruths.ToArray(),
                keptRecords.ToArray(),
                minChargeGlobal,
                maxChargeGlobal,
                sigmaMz,
                widthModel);
        }

        /// <summary>
        /// Runs the joint non-negative abundance refit over one set of fitted proteoforms and
        /// returns their models. The refit is conditioned on exactly the set passed in.
        /// </summary>
        private static ProteoformModel[] ApplyGlobalRefit(
            FittedProteoform[] fits,
            ProteoformGroundTruth[] truths,
            int minCharge,
            int maxCharge,
            IPeakWidthModel widthModel,
            bool useGlobalAbundanceRefit,
            int globalRefitMaxModels,
            string label)
        {
            if (!useGlobalAbundanceRefit || fits.Length <= 1)
                return fits.Select(f => f.Model).ToArray();

            if (fits.Length > globalRefitMaxModels)
            {
                Console.WriteLine($"Skipping global abundance refit for {label}: {fits.Length} models exceeds max {globalRefitMaxModels}.");
                return fits.Select(f => f.Model).ToArray();
            }

            var sw = Stopwatch.StartNew();
            var refitter = new GlobalAbundanceRefitter(new GlobalAbundanceRefitOptions(
                MaxIterations: 8,
                ConvergenceTolerance: 1e-3,
                MinimumAbundance: 0,
                Verbose: true));

            var refitResult = refitter.Refit(fits, truths, minCharge, maxCharge, widthModel);
            Console.WriteLine($"Global abundance refit ({label}) over {fits.Length} models took {sw.Elapsed}");
            Console.WriteLine($"  iterations: {refitResult.IterationsCompleted}, converged: {refitResult.Converged}");
            Console.WriteLine($"  unexplained energy fraction: {refitResult.InitialResidualFraction:G6} -> {refitResult.FinalResidualFraction:G6}");

            // The two numbers describe opposite failures: signal we could not explain, and signal we
            // predicted where the instrument recorded nothing. The unexplained fraction is over
            // distinct experimental peaks, so a peak claimed by several overlapping proteoforms
            // counts once. The overpredicted one is per sample — an unmatched sample has no peak to
            // key on — so it still scales with crowding and is not comparable across runs.
            var final = refitResult.FinalResiduals;
            if (final is not null)
            {
                Console.WriteLine($"  overpredicted energy fraction (per sample, not deduplicated): " +
                                  $"{refitResult.InitialResiduals!.OverpredictedFraction:G6} -> {final.OverpredictedFraction:G6}");
                Console.WriteLine($"  distinct observed peaks: {final.DistinctPeaks}, samples with no matching peak: {final.UnmatchedSamples}");
            }

            return refitResult.FittedProteoforms.Select(f => f.Model).ToArray();
        }

        /// <summary>Species-level grouping; see <see cref="SpeciesGrouper"/>.</summary>
        private static IdentifiedSpecies[] DeduplicateBySpecies(IReadOnlyList<MmResultRecord> records) =>
            SpeciesGrouper.Group(records);

        /// <summary>
        /// Each record as its own species, extracted over its precursor charge ± 2. What the
        /// pipeline did before species-level deduplication, kept for MZLIB_TOPDOWN_SIM_NO_DEDUP.
        /// </summary>
        private static IdentifiedSpecies[] AsSpecies(IEnumerable<MmResultRecord> records) =>
            SpeciesGrouper.Ungrouped(records);

        private static string BuildIdentifier(PsmFromTsv psm)
        {
            if (!string.IsNullOrWhiteSpace(psm.Accession))
                return $"{psm.Accession}:{psm.FullSequence}:{psm.Ms2ScanNumber}";

            return $"{psm.FileNameWithoutExtension}:{psm.Ms2ScanNumber}";
        }


        private static string ResolveLocalPath(string preferredPath)
        {
            if (File.Exists(preferredPath) || Directory.Exists(preferredPath))
                return preferredPath;

            if (preferredPath.Length >= 3 && preferredPath[1] == ':' && (preferredPath[2] == '\\' || preferredPath[2] == '/'))
            {
                char drive = char.ToLowerInvariant(preferredPath[0]);
                string remainder = preferredPath.Substring(3).Replace('\\', '/');
                string wslPath = $"/mnt/{drive}/{remainder}";
                if (File.Exists(wslPath) || Directory.Exists(wslPath))
                    return wslPath;
            }

            return preferredPath;
        }

    }
}
