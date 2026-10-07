#nullable enable
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using MassSpectrometry;
using Readers;
using TopDownSimulator.Model;
using TopDownSimulator.Noise;

namespace TopDownSimulator.Simulation;

/// <param name="Grid">
/// The profile grid the scans were rendered on, or null for a centroided simulation, where every
/// scan has its own peak list rather than samples of one shared axis.
/// </param>
public sealed record SimulationResult(RasterizedScanGrid? Grid, MsDataScan[] Scans, GenericMsDataFile DataFile);

/// <summary>
/// Outcome of writing a simulation to disk. <see cref="PeakCount"/> is counted after reduction,
/// so it reflects what actually landed in the file.
/// </summary>
public sealed record SimulationExportResult(
    string MzmlPath,
    string? GroundTruthPath,
    int ScanCount,
    int PeakCount,
    string? FeatureGroundTruthPath = null,
    int FeatureCount = 0,
    NoiseInjectionSummary? Noise = null,
    PrecursorShiftSummary? Shift = null);

/// <summary>
/// What a constant precursor mass shift did to a written run.
/// </summary>
/// <param name="MassShiftDa">The neutral offset applied to every simulated proteoform.</param>
/// <param name="ShiftSidecarPath">Where the original-to-shifted mass table was written.</param>
/// <param name="Ms2ScansWritten">MSn scans that made it into the file with a shifted precursor.</param>
public sealed record PrecursorShiftSummary(
    double MassShiftDa,
    string? ShiftSidecarPath,
    int Ms2ScansWritten,
    int Ms2ScansDroppedWithoutChargeState,
    int Ms2ScansDroppedWithoutPrecursorMz,
    int Ms2ScansDroppedProfileMode,
    int Ms2ScansDroppedOutsideMs1TimeRange,
    int Ms2ScansDroppedWithEmptyPrecursor);

/// <summary>
/// A simulated MS1 run after reduction and noise injection but before it is written, which is the
/// point at which MSn scans can still be interleaved into it.
/// </summary>
/// <param name="ScanTimes">Retention time of each MS1 scan.</param>
/// <param name="SignalMz">
/// Per scan, the m/z of every signal centroid that survived reduction, before jitter and noise.
/// The feature truth quotes positions from these, so it names peaks that are really in the file.
/// </param>
/// <param name="Floors">The detection limit each scan was reduced against.</param>
public sealed record PreparedSimulation(
    MsDataScan[] Ms1Scans,
    double[] ScanTimes,
    double[][] SignalMz,
    IIntensityFloor[] Floors,
    NoiseInjectionSummary? Noise);

/// <summary>
/// High-level Phase 5 entrypoint for simulating synthetic MS1 data and writing it as mzML.
/// </summary>
public sealed class Simulator
{
    private readonly GridRasterizer _rasterizer;
    private readonly ScanBuilder _scanBuilder;

    public Simulator(GridRasterizer? rasterizer = null, ScanBuilder? scanBuilder = null)
    {
        _rasterizer = rasterizer ?? new GridRasterizer();
        _scanBuilder = scanBuilder ?? new ScanBuilder();
    }

    public SimulationResult Simulate(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        double[] scanTimes,
        bool isCentroid = false,
        int pointsPerSigma = 3,
        double mzPaddingInSigmas = 6.0)
    {
        var grid = _rasterizer.Rasterize(proteoforms, minCharge, maxCharge, widthModel, scanTimes, pointsPerSigma, mzPaddingInSigmas);
        var scans = _scanBuilder.BuildMs1Scans(grid, isCentroid: isCentroid);
        var file = _scanBuilder.BuildMsDataFile(scans);
        return new SimulationResult(grid, scans, file);
    }

    public SimulationResult Simulate(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes,
        bool isCentroid = false,
        int pointsPerSigma = 3,
        double mzPaddingInSigmas = 6.0) =>
        Simulate(proteoforms, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz), scanTimes, isCentroid, pointsPerSigma, mzPaddingInSigmas);

    /// <summary>
    /// Simulates centroid spectra: one centroid per local maximum of the summed profile, at that
    /// maximum's height, which is what an instrument reports. See <see cref="ProfileCentroider"/>.
    /// </summary>
    public SimulationResult SimulateCentroided(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        double[] scanTimes)
    {
        var spectra = new ProfileCentroider(proteoforms, minCharge, maxCharge, widthModel).Centroid(scanTimes);
        var scans = _scanBuilder.BuildCentroidedMs1Scans(scanTimes, spectra);
        var file = _scanBuilder.BuildMsDataFile(scans);
        return new SimulationResult(null, scans, file);
    }

    public SimulationResult SimulateCentroided(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes) =>
        SimulateCentroided(proteoforms, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz), scanTimes);

    /// <summary>
    /// Simulates centroid spectra and writes them as mzML, together with two sidecars: the
    /// parameter vector θ_p that generated the run, and the feature-level truth a benchmark scores
    /// against. Both are suppressed by <paramref name="writeGroundTruthSidecar"/> = false.
    /// </summary>
    /// <remarks>
    /// Output is centroided, which is not merely a preference: mzLib's own mzML reader rejects
    /// profile-mode files, so a profile simulation could be written but never read back by mzLib
    /// or MetaMorpheus.
    /// <para>
    /// Supplying <paramref name="noise"/> both injects a sampled noise floor and switches reduction
    /// to that model's detection limit, so signal and noise are thresholded under one rule. Note
    /// that at the measured density this is expensive — roughly 17 000 extra peaks per scan, which
    /// is what a real run carries but which takes the written file from megabytes to hundreds of
    /// megabytes. Use <see cref="NoiseFloorModel.DensityScale"/> to trade fidelity for size.
    /// </para>
    /// </remarks>
    /// <param name="noise">Noise floor to inject, or null for a clean analytical simulation.</param>
    /// <param name="noiseSeed">Seed making an injected noise floor reproducible.</param>
    /// <param name="jitter">
    /// Scan-to-scan measurement error on the true peaks. Requires <paramref name="noise"/>, which
    /// defines the S/N scale the jitter laws are expressed against; set that model's
    /// <see cref="NoiseFloorModel.DensityScale"/> to 0 for jitter without an added noise floor.
    /// Defaults to <see cref="PeakJitterModel.Orbitrap"/> whenever noise is supplied, since a noisy
    /// file with perfectly smooth XICs is not a realistic combination.
    /// </param>
    /// <param name="scanNoise">
    /// One noise model per scan, replacing <paramref name="noise"/>, for a floor whose amplitude and
    /// density follow the run. Build it with <see cref="ScanNoiseConditions.FromSourceScans"/>.
    /// </param>
    /// <param name="background">
    /// Proteoforms rendered into the spectra but left out of both sidecars: analytes the instrument
    /// saw that no identification accounts for. A benchmark is not scored on finding them.
    /// </param>
    /// <param name="scanNoiseFromSignal">
    /// Builds the per-scan noise models from the rendered signal scans, before noise is added, for a
    /// floor that follows the simulated ion load (see <see cref="AutomaticGainControl"/>). Replaces
    /// <paramref name="scanNoise"/>.
    /// </param>
    public SimulationExportResult WriteMzml(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        double[] scanTimes,
        string outputPath,
        ScanReductionOptions? reduction = null,
        bool writeIndexed = true,
        bool writeGroundTruthSidecar = true,
        NoiseFloorModel? noise = null,
        int noiseSeed = 0,
        PeakJitterModel? jitter = null,
        IReadOnlyList<NoiseFloorModel>? scanNoise = null,
        IReadOnlyList<ProteoformModel>? background = null,
        Func<MsDataScan[], IReadOnlyList<NoiseFloorModel>>? scanNoiseFromSignal = null)
    {
        if (string.IsNullOrWhiteSpace(outputPath))
            throw new ArgumentException("An output path is required.", nameof(outputPath));

        var rendered = background is { Count: > 0 } ? proteoforms.Concat(background).ToArray() : proteoforms;
        var prepared = PrepareMs1(
            rendered, minCharge, maxCharge, widthModel, scanTimes, reduction, noise, noiseSeed, jitter, scanNoise,
            scanNoiseFromSignal);

        var scans = prepared.Ms1Scans;
        var scanNumbers = new int[scans.Length];
        for (int s = 0; s < scans.Length; s++)
            scanNumbers[s] = scans[s].OneBasedScanNumber;

        return Write(
            scans, prepared, proteoforms, minCharge, maxCharge, widthModel,
            scanNumbers, outputPath, writeIndexed, writeGroundTruthSidecar, shift: null);
    }

    /// <summary>
    /// Simulates centroided MS1 spectra with every proteoform's monoisotopic mass moved by
    /// <paramref name="massShiftDa"/>, interleaves the MSn scans of <paramref name="sourceMs2Scans"/>
    /// with their isolation windows moved by <paramref name="massShiftDa"/>/z, and writes the result
    /// as one mzML.
    /// </summary>
    /// <remarks>
    /// The MSn peak lists are carried over untouched, so the fragments still describe the unshifted
    /// proteoform while the precursor no longer does — the file reads as an unlocalized modification
    /// of mass <paramref name="massShiftDa"/> on every identification. See
    /// <see cref="PrecursorMassShift"/> for why that is the useful arrangement.
    /// <para>
    /// <paramref name="proteoforms"/> are the <em>unshifted</em> models; the shift is applied here so
    /// that the original masses are still available to write the <c>.massshift.tsv</c> sidecar. The
    /// parameter and feature sidecars describe the file as written, and so quote shifted masses.
    /// </para>
    /// <para>
    /// MSn scans outside the retention-time span of <paramref name="scanTimes"/>, or with no charge
    /// state to convert the neutral shift through, are dropped and counted in
    /// <see cref="SimulationExportResult.Shift"/>.
    /// </para>
    /// </remarks>
    /// <param name="sourceMs2Scans">MSn scans from the run being mimicked, in any order.</param>
    /// <param name="massShiftDa">The neutral offset, in daltons. May be negative.</param>
    /// <param name="dropProfileMs2Scans">
    /// Drop profile-mode MSn scans rather than refusing to write. A profile spectrum anywhere in the
    /// file makes the whole thing unreadable by mzLib, so the default is to fail loudly.
    /// </param>
    public SimulationExportResult WriteShiftedMzml(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        double[] scanTimes,
        IReadOnlyList<MsDataScan> sourceMs2Scans,
        double massShiftDa,
        string outputPath,
        ScanReductionOptions? reduction = null,
        bool writeIndexed = true,
        bool writeGroundTruthSidecar = true,
        NoiseFloorModel? noise = null,
        int noiseSeed = 0,
        PeakJitterModel? jitter = null,
        bool dropProfileMs2Scans = false,
        IReadOnlyList<NoiseFloorModel>? scanNoise = null)
    {
        if (proteoforms is null) throw new ArgumentNullException(nameof(proteoforms));
        if (sourceMs2Scans is null) throw new ArgumentNullException(nameof(sourceMs2Scans));
        if (string.IsNullOrWhiteSpace(outputPath))
            throw new ArgumentException("An output path is required.", nameof(outputPath));

        var shiftedModels = PrecursorMassShift.ApplyToModels(proteoforms, massShiftDa);
        var (shiftedMs2, ms2Summary) = PrecursorMassShift.ApplyToMs2Scans(
            sourceMs2Scans, massShiftDa, dropProfileMs2Scans);

        var prepared = PrepareMs1(
            shiftedModels, minCharge, maxCharge, widthModel, scanTimes, reduction, noise, noiseSeed, jitter, scanNoise);

        var merged = ScanListMerger.Merge(prepared.Ms1Scans, shiftedMs2);

        string? shiftPath = null;
        if (writeGroundTruthSidecar)
        {
            shiftPath = Path.ChangeExtension(outputPath, ".massshift.tsv");
            PrecursorMassShift.WriteSidecar(
                PrecursorMassShift.Describe(proteoforms, massShiftDa), massShiftDa, shiftPath);
        }

        var shift = new PrecursorShiftSummary(
            MassShiftDa: massShiftDa,
            ShiftSidecarPath: shiftPath,
            Ms2ScansWritten: shiftedMs2.Length
                - merged.Ms2ScansDroppedOutsideMs1TimeRange
                - merged.Ms2ScansDroppedWithEmptyPrecursor,
            Ms2ScansDroppedWithoutChargeState: ms2Summary.DroppedWithoutChargeState,
            Ms2ScansDroppedWithoutPrecursorMz: ms2Summary.DroppedWithoutPrecursorMz,
            Ms2ScansDroppedProfileMode: ms2Summary.DroppedProfileMode,
            Ms2ScansDroppedOutsideMs1TimeRange: merged.Ms2ScansDroppedOutsideMs1TimeRange,
            Ms2ScansDroppedWithEmptyPrecursor: merged.Ms2ScansDroppedWithEmptyPrecursor);

        return Write(
            merged.Scans, prepared, shiftedModels, minCharge, maxCharge, widthModel,
            merged.Ms1ScanNumbers, outputPath, writeIndexed, writeGroundTruthSidecar, shift);
    }

    public SimulationExportResult WriteShiftedMzml(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes,
        IReadOnlyList<MsDataScan> sourceMs2Scans,
        double massShiftDa,
        string outputPath,
        ScanReductionOptions? reduction = null,
        bool writeIndexed = true,
        bool writeGroundTruthSidecar = true,
        NoiseFloorModel? noise = null,
        int noiseSeed = 0,
        PeakJitterModel? jitter = null,
        bool dropProfileMs2Scans = false) =>
        WriteShiftedMzml(proteoforms, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz), scanTimes,
            sourceMs2Scans, massShiftDa, outputPath, reduction, writeIndexed, writeGroundTruthSidecar,
            noise, noiseSeed, jitter, dropProfileMs2Scans);

    /// <summary>
    /// Runs the simulation up to the point where the scans are ready to write, so that MSn scans can
    /// still be interleaved into them.
    /// </summary>
    public PreparedSimulation PrepareMs1(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        double[] scanTimes,
        ScanReductionOptions? reduction = null,
        NoiseFloorModel? noise = null,
        int noiseSeed = 0,
        PeakJitterModel? jitter = null,
        IReadOnlyList<NoiseFloorModel>? scanNoise = null,
        Func<MsDataScan[], IReadOnlyList<NoiseFloorModel>>? scanNoiseFromSignal = null)
    {
        var simulation = SimulateCentroided(proteoforms, minCharge, maxCharge, widthModel, scanTimes);
        if (scanNoiseFromSignal is not null)
            scanNoise = scanNoiseFromSignal(simulation.Scans);

        if (scanNoise is not null && scanNoise.Count != scanTimes.Length)
            throw new ArgumentException("There must be one noise model per scan.", nameof(scanNoise));

        reduction ??= new ScanReductionOptions();
        if (noise is not null && scanNoise is null && reduction.Floor is null)
            reduction = reduction with { Floor = new NoiseRelativeFloor(noise) };

        // Per-scan noise brings a per-scan detection limit, unless the caller fixed one.
        IIntensityFloor[] floors;
        if (scanNoise is not null && reduction.Floor is null)
        {
            floors = new IIntensityFloor[scanTimes.Length];
            for (int s = 0; s < floors.Length; s++)
                floors[s] = new NoiseRelativeFloor(scanNoise[s]);
        }
        else
        {
            floors = new IIntensityFloor[scanTimes.Length];
            Array.Fill(floors, SimulatedScanReducer.ComputeFloor(simulation.Scans, reduction));
        }

        var scans = SimulatedScanReducer.Reduce(simulation.Scans, floors);
        var signalMz = new double[scans.Length][];
        for (int s = 0; s < scans.Length; s++)
            signalMz[s] = scans[s].MassSpectrum.XArray;

        NoiseInjectionSummary? noiseSummary = null;
        if (scanNoise is not null || noise is not null)
        {
            var injector = scanNoise is not null
                ? new NoiseInjector(scanNoise, widthModel, noiseSeed,
                    jitter: jitter ?? PeakJitterModel.Orbitrap, scanFloors: floors)
                : new NoiseInjector(noise!, widthModel, noiseSeed,
                    jitter: jitter ?? PeakJitterModel.Orbitrap, floor: floors.Length > 0 ? floors[0] : null);
            (scans, noiseSummary) = injector.Apply(scans);
        }

        return new PreparedSimulation(scans, (double[])scanTimes.Clone(), signalMz, floors, noiseSummary);
    }

    /// <summary>
    /// Writes <paramref name="scans"/> as mzML plus the sidecars, where
    /// <paramref name="ms1ScanNumbers"/> gives the number each MS1 scan carries in the written file.
    /// </summary>
    private SimulationExportResult Write(
        MsDataScan[] scans,
        PreparedSimulation prepared,
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        int[] ms1ScanNumbers,
        string outputPath,
        bool writeIndexed,
        bool writeGroundTruthSidecar,
        PrecursorShiftSummary? shift)
    {
        var file = _scanBuilder.BuildMsDataFile(scans);
        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(file, outputPath, writeIndexed);

        int peakCount = 0;
        foreach (var scan in scans)
            peakCount += scan.MassSpectrum.XArray.Length;

        string? groundTruthPath = null;
        string? featurePath = null;
        int featureCount = 0;
        if (writeGroundTruthSidecar)
        {
            groundTruthPath = Path.ChangeExtension(outputPath, ".groundtruth.tsv");
            WriteGroundTruth(proteoforms, minCharge, maxCharge, widthModel, groundTruthPath);

            var features = FeatureGroundTruth.Build(
                proteoforms, minCharge, maxCharge, widthModel,
                prepared.ScanTimes, ms1ScanNumbers, prepared.SignalMz, prepared.Floors);

            featurePath = Path.ChangeExtension(outputPath, ".features.tsv");
            FeatureGroundTruth.Write(features, featurePath);
            featureCount = features.Length;
        }

        return new SimulationExportResult(
            outputPath, groundTruthPath, scans.Length, peakCount, featurePath, featureCount,
            prepared.Noise, shift);
    }

    public SimulationExportResult WriteMzml(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        double[] scanTimes,
        string outputPath,
        ScanReductionOptions? reduction = null,
        bool writeIndexed = true,
        bool writeGroundTruthSidecar = true,
        NoiseFloorModel? noise = null,
        int noiseSeed = 0,
        PeakJitterModel? jitter = null) =>
        WriteMzml(proteoforms, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz), scanTimes, outputPath,
            reduction, writeIndexed, writeGroundTruthSidecar, noise, noiseSeed, jitter);

    /// <summary>
    /// Writes the parameter vector for every simulated proteoform so a run can be reproduced.
    /// </summary>
    /// <remarks>
    /// This file is not scorable on its own — see <see cref="FeatureGroundTruth"/> for the
    /// feature-level truth. <c>SigmaMz</c> is the width at that proteoform's most abundant charge,
    /// which under a non-constant width law differs between rows; <c>PeakWidthModel</c> records the
    /// law itself.
    /// </remarks>
    public static void WriteGroundTruth(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        string outputPath)
    {
        var sb = new StringBuilder();
        sb.AppendLine(string.Join('\t',
            "Identifier", "MonoisotopicMass", "Abundance",
            "RtMu", "RtSigma", "RtTau",
            "ChargeMu", "ChargeSigma", "MinCharge", "MaxCharge", "SigmaMz", "PeakWidthModel"));

        foreach (var p in proteoforms)
        {
            var gaussian = p.ChargeDistribution as GaussianChargeDistribution;

            sb.AppendLine(string.Join('\t',
                p.Identifier ?? string.Empty,
                Num(p.MonoisotopicMass),
                Num(p.Abundance),
                Num(p.RtProfile.Mu),
                Num(p.RtProfile.Sigma),
                Num(p.RtProfile.Tau),
                gaussian is null ? string.Empty : Num(gaussian.MuZ),
                gaussian is null ? string.Empty : Num(gaussian.SigmaZ),
                minCharge.ToString(CultureInfo.InvariantCulture),
                maxCharge.ToString(CultureInfo.InvariantCulture),
                Num(RepresentativeSigma(p, minCharge, maxCharge, widthModel)),
                widthModel.ToString()));
        }

        File.WriteAllText(outputPath, sb.ToString());
    }

    /// <summary>
    /// Reads back a sidecar written by <see cref="WriteGroundTruth(IReadOnlyList{ProteoformModel}, int, int, IPeakWidthModel, string)"/>,
    /// which is enough to re-simulate the run without refitting.
    /// </summary>
    /// <remarks>
    /// The width law is parsed from its <c>ToString</c> form, so only <see cref="ConstantPeakWidth"/>
    /// and <see cref="OrbitrapPeakWidth"/> round-trip. Charge distributions other than Gaussian are
    /// not written by the sidecar and cannot be read.
    /// </remarks>
    public static (ProteoformModel[] Models, int MinCharge, int MaxCharge, IPeakWidthModel WidthModel) ReadGroundTruth(string path)
    {
        var lines = File.ReadAllLines(path);
        if (lines.Length == 0)
            throw new InvalidDataException($"{path} is empty.");

        var header = lines[0].Split('\t');
        int Column(string name)
        {
            int index = Array.IndexOf(header, name);
            return index >= 0 ? index : throw new InvalidDataException($"{path} has no {name} column.");
        }

        int id = Column("Identifier"), mass = Column("MonoisotopicMass"), abundance = Column("Abundance");
        int mu = Column("RtMu"), sigma = Column("RtSigma"), tau = Column("RtTau");
        int muZ = Column("ChargeMu"), sigmaZ = Column("ChargeSigma");
        int minZ = Column("MinCharge"), maxZ = Column("MaxCharge"), width = Column("PeakWidthModel");

        static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

        var models = new List<ProteoformModel>();
        int minCharge = int.MaxValue, maxCharge = int.MinValue;
        string? widthText = null;
        foreach (var line in lines.Skip(1))
        {
            if (line.Length == 0) continue;
            var f = line.Split('\t');
            models.Add(new ProteoformModel(
                D(f[mass]), D(f[abundance]),
                new EmgProfile(D(f[mu]), D(f[sigma]), D(f[tau])),
                new GaussianChargeDistribution(D(f[muZ]), D(f[sigmaZ])),
                f[id].Length > 0 ? f[id] : null));
            minCharge = Math.Min(minCharge, int.Parse(f[minZ], CultureInfo.InvariantCulture));
            maxCharge = Math.Max(maxCharge, int.Parse(f[maxZ], CultureInfo.InvariantCulture));
            widthText = f[width];
        }

        if (widthText is null)
            throw new InvalidDataException($"{path} holds no proteoforms.");

        return (models.ToArray(), minCharge, maxCharge, ParseWidthModel(widthText));
    }

    private static IPeakWidthModel ParseWidthModel(string text)
    {
        var constant = System.Text.RegularExpressions.Regex.Match(text, @"^Constant\(sigma=([^)]+)\)$");
        if (constant.Success)
            return new ConstantPeakWidth(double.Parse(constant.Groups[1].Value, CultureInfo.InvariantCulture));

        var orbitrap = System.Text.RegularExpressions.Regex.Match(text, @"^Orbitrap\(k=([^)]+)\)$");
        if (orbitrap.Success)
            return new OrbitrapPeakWidth(double.Parse(orbitrap.Groups[1].Value, CultureInfo.InvariantCulture));

        throw new InvalidDataException($"Cannot rebuild the peak width model '{text}'.");
    }

    public static void WriteGroundTruth(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        double sigmaMz,
        string outputPath) =>
        WriteGroundTruth(proteoforms, minCharge, maxCharge, new ConstantPeakWidth(sigmaMz), outputPath);

    /// <summary>σ at the m/z of the proteoform's most abundant charge state.</summary>
    private static double RepresentativeSigma(
        ProteoformModel proteoform, int minCharge, int maxCharge, IPeakWidthModel widthModel)
    {
        int bestCharge = minCharge;
        double bestWeight = double.NegativeInfinity;
        for (int z = minCharge; z <= maxCharge; z++)
        {
            double fz = proteoform.ChargeDistribution.Evaluate(z);
            if (fz > bestWeight)
            {
                bestWeight = fz;
                bestCharge = z;
            }
        }

        var centroids = new IsotopeEnvelopeKernel(proteoform.MonoisotopicMass).CentroidMzs(bestCharge);
        return centroids.Length == 0 ? 0 : widthModel.SigmaAt(centroids[0]);
    }

    private static string Num(double value) => value.ToString("R", CultureInfo.InvariantCulture);
}
