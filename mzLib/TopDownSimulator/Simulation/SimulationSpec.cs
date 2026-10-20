#nullable enable
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;
using System.Text.Json.Serialization;
using TopDownSimulator.Model;
using TopDownSimulator.Noise;

namespace TopDownSimulator.Simulation;

/// <summary>
/// SimSpec v1 — the self-describing simulation input the Rust <c>sim-core</c> port reads, and the
/// C# simulator emits. It supersedes the <c>.groundtruth.tsv</c> as the cross-language input because
/// that TSV omits the scan-time axis and the noise/reduction config (see
/// <c>agent_info/TopDown-SimSpec-Contract-v1.md</c> in the ProteoStar repo for the full contract).
/// </summary>
/// <remarks>
/// Serialized with <see cref="System.Text.Json"/> in camelCase + indented, matching the project's
/// existing JSON idiom (<see cref="Noise.AcquisitionProfile"/>). Doubles use .NET's shortest
/// round-trippable form, which <c>serde_json</c> parses back to the identical <c>f64</c>.
/// </remarks>
public sealed record SimSpec(
    int SchemaVersion,
    SimInstrument Instrument,
    SimChargeRange ChargeRange,
    SimPeakWidthModel PeakWidthModel,
    double[] ScanTimes,
    IReadOnlyList<SimProteoform> Proteoforms)
{
    public const int CurrentSchemaVersion = 1;

    public SimGenerator? Generator { get; init; }
    public SimSourceRun? SourceRun { get; init; }
    public SimReduction? Reduction { get; init; }
    public SimNoise? Noise { get; init; }
    public IReadOnlyList<SimProteoform> Background { get; init; } = Array.Empty<SimProteoform>();

    private static readonly JsonSerializerOptions JsonOptions = new()
    {
        PropertyNamingPolicy = JsonNamingPolicy.CamelCase,
        DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull,
        WriteIndented = true
    };

    /// <summary>Builds a spec from the same arguments the simulator's entry points take.</summary>
    public static SimSpec FromModels(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        double[] scanTimes,
        ScanReductionOptions? reduction = null,
        SimNoise? noise = null,
        IReadOnlyList<ProteoformModel>? background = null,
        SimGenerator? generator = null,
        SimSourceRun? sourceRun = null)
    {
        if (proteoforms is null) throw new ArgumentNullException(nameof(proteoforms));
        if (scanTimes is null) throw new ArgumentNullException(nameof(scanTimes));

        return new SimSpec(
            SchemaVersion: CurrentSchemaVersion,
            Instrument: SimInstrument.Default,
            ChargeRange: new SimChargeRange(minCharge, maxCharge),
            PeakWidthModel: SimPeakWidthModel.FromModel(widthModel),
            ScanTimes: scanTimes,
            Proteoforms: proteoforms.Select(SimProteoform.FromModel).ToArray())
        {
            Generator = generator,
            SourceRun = sourceRun,
            Reduction = reduction is null ? null : SimReduction.FromOptions(reduction),
            Noise = noise,
            Background = (background ?? Array.Empty<ProteoformModel>())
                .Select(SimProteoform.FromModel).ToArray()
        };
    }

    public void Write(string path) => File.WriteAllText(path, JsonSerializer.Serialize(this, JsonOptions));

    public static SimSpec Read(string path)
    {
        var spec = JsonSerializer.Deserialize<SimSpec>(File.ReadAllText(path), JsonOptions)
            ?? throw new InvalidDataException($"{path} did not deserialize to a SimSpec.");
        if (spec.SchemaVersion != CurrentSchemaVersion)
            throw new InvalidDataException(
                $"{path} is SimSpec v{spec.SchemaVersion}; this build reads v{CurrentSchemaVersion}.");
        return spec;
    }

    /// <summary>
    /// Reconstructs the arguments a <see cref="Simulator"/> entry point needs. Noise reconstruction
    /// is deferred (v1 goldens run <c>noise: null</c>); callers inspect <see cref="Noise"/> directly.
    /// </summary>
    public SimulationInputs ToSimulationInputs() => new(
        Proteoforms.Select(p => p.ToModel()).ToArray(),
        Background.Select(p => p.ToModel()).ToArray(),
        ChargeRange.Min,
        ChargeRange.Max,
        PeakWidthModel.ToModel(),
        ScanTimes,
        Reduction?.ToOptions());
}

/// <summary>The decoded, model-level view of a <see cref="SimSpec"/>.</summary>
public sealed record SimulationInputs(
    ProteoformModel[] Proteoforms,
    ProteoformModel[] Background,
    int MinCharge,
    int MaxCharge,
    IPeakWidthModel WidthModel,
    double[] ScanTimes,
    ScanReductionOptions? Reduction);

public sealed record SimGenerator(string Tool, string? MzLibCommit = null, string? WrittenUtc = null)
{
    public static SimGenerator Now(string? mzLibCommit = null) =>
        new("TopDownSimulator", mzLibCommit, DateTime.UtcNow.ToString("O"));
}

public sealed record SimSourceRun(string? FileName = null, string? Description = null);

public sealed record SimInstrument(string Polarity, string Analyzer, int MsLevel)
{
    public static SimInstrument Default => new("Positive", "Orbitrap", 1);
}

public sealed record SimChargeRange(int Min, int Max);

/// <summary>Discriminated by <see cref="Kind"/>: <c>"Orbitrap"</c> (k) or <c>"Constant"</c> (sigma).</summary>
public sealed record SimPeakWidthModel(string Kind, double? K = null, double? Sigma = null)
{
    public static SimPeakWidthModel FromModel(IPeakWidthModel model) => model switch
    {
        OrbitrapPeakWidth o => new SimPeakWidthModel("Orbitrap", K: o.K),
        ConstantPeakWidth c => new SimPeakWidthModel("Constant", Sigma: c.SigmaMz),
        _ => throw new NotSupportedException($"SimSpec cannot represent width model {model.GetType().Name}.")
    };

    public IPeakWidthModel ToModel() => Kind switch
    {
        "Orbitrap" => new OrbitrapPeakWidth(K ?? throw Missing("k")),
        "Constant" => new ConstantPeakWidth(Sigma ?? throw Missing("sigma")),
        _ => throw new InvalidDataException($"Unknown peak width model kind '{Kind}'.")
    };

    private InvalidDataException Missing(string field) =>
        new($"Peak width model kind '{Kind}' requires '{field}'.");
}

public sealed record SimReduction(double RelativeIntensityThreshold, double AbsoluteIntensityFloor)
{
    public static SimReduction FromOptions(ScanReductionOptions o) =>
        new(o.RelativeIntensityThreshold, o.MinimumIntensity);

    public ScanReductionOptions ToOptions() => new()
    {
        RelativeIntensityThreshold = RelativeIntensityThreshold,
        MinimumIntensity = AbsoluteIntensityFloor
    };
}

/// <summary>
/// The full, live-editable noise configuration (SimSpec contract §2a). Round-trips to the C# noise
/// models so a spec authored in the UI reconstructs the exact <see cref="NoiseFloorModel"/> and
/// <see cref="PeakJitterModel"/> the simulator runs. The per-scan-conditioned floor
/// (<see cref="ScanNoiseConditions"/>) and its bin tables are not v1 and are not represented here.
/// </summary>
public sealed record SimNoise(long Seed, SimNoiseFloor Floor, SimNoiseJitter Jitter)
{
    public static SimNoise FromModels(NoiseFloorModel floor, PeakJitterModel? jitter, long seed)
    {
        if (floor is null) throw new ArgumentNullException(nameof(floor));
        return new SimNoise(seed, SimNoiseFloor.FromModel(floor), SimNoiseJitter.FromModel(jitter));
    }

    /// <summary>The measured Orbitrap floor + jitter defaults, enabled.</summary>
    public static SimNoise OrbitrapDefaults(long seed = 0) =>
        FromModels(new NoiseFloorModel(), PeakJitterModel.Orbitrap, seed);

    public NoiseFloorModel ToFloorModel() => Floor.ToModel();

    /// <summary>The jitter model, or null when disabled (no jitter is applied).</summary>
    public PeakJitterModel? ToJitterModel() => Jitter.ToModel();
}

/// <summary>Added FT noise floor — mirrors <see cref="NoiseFloorModel"/>'s scalar knobs.</summary>
public sealed record SimNoiseFloor(
    double NoiseLevelAtReferenceMz,
    double DensityScale,
    double LevelDispersion,
    SimSignalToNoise SignalToNoise)
{
    public static SimNoiseFloor FromModel(NoiseFloorModel m) => new(
        m.NoiseLevelAtReferenceMz, m.DensityScale, m.LevelDispersion,
        new SimSignalToNoise(m.SignalToNoise.Threshold, m.SignalToNoise.MeanExcess));

    public NoiseFloorModel ToModel() => new(
        noiseLevelAtReferenceMz: NoiseLevelAtReferenceMz,
        densityScale: DensityScale,
        signalToNoise: SignalToNoise.ToModel(),
        levelDispersion: LevelDispersion);
}

/// <summary>S/N = Threshold + Exp(MeanExcess) — mirrors <see cref="SignalToNoiseDistribution"/>.</summary>
public sealed record SimSignalToNoise(double Threshold, double MeanExcess)
{
    public SignalToNoiseDistribution ToModel() => new(Threshold, MeanExcess);
}

/// <summary>Scan-to-scan measurement error — mirrors <see cref="PeakJitterModel"/>.</summary>
public sealed record SimNoiseJitter(
    bool Enabled,
    double MzPpmAtUnitSignalToNoise,
    double MzPpmFloor,
    double CommonModeMzPpm,
    double LogIntensitySigmaAtUnitSignalToNoise,
    double LogIntensitySigmaFloor,
    double CommonModeLogIntensitySigma,
    bool DropPeaksFallingBelowFloor)
{
    public static readonly SimNoiseJitter Disabled = new(false, 0, 0, 0, 0, 0, 0, false);

    public static SimNoiseJitter FromModel(PeakJitterModel? m) =>
        m is null || !m.IsActive
            ? Disabled
            : new SimNoiseJitter(
                true, m.MzPpmAtUnitSignalToNoise, m.MzPpmFloor, m.CommonModeMzPpm,
                m.LogIntensitySigmaAtUnitSignalToNoise, m.LogIntensitySigmaFloor,
                m.CommonModeLogIntensitySigma, m.DropPeaksFallingBelowFloor);

    public PeakJitterModel? ToModel() => !Enabled
        ? null
        : new PeakJitterModel
        {
            MzPpmAtUnitSignalToNoise = MzPpmAtUnitSignalToNoise,
            MzPpmFloor = MzPpmFloor,
            CommonModeMzPpm = CommonModeMzPpm,
            LogIntensitySigmaAtUnitSignalToNoise = LogIntensitySigmaAtUnitSignalToNoise,
            LogIntensitySigmaFloor = LogIntensitySigmaFloor,
            CommonModeLogIntensitySigma = CommonModeLogIntensitySigma,
            DropPeaksFallingBelowFloor = DropPeaksFallingBelowFloor,
        };
}

public sealed record SimProteoform(
    double MonoisotopicMass,
    double Abundance,
    SimRtProfile Rt,
    SimChargeDistribution Charge,
    string? Identifier = null)
{
    public static SimProteoform FromModel(ProteoformModel p)
    {
        if (p.ChargeDistribution is not GaussianChargeDistribution g)
            throw new NotSupportedException(
                $"SimSpec v1 only serializes Gaussian charge distributions, not {p.ChargeDistribution.GetType().Name}.");

        return new SimProteoform(
            p.MonoisotopicMass,
            p.Abundance,
            new SimRtProfile(p.RtProfile.Mu, p.RtProfile.Sigma, p.RtProfile.Tau),
            new SimChargeDistribution("Gaussian", g.MuZ, g.SigmaZ),
            string.IsNullOrEmpty(p.Identifier) ? null : p.Identifier);
    }

    public ProteoformModel ToModel() => new(
        MonoisotopicMass,
        Abundance,
        new EmgProfile(Rt.Mu, Rt.Sigma, Rt.Tau),
        Charge.ToModel(),
        Identifier);
}

public sealed record SimRtProfile(double Mu, double Sigma, double Tau);

public sealed record SimChargeDistribution(string Kind, double MuZ, double SigmaZ)
{
    public IChargeStateDistribution ToModel() => Kind switch
    {
        "Gaussian" => new GaussianChargeDistribution(MuZ, SigmaZ),
        _ => throw new InvalidDataException($"Unknown charge distribution kind '{Kind}'.")
    };
}
