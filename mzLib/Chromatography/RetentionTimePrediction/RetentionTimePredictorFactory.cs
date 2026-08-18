using System.Reflection;
using Chromatography.RetentionTimePrediction.CZE;
using Chromatography.RetentionTimePrediction.SSRCalc;
using MzLibUtil;

namespace Chromatography.RetentionTimePrediction;

/// <summary>
/// Discriminator identifying which concrete retention time predictor to build
/// via <see cref="RetentionTimePredictorFactory.Create"/>.
/// </summary>
/// <remarks>
/// This enum was previously nested inside <see cref="IRetentionTimePredictor"/>
/// (as <c>IRetentionTimePredictor.PredictorType</c>). It was promoted to a top-level
/// type alongside the factory so the interface no longer needs to reference any
/// concrete predictor. Callers that used the old nested name should update to the
/// top-level <see cref="PredictorType"/>.
/// </remarks>
public enum PredictorType
{
    /// <summary>
    /// Krokhin/SSRCalc v3 — analytical hydrophobicity predictor for reversed-phase HPLC.
    /// Lightweight, no external model files or native dependencies.
    /// </summary>
    SSRCalc3,

    /// <summary>
    /// Capillary zone electrophoresis (CZE) migration-time predictor.
    /// Lightweight, no external model files or native dependencies.
    /// </summary>
    CZE,

    /// <summary>
    /// Chronologer deep-learning retention time predictor (TorchSharp-backed).
    /// Supplied by the separate <c>Chromatography.Chronologer</c> assembly
    /// (<c>mzLib.Chronologer</c> package), which pulls in TorchSharp and large native
    /// runtime binaries; reference it only when model-based prediction is actually required.
    /// </summary>
    Chronologer
}

/// <summary>
/// Factory for <see cref="IRetentionTimePredictor"/> instances. Given a
/// <see cref="PredictorType"/>, returns a fully initialized predictor ready for use.
/// </summary>
/// <remarks>
/// <para>
/// <b>Why this is a standalone class rather than a static method on the interface:</b>
/// when construction lived as a <c>static</c> member of <see cref="IRetentionTimePredictor"/>,
/// any project that merely referenced the interface was forced to compile against — and ship —
/// the full predictor dependency graph, even if it only described or consumed predictors
/// (e.g. dependency-injection hosts, tests with fakes, lightweight tooling).
/// </para>
/// <para>
/// <b>Why it stays in <c>Chromatography</c> rather than moving to
/// <c>Chromatography.Chronologer</c>:</b> two of the three predictors are pure C# with no
/// native dependency. Moving the factory and <see cref="PredictorType"/> into the TorchSharp
/// assembly would have gated the documented construction path for SSRCalc3 and CZE behind the
/// 318 MB opt-in package, which inverts the goal of splitting it out in the first place.
/// Chronologer is instead reached through <see cref="ChronologerFactory"/>, so this assembly
/// never references the Torch one — the dependency runs the other way.
/// </para>
/// <para>
/// This mirrors how <c>Deconvoluter.CreateAlgorithm</c> reaches
/// <c>FromFileDeconvolutionAlgorithm</c>, which likewise lives in an assembly
/// <c>MassSpectrometry</c> does not reference: a polymorphic hook first, the in-project
/// switch as the fallback, and an <see cref="MzLibException"/> naming the missing project
/// when neither is available.
/// </para>
/// <para>
/// <b>Extending:</b> to add a new predictor that lives in this assembly, add a value to
/// <see cref="PredictorType"/> and a branch in <see cref="Create"/>. To add one that lives
/// elsewhere, follow the Chronologer pattern.
/// </para>
/// </remarks>
public static class RetentionTimePredictorFactory
{
    /// <summary>
    /// Assembly-qualified name of the Chronologer predictor, used to auto-wire
    /// <see cref="PredictorType.Chronologer"/> when the <c>mzLib.Chronologer</c> package is
    /// present but nothing has called <see cref="RegisterChronologerFactory"/>.
    /// </summary>
    /// <remarks>
    /// Resolved by name rather than by reference so that <c>Chromatography</c> does not
    /// depend on <c>Chromatography.Chronologer</c>. Both assemblies ship from the same
    /// release tag and land in the same output folder, so the load succeeds whenever the
    /// package is referenced. A rename on either side is caught by
    /// <c>ChronologerFactoryTests.ChronologerTypeNameResolves</c>, which fails the build's
    /// test pass rather than deferring the break to a consumer's runtime.
    /// </remarks>
    internal const string ChronologerPredictorTypeName =
        "Chromatography.RetentionTimePrediction.Chronologer.ChronologerRetentionTimePredictor, Chromatography.Chronologer";

    private static Func<IRetentionTimePredictor>? _chronologerFactory;

    /// <summary>
    /// Seam over <see cref="Type.GetType(string, bool)"/> so the "package not referenced" path can
    /// be exercised. The test assembly references every project, so the type always resolves there
    /// and the <see cref="MzLibException"/> branch would otherwise be unreachable.
    /// </summary>
    internal static Func<string, Type?> ChronologerTypeResolver { get; set; } =
        name => Type.GetType(name, throwOnError: false);

    /// <summary>
    /// Registers the delegate used to build <see cref="PredictorType.Chronologer"/>.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Optional. When the <c>mzLib.Chronologer</c> package is referenced the predictor is
    /// found by name without any registration, so most callers never need this. It exists for
    /// composition roots that want construction to be explicit rather than reflective — for
    /// instance to inject a pre-warmed or differently-configured Chronologer instance — and it
    /// is how the "package not installed" path is exercised in tests, which necessarily
    /// reference every assembly.
    /// </para>
    /// <para>
    /// Pass <c>null</c> to clear a previously registered delegate and fall back to name
    /// resolution. Registration is process-wide and is not thread-safe against a concurrent
    /// <see cref="Create"/>; set it during startup.
    /// </para>
    /// </remarks>
    /// <param name="factory">Delegate producing a new Chronologer predictor, or <c>null</c> to clear.</param>
    public static void RegisterChronologerFactory(Func<IRetentionTimePredictor>? factory) =>
        _chronologerFactory = factory;

    /// <summary>
    /// Creates a new <see cref="IRetentionTimePredictor"/> of the specified
    /// <paramref name="type"/>, constructed with default parameters.
    /// </summary>
    /// <param name="type">The predictor variant to instantiate.</param>
    /// <returns>A fully initialized predictor ready for prediction calls.</returns>
    /// <remarks>
    /// <para>
    /// Each call returns a fresh instance — instances are not cached or pooled.
    /// </para>
    /// <para>
    /// The caller owns the lifetime of the returned predictor and is responsible
    /// for disposing predictors that implement <see cref="IDisposable"/>. In
    /// particular, the Chronologer predictor holds an unmanaged TorchSharp model and
    /// must be disposed; SSRCalc3 and CZE hold no unmanaged resources but should still
    /// be disposed for consistency.
    /// </para>
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException">
    /// Thrown when <paramref name="type"/> is not a defined <see cref="PredictorType"/>
    /// value (for example, when an out-of-range integer is cast to the enum).
    /// </exception>
    /// <exception cref="MzLibException">
    /// Thrown for <see cref="PredictorType.Chronologer"/> when the <c>mzLib.Chronologer</c>
    /// package is not referenced.
    /// </exception>
    public static IRetentionTimePredictor Create(PredictorType type) => type switch
    {
        PredictorType.SSRCalc3 => new SSRCalc3RetentionTimePredictor(),
        PredictorType.CZE => new CZERetentionTimePredictor(),
        PredictorType.Chronologer => CreateChronologer(),
        _ => throw new ArgumentOutOfRangeException(nameof(type), type, null)
    };

    /// <summary>
    /// Builds the Chronologer predictor from the registered delegate if there is one, and
    /// otherwise from the <c>Chromatography.Chronologer</c> assembly resolved by name.
    /// </summary>
    private static IRetentionTimePredictor CreateChronologer()
    {
        if (_chronologerFactory is not null)
            return _chronologerFactory();

        Type? predictorType = ChronologerTypeResolver(ChronologerPredictorTypeName);
        if (predictorType is null)
            throw new MzLibException(
                "Chronologer retention time prediction requires the mzLib.Chronologer package, which " +
                "carries TorchSharp and the libtorch native runtime. Add " +
                "<PackageReference Include=\"mzLib.Chronologer\" /> (or a ProjectReference to " +
                "Chromatography.Chronologer), register a factory with " +
                "RetentionTimePredictorFactory.RegisterChronologerFactory, or use " +
                "PredictorType.SSRCalc3 or PredictorType.CZE, which need no native dependencies.");

        // OptionalParamBinding rather than the one-argument Activator.CreateInstance overload:
        // ChronologerRetentionTimePredictor's only constructor takes two *optional* parameters, and
        // optional parameters do not synthesise a parameterless constructor, so the simple overload
        // would throw MissingMethodException. ChronologerFactoryTests covers this path, so a change
        // to that constructor fails the test pass rather than a consumer's first prediction.
        object instance = Activator.CreateInstance(
            predictorType,
            BindingFlags.Instance | BindingFlags.Public | BindingFlags.CreateInstance | BindingFlags.OptionalParamBinding,
            binder: null,
            args: [],
            culture: null)!;

        return (IRetentionTimePredictor)instance;
    }
}
