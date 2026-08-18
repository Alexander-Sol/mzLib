using Omics.RetentionTimePrediction;

// This file intentionally declares a namespace that does not match its assembly. It preserves the
// pre-Chromatography type locations for downstream consumers (notably MetaMorpheus) that still
// reference Proteomics.RetentionTimePrediction.SSRCalc3. It lives in the Chromatography project
// rather than Proteomics so that Proteomics no longer has to reference Chromatography — that edge
// was what dragged TorchSharp and libtorch into the dependency graph of every project built on
// Proteomics, which is to say almost all of them.
namespace Proteomics.RetentionTimePrediction;

/// <summary>
/// Backwards-compatibility shims for the type layout that predated the Chromatography namespace.
/// </summary>
/// <remarks>
/// Prefer <see cref="Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3"/> and its
/// <c>ScoreSequence(string)</c> overload directly. These shims are scheduled for removal in the
/// next major version.
/// </remarks>
[Obsolete("Use Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3.ScoreSequence(string) directly. This shim will be removed in the next major version.")]
public static class ChromatographyExtensions
{
    /// <summary>
    /// Scores a peptide's base sequence with SSRCalc3.
    /// </summary>
    /// <remarks>
    /// The parameter was widened from <c>PeptideWithSetModifications</c> to
    /// <see cref="IRetentionPredictable"/> so this shim depends only on Omics. Existing call sites
    /// that pass a <c>PeptideWithSetModifications</c> continue to compile unchanged, because that
    /// type implements the interface.
    /// <para>
    /// Deliberately <see cref="IRetentionPredictable"/> rather than <c>IBioPolymerWithSetMods</c>.
    /// SSRCalc3 is a Krokhin reversed-phase model over the 20 proteinogenic residues, and
    /// <c>OligoWithSetMods</c> satisfies <c>IBioPolymerWithSetMods</c> — so the wider signature
    /// would accept an RNA oligo and, because A, C, G and U are also valid amino-acid one-letter
    /// codes, return a plausible-looking number instead of throwing. Nucleic acids deliberately do
    /// not implement <see cref="IRetentionPredictable"/>, so this signature excludes them
    /// structurally rather than by convention.
    /// </para>
    /// </remarks>
    public static double ScoreSequence(this Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3 predictor,
        IRetentionPredictable peptide) => predictor.ScoreSequence(peptide.BaseSequence);
}

/// <inheritdoc cref="ChromatographyExtensions"/>
[Obsolete("Use Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3 directly. This shim will be removed in the next major version.")]
public class SSRCalc3 : Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3
{
    public SSRCalc3(string name, Column column) : base(name, column) { }
}
