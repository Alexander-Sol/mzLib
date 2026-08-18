namespace Omics.RetentionTimePrediction;

/// <summary>
/// Represents an entity for which retention time can be predicted.
/// Designed for minimal allocation and no dependencies on higher layers.
/// </summary>
/// <remarks>
/// This interface lives in Omics rather than Chromatography deliberately. Every real
/// implementer is an Omics type (see <see cref="IBioPolymerWithSetMods"/>), so hosting it
/// here lets biopolymer projects describe themselves as retention-predictable without
/// referencing Chromatography — which in turn keeps the deep-learning predictors and their
/// native runtime out of the dependency graph of every project that merely handles peptides.
/// <para>
/// <b>Migration:</b> this interface was declared as
/// <c>Chromatography.RetentionTimePrediction.IRetentionPredictable</c> from its introduction in
/// December 2025 until the Chronologer split, so the move is a deliberate source break for
/// external implementers — the one relocation in that change that did not preserve its namespace.
/// Replace <c>using Chromatography.RetentionTimePrediction;</c> with
/// <c>using Omics.RetentionTimePrediction;</c> where this type is named.
/// </para>
/// <para>
/// No <c>[Obsolete]</c> shim is provided in the old namespace, and deliberately so: an empty
/// derived interface of the same name would sit in the namespace that <c>Chromatography</c>'s own
/// predictors declare, so every consumer holding both usings — <c>PredictionClients</c> among them
/// — would get CS0104 ambiguity instead of the CS0246 it was meant to prevent. The shim would
/// break the callers it was written to protect.
/// </para>
/// </remarks>
public interface IRetentionPredictable
{
    /// <summary>
    /// Gets the base (unmodified) sequence
    /// </summary>
    string BaseSequence { get; }

    /// <summary>
    /// Gets the full sequence representation with modification identifiers
    /// e.g., "PEPTIDE[Variable:Oxidation on M]K[Variable:Acetylation on K]"
    /// </summary>
    string FullSequence { get; }

    /// <summary>
    /// Gets the monoisotopic mass of the peptide.
    /// Required for CZE electrophoretic mobility predictions.
    /// </summary>
    double MonoisotopicMass { get; }

    /// <summary>
    /// Builds a sequence string with mass shifts for modifications.
    /// Format: "PEPTIDE[+15.995]K[+42.011]" or "PEPTIDE[-17.026]K"
    /// This is used by predictors that work with mass-based representations (e.g., Chronologer).
    /// </summary>
    /// <returns>Sequence with mass shift annotations, or null if not applicable</returns>
    string FullSequenceWithMassShifts { get; }
}
