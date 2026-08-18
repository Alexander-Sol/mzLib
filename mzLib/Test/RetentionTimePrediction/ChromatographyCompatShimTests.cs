using NUnit.Framework;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Proteomics.RetentionTimePrediction;

// The types under test are the deliberately-obsolete back-compat shims, so calling them is the point.
#pragma warning disable CS0618

namespace Test.RetentionTimePrediction
{
    /// <summary>
    /// Covers the <c>Proteomics.RetentionTimePrediction</c> shims that preserve the type layout
    /// predating the Chromatography namespace. Nothing inside mzLib routes through them — they
    /// exist for MetaMorpheus, whose PlotModelStat.cs constructs
    /// <c>Proteomics.RetentionTimePrediction.SSRCalc3</c> and calls <c>ScoreSequence</c> on a
    /// PeptideWithSetModifications — so without these tests the shims ship uncovered and a change
    /// to them is caught only downstream.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class ChromatographyCompatShimTests
    {
        /// <summary>
        /// The shim class must remain constructible under its old name and delegate to the
        /// relocated SSRCalc3, matching MetaMorpheus's three call sites exactly.
        /// </summary>
        [Test]
        public void ObsoleteSSRCalc3Shim_ScoresTheSameAsTheRelocatedType()
        {
            var shim = new global::Proteomics.RetentionTimePrediction.SSRCalc3(
                "SSRCalc 3.0 (300A)", Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3.Column.A300);
            var relocated = new Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3(
                "SSRCalc 3.0 (300A)", Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3.Column.A300);

            const string sequence = "PEPTIDEK";

            Assert.That(shim.ScoreSequence(sequence), Is.EqualTo(relocated.ScoreSequence(sequence)));
        }

        /// <summary>
        /// The extension overload takes IRetentionPredictable rather than IBioPolymerWithSetMods so
        /// a nucleic acid cannot reach this peptide-only model. A PeptideWithSetModifications
        /// implements it, so the pre-existing call shape still binds and still scores the base
        /// sequence.
        /// </summary>
        [Test]
        public void ObsoleteScoreSequenceExtension_ScoresAPeptidesBaseSequence()
        {
            var predictor = new Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3(
                "SSRCalc 3.0 (300A)", Chromatography.RetentionTimePrediction.SSRCalc.SSRCalc3.Column.A300);
            var peptide = new PeptideWithSetModifications("PEPTIDEK", new System.Collections.Generic.Dictionary<string, Modification>());

            double viaExtension = predictor.ScoreSequence(peptide);

            Assert.That(viaExtension, Is.EqualTo(predictor.ScoreSequence(peptide.BaseSequence)));
        }
    }
}
