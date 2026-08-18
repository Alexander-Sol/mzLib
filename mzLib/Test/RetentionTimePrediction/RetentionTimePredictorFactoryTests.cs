using System;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.Chronologer;
using Chromatography.RetentionTimePrediction.CZE;
using Chromatography.RetentionTimePrediction.SSRCalc;
using MzLibUtil;
using NUnit.Framework;
using Omics.RetentionTimePrediction;

namespace Test.RetentionTimePrediction
{
    /// <summary>
    /// Covers the seam between Chromatography and the opt-in Chromatography.Chronologer assembly.
    /// The factory lives in Chromatography so SSRCalc3 and CZE stay constructible without libtorch,
    /// and reaches Chronologer either through a registered delegate or by resolving the type by
    /// name — neither of which the compiler can check, hence these tests.
    /// </summary>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public class RetentionTimePredictorFactoryTests
    {
        [TearDown]
        public void ResetFactoryState()
        {
            RetentionTimePredictorFactory.RegisterChronologerFactory(null);
            RetentionTimePredictorFactory.ChronologerTypeResolver =
                name => Type.GetType(name, throwOnError: false);
        }

        [Test]
        public void Create_SSRCalc3_ReturnsSsrCalcPredictor()
        {
            using IRetentionTimePredictor predictor = RetentionTimePredictorFactory.Create(PredictorType.SSRCalc3);
            Assert.That(predictor, Is.InstanceOf<SSRCalc3RetentionTimePredictor>());
        }

        [Test]
        public void Create_Cze_ReturnsCzePredictor()
        {
            using IRetentionTimePredictor predictor = RetentionTimePredictorFactory.Create(PredictorType.CZE);
            Assert.That(predictor, Is.InstanceOf<CZERetentionTimePredictor>());
        }

        /// <summary>
        /// The assembly-qualified name is a string, so a rename or a folder move on the Chronologer
        /// side would otherwise surface as a runtime failure in a consumer rather than a red build.
        /// </summary>
        [Test]
        public void ChronologerTypeName_ResolvesToTheChronologerPredictor()
        {
            Type? resolved = Type.GetType(
                RetentionTimePredictorFactory.ChronologerPredictorTypeName, throwOnError: false);

            Assert.That(resolved, Is.Not.Null,
                "Chromatography.Chronologer's predictor could not be resolved by name. If it was " +
                "renamed or moved, update RetentionTimePredictorFactory.ChronologerPredictorTypeName.");
            Assert.That(typeof(IRetentionTimePredictor).IsAssignableFrom(resolved), Is.True);
        }

        /// <summary>
        /// The predictor's only constructor takes optional parameters, which do NOT synthesise a
        /// parameterless constructor — so the factory must activate it with OptionalParamBinding.
        /// Adding a required parameter would break this without breaking compilation.
        /// </summary>
        [Test]
        public void Create_Chronologer_ResolvesFromTheSeparateAssembly()
        {
            using IRetentionTimePredictor predictor = RetentionTimePredictorFactory.Create(PredictorType.Chronologer);

            Assert.That(predictor, Is.InstanceOf<ChronologerRetentionTimePredictor>());
            Assert.That(predictor.PredictorName, Is.EqualTo("Chronologer"));
        }

        [Test]
        public void Create_Chronologer_PrefersARegisteredFactory()
        {
            using var expected = new CZERetentionTimePredictor();
            RetentionTimePredictorFactory.RegisterChronologerFactory(() => expected);

            IRetentionTimePredictor predictor = RetentionTimePredictorFactory.Create(PredictorType.Chronologer);

            Assert.That(predictor, Is.SameAs(expected));
        }

        [Test]
        public void Create_Chronologer_WithoutThePackage_ThrowsAnActionableMzLibException()
        {
            RetentionTimePredictorFactory.ChronologerTypeResolver = _ => null;

            var exception = Assert.Throws<MzLibException>(
                () => RetentionTimePredictorFactory.Create(PredictorType.Chronologer));

            Assert.That(exception!.Message, Does.Contain("mzLib.Chronologer"));
            Assert.That(exception.Message, Does.Contain("SSRCalc3"));
        }

        [Test]
        public void Create_UndefinedPredictorType_Throws()
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => RetentionTimePredictorFactory.Create((PredictorType)(-1)));
        }

        /// <summary>
        /// The compat shim takes <see cref="IRetentionPredictable"/> rather than
        /// <c>IBioPolymerWithSetMods</c> precisely so a nucleic acid cannot reach SSRCalc3's
        /// peptide-only model. OligoWithSetMods implements the latter but not the former, so this
        /// exclusion is structural — if it were ever made to implement IRetentionPredictable, this
        /// assertion is where that shows up.
        /// </summary>
        [Test]
        public void OligoWithSetMods_IsNotRetentionPredictable()
        {
            Assert.That(
                typeof(IRetentionPredictable).IsAssignableFrom(typeof(global::Transcriptomics.Digestion.OligoWithSetMods)),
                Is.False);
        }
    }
}
