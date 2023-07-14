using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Easy.Common.Extensions;
using Test.FileReadingTests;
using UsefulProteomicsDatabases;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;
using Stopwatch = System.Diagnostics.Stopwatch;
using NUnit.Framework.Constraints;
using System.Reflection;

namespace TestFlashLFQ
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal static class TestVerboseOutput
    {

        [Test]
        public static void TestWritePeakDictionary()
        {
            IndexedMassSpectralPeak peak1 = new IndexedMassSpectralPeak(mz: 1 + 0 * Constants.C13MinusC12, 1, -1, -1);
            IndexedMassSpectralPeak peak2 = new IndexedMassSpectralPeak(mz: 1 + 1 * Constants.C13MinusC12, 2, -1, -1);
            IndexedMassSpectralPeak peak3 = new IndexedMassSpectralPeak(mz: 1 + 2 * Constants.C13MinusC12, 3, -1, -1);
            IndexedMassSpectralPeak peak4 = new IndexedMassSpectralPeak(mz: 1 + 3 * Constants.C13MinusC12, 4, -1, -1);
            IndexedMassSpectralPeak peak5 = new IndexedMassSpectralPeak(mz: 1 + 4 * Constants.C13MinusC12, 5, -1, -1);

            // Simplest case, we start with the monoisotopic peak, then have 4 isotope peaks
            var peakDict = VerboseIsotopicEnvelope.WritePeakDictionary(
                peaks: new List<IndexedMassSpectralPeak> { peak1, peak2, peak3, peak4, peak5 },
                monoisotopicMz: 1,
                isotopePpmTolerance: 5);

            Assert.AreEqual(new double[] { 1, 2, 3, 4, 5 }, peakDict.Values.Select(p => p.Intensity).ToArray());
            Assert.AreEqual(new int[] { 0, 1, 2, 3, 4 }, peakDict.Keys.ToArray());

            // Missing monoisotopic peak
            peakDict = VerboseIsotopicEnvelope.WritePeakDictionary(
                peaks: new List<IndexedMassSpectralPeak> { peak2, peak3, peak4, peak5 },
                monoisotopicMz: 1,
                isotopePpmTolerance: 5);

            Assert.AreEqual(new double[] { 2, 3, 4, 5 }, peakDict.Values.Select(p => p.Intensity).ToArray());
            Assert.AreEqual(new int[] { 1, 2, 3, 4 }, peakDict.Keys.ToArray());

            // Only one peak
            peakDict = VerboseIsotopicEnvelope.WritePeakDictionary(
                peaks: new List<IndexedMassSpectralPeak> { peak4 },
                monoisotopicMz: 1,
                isotopePpmTolerance: 5);

            Assert.AreEqual(new double[] { 4 }, peakDict.Values.Select(p => p.Intensity).ToArray());
            Assert.AreEqual(new int[] { 3 }, peakDict.Keys.ToArray());


            // Test recursive tolerance expansion
            double monoisotopicMz = 1 - 1.5e-5; // I believe this is a 15 ppm error

            // Quick check to make sure I understand how tolerance calculations work...
            PpmTolerance tolTest = new PpmTolerance(5);
            Assert.False(tolTest.Within(monoisotopicMz, 1));
            tolTest = new PpmTolerance(15);
            Assert.True(tolTest.Within(monoisotopicMz, 1));

            // All peaks with -15 ppm mass error
            peak1 = new IndexedMassSpectralPeak(mz: monoisotopicMz + 0 * Constants.C13MinusC12, 1, -1, -1);
            peak2 = new IndexedMassSpectralPeak(mz: monoisotopicMz + 1 * Constants.C13MinusC12, 2, -1, -1);
            peak3 = new IndexedMassSpectralPeak(mz: monoisotopicMz + 2 * Constants.C13MinusC12, 3, -1, -1);
            peak4 = new IndexedMassSpectralPeak(mz: monoisotopicMz + 3 * Constants.C13MinusC12, 4, -1, -1);
            peak5 = new IndexedMassSpectralPeak(mz: monoisotopicMz + 4 * Constants.C13MinusC12, 5, -1, -1);

            peakDict = VerboseIsotopicEnvelope.WritePeakDictionary(
                peaks: new List<IndexedMassSpectralPeak> { peak1, peak2, peak3, peak4, peak5 },
                monoisotopicMz: 1,
                isotopePpmTolerance: 5);

            Assert.AreEqual(new double[] { 1, 2, 3, 4, 5 }, peakDict.Values.Select(p => p.Intensity).ToArray());
            Assert.AreEqual(new int[] { 0, 1, 2, 3, 4 }, peakDict.Keys.ToArray());

            // All peaks with +15 ppm mass error
            monoisotopicMz = 1 + 1.5e-5;
            peak1 = new IndexedMassSpectralPeak(mz: monoisotopicMz + 0 * Constants.C13MinusC12, 1, -1, -1);
            peak2 = new IndexedMassSpectralPeak(mz: monoisotopicMz + 1 * Constants.C13MinusC12, 2, -1, -1);
            peak3 = new IndexedMassSpectralPeak(mz: monoisotopicMz + 2 * Constants.C13MinusC12, 3, -1, -1);
            peak4 = new IndexedMassSpectralPeak(mz: monoisotopicMz + 3 * Constants.C13MinusC12, 4, -1, -1);
            peak5 = new IndexedMassSpectralPeak(mz: monoisotopicMz + 4 * Constants.C13MinusC12, 5, -1, -1);

            peakDict = VerboseIsotopicEnvelope.WritePeakDictionary(
                peaks: new List<IndexedMassSpectralPeak> { peak1, peak2, peak3, peak4, peak5 },
                monoisotopicMz: 1,
                isotopePpmTolerance: 5);

            Assert.AreEqual(new double[] { 1, 2, 3, 4, 5 }, peakDict.Values.Select(p => p.Intensity).ToArray());
            Assert.AreEqual(new int[] { 0, 1, 2, 3, 4 }, peakDict.Keys.ToArray());

            // Missing monoisotopic peak with +15 ppm error
            peakDict = VerboseIsotopicEnvelope.WritePeakDictionary(
                peaks: new List<IndexedMassSpectralPeak> { peak2, peak3, peak4, peak5 },
                monoisotopicMz: 1,
                isotopePpmTolerance: 5);

            Assert.AreEqual(new double[] { 2, 3, 4, 5 }, peakDict.Values.Select(p => p.Intensity).ToArray());
            Assert.AreEqual(new int[] { 1, 2, 3, 4 }, peakDict.Keys.ToArray());
        }

        [Test]
        public static void TestChromPeakVerboseOutput()
        {
            IndexedMassSpectralPeak peak1 = new IndexedMassSpectralPeak(mz: 1 + 0 * Constants.C13MinusC12, 1, -1, 1);
            IndexedMassSpectralPeak peak2 = new IndexedMassSpectralPeak(mz: 1 + 1 * Constants.C13MinusC12, 2, -1, 1);
            IndexedMassSpectralPeak peak3 = new IndexedMassSpectralPeak(mz: 1 + 2 * Constants.C13MinusC12, 3, -1, 1);
            IndexedMassSpectralPeak peak4 = new IndexedMassSpectralPeak(mz: 1 + 3 * Constants.C13MinusC12, 4, -1, 1);
            IndexedMassSpectralPeak peak5 = new IndexedMassSpectralPeak(mz: 1 + 4 * Constants.C13MinusC12, 5, -1, 1);

            // Simplest case, we start with the monoisotopic peak, then have 4 isotope peaks
            var envelope1 = new VerboseIsotopicEnvelope(
                peak1,
                allPeaks: new List<IndexedMassSpectralPeak> { peak1, peak2, peak3, peak4, peak5 },
                monoisotopicMass: 1,
                chargeState: 1,
                isotopePpmTolerance: 5);

            SpectraFileInfo fileInfo = new("fake", "", 1, 1, 1);
            ChromatographicPeak chromPeak = new ChromatographicPeak(
                new Identification(fileInfo, "a", "a", 1, 1, 1, new List<ProteinGroup> { new ProteinGroup("a", "a", "a") } ),
                false, fileInfo);
            chromPeak.IsotopicEnvelopes = new List<FlashLFQ.IsotopicEnvelope> { envelope1 };

            string[] toStringArray = chromPeak.ToString(verbose: true).Split('\t');
            var intensityCell = toStringArray[22];
            StringAssert.Contains("i0z1: 1\r\ni1z1: 2", intensityCell);

            // Add a second envelope with a different charge state
            var envelope2 = new VerboseIsotopicEnvelope(
                peak1,
                allPeaks: new List<IndexedMassSpectralPeak> { peak1, peak2, peak3, peak4, peak5 },
                monoisotopicMass: 1,
                chargeState: 2,
                isotopePpmTolerance: 5);
            chromPeak.IsotopicEnvelopes = new List<FlashLFQ.IsotopicEnvelope> { envelope1, envelope2 };
            intensityCell = chromPeak.ToString(verbose: true).Split('\t')[22];
            StringAssert.Contains("i4z1: 5\r\ni0z2: 1\r\ni1z2: 2", intensityCell);

            // Add a third envelope with charge of 2 and retention time of 2
            IndexedMassSpectralPeak peak1time2 = new IndexedMassSpectralPeak(mz: 1 + 0 * Constants.C13MinusC12, 1, -1, 2);
            IndexedMassSpectralPeak peak2time2 = new IndexedMassSpectralPeak(mz: 1 + 1 * Constants.C13MinusC12, 2, -1, 2);
            IndexedMassSpectralPeak peak3time2 = new IndexedMassSpectralPeak(mz: 1 + 2 * Constants.C13MinusC12, 3, -1, 2);
            IndexedMassSpectralPeak peak4time2 = new IndexedMassSpectralPeak(mz: 1 + 3 * Constants.C13MinusC12, 4, -1, 2);
            IndexedMassSpectralPeak peak5time2 = new IndexedMassSpectralPeak(mz: 1 + 4 * Constants.C13MinusC12, 5, -1, 2);

            var envelope2time2 = new VerboseIsotopicEnvelope(
                peak1time2,
                allPeaks: new List<IndexedMassSpectralPeak> { peak1time2, peak2time2, peak3time2, peak4time2, peak5time2 },
                monoisotopicMass: 1,
                chargeState: 2,
                isotopePpmTolerance: 5);

            chromPeak.IsotopicEnvelopes = new List<FlashLFQ.IsotopicEnvelope> { envelope1, envelope2, envelope2time2 };
            intensityCell = chromPeak.ToString(verbose: true).Split('\t')[22];
            var timeCell = chromPeak.ToString(verbose: true).Split('\t')[24];
            StringAssert.Contains("i4z1: 5, -\r\ni0z2: 1, 1\r\ni1z2: 2, 2", intensityCell);
            StringAssert.Contains("1, 2", timeCell);
        }

        public static void SetChromatographicPeakProperties(this ChromatographicPeak peak, string propName, Object newValue)
        {
            PropertyInfo propertyInfo = typeof(ChromatographicPeak).GetProperty(propName);
            if (propertyInfo == null || propertyInfo.PropertyType != newValue.GetType()) return;
            propertyInfo.SetValue(peak, newValue);
        }
    }
}
