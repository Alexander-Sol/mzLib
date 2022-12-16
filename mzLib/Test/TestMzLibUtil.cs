﻿using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MzLibUtil;
using MassSpectrometry;
using System.Globalization;
using System.IO;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public sealed class TestMzLibUtil
    {
        [Test]
        [TestCase(@"C:\Users\bubba\Documents\Projects\K562\K562_2\20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase("C:\\Users\bubba\\Documents\\Projects\\K562\\K562_2\\20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase("20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase("20100730_Velos1_TaGe_SA_K565_4", "20100730_Velos1_TaGe_SA_K565_4")]
        [TestCase(@"C:\Users\bubba\Documents\Projects\K562\K562_2\20100730_Velos1_TaGe_SA_K565_4", "20100730_Velos1_TaGe_SA_K565_4")]
        //test extra period in folder name of path
        [TestCase(@"C:\Users\bubba\Documents.docs\Projects\K562\K562_2\20100730_Velos1_TaGe_SA_K565_4.raw", "20100730_Velos1_TaGe_SA_K565_4")]
        //test extra period in filename
        [TestCase(@"C:\Users\bubba\Documents\Projects\K562\K562_2\20100730_Velos1_.TaGe_SA_K565_4.raw", "20100730_Velos1_.TaGe_SA_K565_4")]
        [TestCase("/home/seth/Pictures/penguin.jpg","penguin")]
        [TestCase("/home/seth/Pictures/penguin", "penguin")]
        [TestCase("penguin.jpg", "penguin")]
        [TestCase("penguin", "penguin")]
        [TestCase("penguin.jpg.gz", "penguin")]
        [TestCase("penguin.jpg.zip", "penguin")]
        public static void TestPeriodTolerantFilenameWithoutExtension(string filenameAndOrPath, string expectedResult)
        {
            string result = PeriodTolerantFilenameWithoutExtension.GetPeriodTolerantFilenameWithoutExtension(filenameAndOrPath);
            Assert.AreEqual(expectedResult, result);
        }

        [Test]
        public void TestSpectrumTreeInsertion()
        {
            MzSpectrum unorderedSpectrum = new(new double[] { 5, 6, 1, 2, 3, 4, 7, 8 }, new double[] { 10, 8, 2, 4, 6, 8, 6, 4 }, false);
            MzSpectrum orderedSpectrum = new(new double[] { 1, 2, 3, 4, 5, 6, 7, 8 }, new double[] { 2, 4, 6, 8, 10, 8, 6, 4 }, false);
            SpectrumTree testTree = new();
            for (int i = 0; i < 8; i++)
            {
                testTree.Insert(new Node(unorderedSpectrum.XArray[i], unorderedSpectrum.YArray[i]));
            }

            double[] mz = new double[8];
            double[] intensity = new double[8];
            int j = 0;
            foreach (Node node in testTree)
            {
                mz[j] = node.Key;
                intensity[j] = node.Value;
                j++;
            }

            MzSpectrum reconstructedSpectrum = new(mz, intensity, false);
            Assert.That(orderedSpectrum.Equals(reconstructedSpectrum));
        }

        [Test]
        public void TestSpectrumTreeDeletion()
        {
            MzSpectrum unorderedSpectrum = new(new double[] {  5, 6, 1, 2, 3, 4, 7, 8 }, new double[] { 10, 8, 2, 4, 6, 8, 6, 4 }, false);
            SpectrumTree testTree = new();
            List<Node> nodeList = new List<Node>();
            for (int i = 0; i < 8; i++)
            {
                nodeList.Add(new Node(unorderedSpectrum.XArray[i], unorderedSpectrum.YArray[i]));
                testTree.Insert(nodeList[i]);
            }

            testTree.Delete(nodeList[3]);
            testTree.Delete(nodeList[5]);
            double[] mz = new double[6];
            double[] intensity = new double[6];
            int j = 0;
            foreach (Node node in testTree)
            {
                mz[j] = node.Key;
                intensity[j] = node.Value;
                j++;
            }

            MzSpectrum orderedSpectrum = new(new double[] { 1, 3, 5, 6, 7, 8 }, new double[] { 2, 6, 10, 8, 6, 4 }, false);
            MzSpectrum reconstructedSpectrum = new(mz, intensity, false);
            Assert.That(orderedSpectrum.Equals(reconstructedSpectrum));
        }

        [Test]
        public void TestSpectrumTreeBuilder()
        {
            MzSpectrum orderedSpectrum = new(new double[] { 1, 2, 3, 4, 5, 6, 7, 8 }, new double[] { 2, 4, 6, 8, 10, 8, 6, 4 }, false);
            SpectrumTree testTree = new();
            testTree.BuildTreeFromSpectrum(orderedSpectrum.XArray, orderedSpectrum.YArray);

            double[] mz = new double[8];
            double[] intensity = new double[8];
            int j = 0;
            foreach (Node node in testTree)
            {
                mz[j] = node.Key;
                intensity[j] = node.Value;
                j++;
            }

            MzSpectrum reconstructedSpectrum = new(mz, intensity, false);
            Assert.That(orderedSpectrum.Equals(reconstructedSpectrum));

            // Test with actual data
            //txt file, not mgf, because it's an MS1. Most intense proteoform has mass of ~14037.9 Da
            string Ms1SpectrumPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DataFiles\14kDaProteoformMzIntensityMs1.txt");

            string[] spectrumLines = File.ReadAllLines(Ms1SpectrumPath);

            int mzIntensityPairsCount = spectrumLines.Length;
            double[] ms1mzs = new double[mzIntensityPairsCount];
            double[] ms1intensities = new double[mzIntensityPairsCount];

            for (int i = 0; i < mzIntensityPairsCount; i++)
            {
                string[] pair = spectrumLines[i].Split('\t');
                ms1mzs[i] = Convert.ToDouble(pair[0], CultureInfo.InvariantCulture);
                ms1intensities[i] = Convert.ToDouble(pair[1], CultureInfo.InvariantCulture);
            }

            MzSpectrum actualSpectrum = new MzSpectrum(ms1mzs, ms1intensities, false);
            testTree = new();
            testTree.BuildTreeFromSpectrum(actualSpectrum.XArray, actualSpectrum.YArray);

            mz = new double[mzIntensityPairsCount];
            intensity = new double[mzIntensityPairsCount];
            j = 0;
            foreach (Node node in testTree)
            {
                mz[j] = node.Key;
                intensity[j] = node.Value;
                j++;
            }

            reconstructedSpectrum = new(mz, intensity, false);
            Assert.That(actualSpectrum.Equals(reconstructedSpectrum));

        }
    }
}
