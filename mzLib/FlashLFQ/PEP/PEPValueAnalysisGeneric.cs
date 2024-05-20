using Microsoft.ML;
using Microsoft.ML.Data;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Omics;

namespace FlashLFQ.PEP
{
    public static class PEP_Analysis_Cross_Validation
    {
        public static double PipQValueCutoff;

        private static readonly double AbsoluteProbabilityThatDistinguishesPeptides = 0.05;
       
        public static string ComputePEPValuesForAllPeaks(List<ChromatographicPeak> peaks, string outputFolder, int maxThreads, double pepTrainingFraction)
        {
            string[] trainingVariables = ChromatographicPeakData.trainingInfos["standard"];

            //ensure that the order is always stable.
            peaks = peaks.OrderBy(p => p.GetHashCode()).ToList();

            var peakQValues = peaks.Select(peak => peak.MbrQValue).OrderBy(q => q).ToList();
            PipQValueCutoff = peakQValues[(int)Math.Floor(peakQValues.Count * pepTrainingFraction)]; //Select the top 75 percent of all peaks, only use those as positive examples

            //These two dictionaries contain the average and standard deviations of hydrophobicitys measured in 1 minute increments accross each raw
            //file separately. An individully measured hydrobophicty calculated for a specific PSM sequence is compared to these values by computing
            //the z-score. That z-score is used as a feature for machine learning.
            //Separate dictionaries are created for peptides with modifications because SSRcalc doesn't really do a good job predicting hyrophobicity

            //The first string in the dictionary is the filename
            //The value of the dictionary is another dictionary that profiles the hydrophobicity behavior.
            //Each key is a retention time rounded to the nearest minute.
            //The value Tuple is the average and standard deviation, respectively, of the predicted hydrophobicities of the observed peptides eluting at that rounded retention time.

            MLContext mlContext = new MLContext();
            //the number of groups used for cross-validation is hard-coded at four. Do not change this number without changes other areas of effected code.
            const int numGroups = 4;

            List<int>[] peakGroupIndices = Get_PSM_Group_Indices(peaks, numGroups);
            IEnumerable<ChromatographicPeakData>[] ChromatographicPeakDataGroups = new IEnumerable<ChromatographicPeakData>[numGroups];
            for (int i = 0; i < numGroups; i++)
            {
                ChromatographicPeakDataGroups[i] = CreateChromatographicPeakData(peaks, peakGroupIndices[i], maxThreads);
            }

            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[] trainedModels = new TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>>[numGroups];

            var trainer = mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features", numberOfTrees: 400);
            var pipeline = mlContext.Transforms.Concatenate("Features", trainingVariables)
                .Append(trainer);

            List<CalibratedBinaryClassificationMetrics> allMetrics = new List<CalibratedBinaryClassificationMetrics>();
            int sumOfAllAmbiguousPeptidesResolved = 0;

            bool allSetsContainPositiveAndNegativeTrainingExamples = true;
            int groupNumber = 0;
            while (allSetsContainPositiveAndNegativeTrainingExamples == true && groupNumber < numGroups)
            {
                if (ChromatographicPeakDataGroups[groupNumber].Where(p => p.Label == true).Count() == 0 || ChromatographicPeakDataGroups[groupNumber].Where(p => p.Label == false).Count() == 0)
                {
                    allSetsContainPositiveAndNegativeTrainingExamples = false;
                }
                groupNumber++;
            }

            if (allSetsContainPositiveAndNegativeTrainingExamples)
            {
                for (int groupIndexNumber = 0; groupIndexNumber < numGroups; groupIndexNumber++)
                {
                    List<int> allGroupIndexes = Enumerable.Range(0, numGroups).ToList();
                    allGroupIndexes.RemoveAt(groupIndexNumber);

                    //concat doesn't work in a loop, therefore I had to hard code the concat to group 3 out of 4 lists. if the const int numGroups value is changed, then the concat has to be changed accordingly.
                    IDataView dataView = mlContext.Data.LoadFromEnumerable(ChromatographicPeakDataGroups[allGroupIndexes[0]].Concat(ChromatographicPeakDataGroups[allGroupIndexes[1]].Concat(ChromatographicPeakDataGroups[allGroupIndexes[2]])));
                    trainedModels[groupIndexNumber] = pipeline.Fit(dataView);
                    var myPredictions = trainedModels[groupIndexNumber].Transform(mlContext.Data.LoadFromEnumerable(ChromatographicPeakDataGroups[groupIndexNumber]));
                    CalibratedBinaryClassificationMetrics metrics = mlContext.BinaryClassification.Evaluate(data: myPredictions, labelColumnName: "Label", scoreColumnName: "Score");

                    //Parallel operation of the following code requires the method to be stored and then read, once for each thread
                    //if not output directory is specified, the model cannot be stored, and we must force single-threaded operation
                    if (outputFolder != null)
                    {
                        mlContext.Model.Save(trainedModels[groupIndexNumber], dataView.Schema, Path.Combine(outputFolder, "model.zip"));
                    }

                    int whoKnows = Compute_Peak_PEP(peaks, peakGroupIndices[groupIndexNumber], mlContext, trainedModels[groupIndexNumber], outputFolder, maxThreads);

                    allMetrics.Add(metrics);
                    sumOfAllAmbiguousPeptidesResolved += whoKnows;
                }

                return AggregateMetricsForOutput(allMetrics, sumOfAllAmbiguousPeptidesResolved);
            }
            else
            {
                return "Posterior error probability analysis failed. This can occur for small data sets when some sample groups are missing positive or negative training examples.";
            }
        }

        public static string AggregateMetricsForOutput(List<CalibratedBinaryClassificationMetrics> allMetrics, int sumOfAllAmbiguousPeptidesResolved)
        {
            List<double> accuracy = allMetrics.Select(m => m.Accuracy).ToList();
            List<double> areaUnderRocCurve = allMetrics.Select(m => m.AreaUnderRocCurve).ToList();
            List<double> areaUnderPrecisionRecallCurve = allMetrics.Select(m => m.AreaUnderPrecisionRecallCurve).ToList();
            List<double> F1Score = allMetrics.Select(m => m.F1Score).ToList();
            List<double> logLoss = allMetrics.Select(m => m.LogLoss).ToList();
            List<double> logLossReduction = allMetrics.Select(m => m.LogLossReduction).ToList();
            List<double> positivePrecision = allMetrics.Select(m => m.PositivePrecision).ToList();
            List<double> positiveRecall = allMetrics.Select(m => m.PositiveRecall).ToList();
            List<double> negativePrecision = allMetrics.Select(m => m.NegativePrecision).ToList();
            List<double> negativeRecall = allMetrics.Select(m => m.NegativeRecall).ToList();

            // log-loss can stochastically take on a value of infinity.
            // correspondingly, log-loss reduction can be negative infinity.
            // when this happens for one or more of the metrics, it can lead to uninformative numbers.
            // so, unless they are all infinite, we remove them from the average. If they are all infinite, we report that.

            logLoss.RemoveAll(x => x == Double.PositiveInfinity);
            logLossReduction.RemoveAll(x => x == Double.NegativeInfinity);

            double logLossAverage = Double.PositiveInfinity;
            double logLossReductionAverage = Double.NegativeInfinity;

            if ((logLoss != null) && (logLoss.Any()))
            {
                logLossAverage = logLoss.Average();
            }

            if ((logLossReduction != null) && (logLossReduction.Any()))
            {
                logLossReductionAverage = logLossReduction.Average();
            }

            StringBuilder s = new StringBuilder();
            s.AppendLine();
            s.AppendLine("************************************************************");
            s.AppendLine("*       Metrics for Determination of PEP Using Binary Classification      ");
            s.AppendLine("*-----------------------------------------------------------");
            s.AppendLine("*       Accuracy:  " + accuracy.Average().ToString());
            s.AppendLine("*       Area Under Curve:  " + areaUnderRocCurve.Average().ToString());
            s.AppendLine("*       Area under Precision recall Curve:  " + areaUnderPrecisionRecallCurve.Average().ToString());
            s.AppendLine("*       F1Score:  " + F1Score.Average().ToString());
            s.AppendLine("*       LogLoss:  " + logLossAverage.ToString());
            s.AppendLine("*       LogLossReduction:  " + logLossReductionAverage.ToString());
            s.AppendLine("*       PositivePrecision:  " + positivePrecision.Average().ToString());
            s.AppendLine("*       PositiveRecall:  " + positiveRecall.Average().ToString());
            s.AppendLine("*       NegativePrecision:  " + negativePrecision.Average().ToString());
            s.AppendLine("*       NegativeRecall:  " + negativeRecall.Average().ToString());
            s.AppendLine("*       Count of Ambiguous Peptides Removed:  " + sumOfAllAmbiguousPeptidesResolved.ToString());
            s.AppendLine("************************************************************");
            return s.ToString();
        }

        public static int Compute_Peak_PEP(
            List<ChromatographicPeak> peaks,
            List<int> peakIndicies,
            MLContext mLContext,
            TransformerChain<BinaryPredictionTransformer<Microsoft.ML.Calibrators.CalibratedModelParametersBase<Microsoft.ML.Trainers.FastTree.FastTreeBinaryModelParameters, Microsoft.ML.Calibrators.PlattCalibrator>>> trainedModel,
            string outputFolder, int maxThreads)
        {
            object lockObject = new object();

            //the trained model is not threadsafe. Therefore, to use the same model for each thread saved the model to disk. Then each thread reads its own copy of the model back from disk.
            //If there is no output folder specified, then this can't happen. We set maxthreads eqaul to one and use the model that gets passed into the method.
            if (String.IsNullOrEmpty(outputFolder))
            {
                maxThreads = 1;
            }

            Parallel.ForEach(Partitioner.Create(0, peakIndicies.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {

                    ITransformer threadSpecificTrainedModel;
                    if (maxThreads == 1)
                    {
                        threadSpecificTrainedModel = trainedModel;
                    }
                    else
                    {
                        threadSpecificTrainedModel = mLContext.Model.Load(Path.Combine(outputFolder, "model.zip"), out DataViewSchema savedModelSchema);
                    }

                    // one prediction engine per thread, because the prediction engine is not thread-safe
                    var threadPredictionEngine = mLContext.Model.CreatePredictionEngine<ChromatographicPeakData, TruePositivePrediction>(threadSpecificTrainedModel);

                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        ChromatographicPeak peak = peaks[peakIndicies[i]];

                        if (peak != null)
                        {
                            List<int> indiciesOfPeptidesToRemove = new List<int>();
                            List<double> pepValuePredictions = new List<double>();

                            //Here we compute the pepvalue predection for each ambiguous peptide in a PSM. Ambiguous peptides with lower pepvalue predictions are removed from the PSM.

                            List<int> allBmpNotches = new List<int>();
                            List<IBioPolymerWithSetMods> allBmpPeptides = new List<IBioPolymerWithSetMods>();


                            ChromatographicPeakData pd = CreateOneChromatographicPeakDataEntry(peak, label: !peak.RandomRt) ;
                            var pepValuePrediction = threadPredictionEngine.Predict(pd);
                            pepValuePredictions.Add(pepValuePrediction.Probability);
                            peak.PipPep = 1 - pepValuePrediction.Probability;

                            //A score is available using the variable pepvaluePrediction.Score
                        }
                    }
                });
            return 1;
        }

        //we add the indexes of the targets and decoys to the groups separately in the hope that we'll get at least one target and one decoy in each group.
        //then training can possibly be more successful.
        public static List<int>[] Get_PSM_Group_Indices(List<ChromatographicPeak> peaks, int numGroups)
        {
            List<int>[] groupsOfIndicies = new List<int>[numGroups];
            for (int i = 0; i < numGroups; i++)
            {
                groupsOfIndicies[i] = new List<int>();
            }

            List<int> targetPeakIndexes = new List<int>();
            List<int> decoyPeakIndexes = new List<int>();

            for (int i = 0; i < peaks.Count; i++)
            {
                if (peaks[i].RandomRt)
                {
                    decoyPeakIndexes.Add(i);
                }
                else
                {
                    targetPeakIndexes.Add(i);
                }
            }

            int myIndex = 0;

            while (myIndex < decoyPeakIndexes.Count)
            {
                int subIndex = 0;
                while (subIndex < numGroups && myIndex < decoyPeakIndexes.Count)
                {
                    groupsOfIndicies[subIndex].Add(decoyPeakIndexes[myIndex]);
                    subIndex++;
                    myIndex++;
                }
            }

            myIndex = 0;

            while (myIndex < targetPeakIndexes.Count)
            {
                int subIndex = 0;
                while (subIndex < numGroups && myIndex < targetPeakIndexes.Count)
                {
                    groupsOfIndicies[subIndex].Add(targetPeakIndexes[myIndex]);
                    subIndex++;
                    myIndex++;
                }
            }

            return groupsOfIndicies;
        }


        public static IEnumerable<ChromatographicPeakData> CreateChromatographicPeakData(List<ChromatographicPeak> peaks, List<int> peakIndicies, int maxThreads)
        {
            object ChromatographicPeakDataListLock = new object();
            List<ChromatographicPeakData> ChromatographicPeakDataList = new List<ChromatographicPeakData>();
            int[] threads = Enumerable.Range(0, maxThreads).ToArray();

            Parallel.ForEach(Partitioner.Create(0, peakIndicies.Count),
                new ParallelOptions { MaxDegreeOfParallelism = maxThreads },
                (range, loopState) =>
                {
                    List<ChromatographicPeakData> localChromatographicPeakDataList = new List<ChromatographicPeakData>();
                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        var peak = peaks[peakIndicies[i]];
                        ChromatographicPeakData newChromatographicPeakData = new ChromatographicPeakData();

                        bool label;

                        if (peak.RandomRt && !peak.DecoyPeptide)
                        {
                            label = false;
                            newChromatographicPeakData = CreateOneChromatographicPeakDataEntry(peak, label);
                        }
                        else if (!peak.RandomRt && !peak.DecoyPeptide && peak.MbrQValue <= PipQValueCutoff)
                        {
                            label = true;
                            newChromatographicPeakData = CreateOneChromatographicPeakDataEntry(peak, label);
                        }
                        localChromatographicPeakDataList.Add(newChromatographicPeakData);
                    }
                    lock (ChromatographicPeakDataListLock)
                    {
                        ChromatographicPeakDataList.AddRange(localChromatographicPeakDataList);
                    }
                });

            ChromatographicPeakData[] pda = ChromatographicPeakDataList.ToArray();

            return pda.AsEnumerable();
        }

        public static ChromatographicPeakData CreateOneChromatographicPeakDataEntry(ChromatographicPeak peak,bool label)
        {

            peak.PepPeakData = new ChromatographicPeakData
            {
                PpmErrorScore = (float)peak.PpmScore,
                IntensityScore = (float)peak.IntensityScore,
                RtScore = (float)peak.RtScore,
                ScanCountScore = (float)peak.ScanCountScore,

                PpmErrorRaw = (float)Math.Abs(peak.MassError),
                IntensityRaw = (float)Math.Log2(peak.Intensity),
                RtPredictionErrorRaw = (float)Math.Abs(peak.RtPredictionError),
                ScanCountRaw = (float)peak.IsotopicEnvelopes.Count,
                IsotopicEnvelopeCorrelation = (float)(peak.Apex.PearsonCorrelation ?? 0),

                Label = label,

            };

            return peak.PepPeakData;
        }
    }
}