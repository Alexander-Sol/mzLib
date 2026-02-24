using System.Collections.Concurrent;
using System.Threading.Channels;
using MzLibUtil;

namespace Readers
{
    public static class SpectrumMatchTsvReader
    {
        private static readonly char[] Split = { '\t' };

        /// <summary>
        /// File size threshold (in bytes) above which the streaming producer-consumer path is used.
        /// Files below this size use the simpler ReadAllLines + Parallel.For approach.
        /// </summary>
        private const long LargeFileThreshold = 1 * 1024 * 1024; // 1 GB

        private const int BytesPerLineEstimate = 1000; // Conservative estimate of average line length for pre-allocating arrays

        /// <summary>
        /// Reads a TSV file, choosing between a batch approach (small files) and a streaming
        /// producer-consumer approach (large files) to balance throughput and memory usage.
        /// </summary>
        /// <exception cref="MzLibException"></exception>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        public static List<T> ReadTsv<T>(string filePath, out List<string> warnings) where T : SpectrumMatchFromTsv
        {
            MzLibException? parsingException = null;
            SupportedFileType type;
            try
            {
                type = filePath.ParseFileType();
            }
            catch (MzLibException e)
            {
                parsingException = e;
                type = SupportedFileType.psmtsv;
            }

            long fileSize;
            try
            {
                fileSize = new FileInfo(filePath).Length;
            }
            catch (Exception e)
            {
                throw new MzLibException("Could not read file: " + e.Message, e);
            }

            List<T> psms;
            if (fileSize >= LargeFileThreshold)
            {
                psms = ReadTsvStreaming<T>(filePath, type, out warnings);
            }
            else
            {
                psms = ReadTsvBatch<T>(filePath, type, out warnings);
            }

            if (parsingException is not null && psms.Count == 0)
            {
                throw new MzLibException($"No spectral matches found in file: {filePath}", parsingException);
            }

            return psms;
        }

        /// <summary>
        /// Batch approach: reads all lines into memory, then parses in parallel.
        /// Best for small-to-medium files where I/O is fast and memory is not a concern.
        /// </summary>
        private static List<T> ReadTsvBatch<T>(string filePath, SupportedFileType type, out List<string> warnings) where T : SpectrumMatchFromTsv
        {
            string[] lines;
            try
            {
                lines = File.ReadAllLines(filePath);
            }
            catch (Exception e)
            {
                throw new MzLibException("Could not read file: " + e.Message, e);
            }

            Dictionary<string, int> parsedHeader = ParseHeader(lines[0]);
            bool fileIsGlyco = parsedHeader.ContainsKey(SpectrumMatchFromTsvHeader.GlycanMass) && parsedHeader[SpectrumMatchFromTsvHeader.GlycanMass] != -1;
            int lineCount = lines.Length - 1;

            T?[] psmsArray = new T[lineCount];
            var warningsBag = new ConcurrentBag<string>();
            int maxThreads = Math.Max(1, Math.Min(8, Environment.ProcessorCount - 1));
            int chunkSize = (int)Math.Ceiling((double)lineCount / maxThreads);

            Parallel.For(0, maxThreads, new ParallelOptions { MaxDegreeOfParallelism = maxThreads }, threadIdx =>
            {
                int start = 1 + threadIdx * chunkSize;
                int end = Math.Min(lines.Length, start + chunkSize);
                for (int i = start; i < end; i++)
                {
                    try
                    {
                        psmsArray[i - 1] = ParseLine<T>(lines[i], type, fileIsGlyco, parsedHeader);
                    }
                    catch (Exception)
                    {
                        warningsBag.Add("Could not read line: " + (i + 1));
                    }
                }
            });

            var psms = new List<T>(lineCount);
            for (int i = 0; i < psmsArray.Length; i++)
            {
                if (psmsArray[i] != null)
                    psms.Add(psmsArray[i]!);
            }
            warnings = warningsBag.ToList();

            if (lineCount != psms.Count)
            {
                warnings.Add("Warning: " + (lineCount - psms.Count) + " PSMs were not read.");
            }

            return psms;
        }

        /// <summary>
        /// Streaming approach: one producer thread reads lines from disk while multiple consumer
        /// threads parse them in parallel. Overlaps I/O with CPU work and limits memory usage
        /// via a bounded channel. Uses a pre-allocated array indexed by line number to avoid
        /// ConcurrentBag overhead and post-sort.
        /// </summary>
        private static List<T> ReadTsvStreaming<T>(string filePath, SupportedFileType type, out List<string> warnings) where T : SpectrumMatchFromTsv
        {
            // Estimate line count from file size to pre-allocate the results array.
            // Average TSV line is ~200-500 bytes; we use a conservative estimate (100 bytes)
            // so the array is large enough. We'll trim at the end.
            long fileSize;
            try
            {
                fileSize = new FileInfo(filePath).Length;
            }
            catch (Exception e)
            {
                throw new MzLibException("Could not read file: " + e.Message, e);
            }

            int estimatedLines = (int)Math.Min(fileSize / BytesPerLineEstimate, int.MaxValue);

            // Single file open: read header, then produce data lines from the same stream
            StreamReader reader;
            string headerLine;
            try
            {
                reader = new StreamReader(filePath);
                headerLine = reader.ReadLine() ?? throw new MzLibException("File is empty: " + filePath);
            }
            catch (MzLibException) { throw; }
            catch (Exception e)
            {
                throw new MzLibException("Could not read file: " + e.Message, e);
            }

            Dictionary<string, int> parsedHeader = ParseHeader(headerLine);
            bool fileIsGlyco = parsedHeader.ContainsKey(SpectrumMatchFromTsvHeader.GlycanMass) && parsedHeader[SpectrumMatchFromTsvHeader.GlycanMass] != -1;

            var warningsBag = new ConcurrentBag<string>();

            // Pre-allocate results array; will grow if estimate was too low
            T?[] resultsArray = new T[estimatedLines];
            object resizeLock = new object();
            int totalLineCount = 0;

            var channel = Channel.CreateBounded<(int lineNumber, string line)>(
                new BoundedChannelOptions(4096)
                {
                    SingleWriter = true,
                    FullMode = BoundedChannelFullMode.Wait
                });

            // Producer: stream lines from the already-opened reader, using async WriteAsync
            // to yield when the channel is full instead of spinning (avoids thread pool starvation)
            var producerTask = Task.Run(async () =>
            {
                int lineNum = 0;
                try
                {
                    var writer = channel.Writer;
                    string? line;
                    while ((line = reader.ReadLine()) != null)
                    {
                        await writer.WriteAsync((lineNum, line));
                        lineNum++;
                    }
                    Interlocked.Exchange(ref totalLineCount, lineNum);
                }
                finally
                {
                    channel.Writer.Complete();
                    reader.Dispose();
                }
            });

            // Consumers: parse lines and write directly into the pre-allocated array by index
            int maxConsumers = Math.Max(1, Math.Min(8, Environment.ProcessorCount - 1));
            var consumerTasks = new Task[maxConsumers];
            for (int c = 0; c < maxConsumers; c++)
            {
                consumerTasks[c] = Task.Run(async () =>
                {
                    var channelReader = channel.Reader;
                    while (await channelReader.WaitToReadAsync())
                    {
                        while (channelReader.TryRead(out var item))
                        {
                            var (lineNumber, line) = item;

                            // Grow array if our estimate was too low (stable lock object)
                            var arr = Volatile.Read(ref resultsArray);
                            if (lineNumber >= arr.Length)
                            {
                                lock (resizeLock)
                                {
                                    arr = resultsArray;
                                    if (lineNumber >= arr.Length)
                                    {
                                        int newSize = Math.Max(arr.Length * 2, lineNumber + 1);
                                        var newArray = new T?[newSize];
                                        Array.Copy(arr, newArray, arr.Length);
                                        Volatile.Write(ref resultsArray, newArray);
                                    }
                                }
                                arr = Volatile.Read(ref resultsArray);
                            }

                            try
                            {
                                arr[lineNumber] = ParseLine<T>(line, type, fileIsGlyco, parsedHeader);
                            }
                            catch (Exception)
                            {
                                warningsBag.Add("Could not read line: " + (lineNumber + 2)); // +2: 1-based + header
                            }
                        }
                    }
                });
            }
            
            Task.WaitAll(producerTask);
            Task.WaitAll(consumerTasks);

            // Collect non-null results in order (array is already indexed by line number)
            int lineCount = totalLineCount;
            var psms = new List<T>(lineCount);
            for (int i = 0; i < lineCount; i++)
            {
                if (resultsArray[i] != null)
                    psms.Add(resultsArray[i]!);
            }

            warnings = warningsBag.ToList();

            if (lineCount != psms.Count)
            {
                warnings.Add("Warning: " + (lineCount - psms.Count) + " PSMs were not read.");
            }

            return psms;
        }

        /// <summary>
        /// Parses a single TSV line into the appropriate SpectrumMatchFromTsv subtype.
        /// </summary>
        private static T ParseLine<T>(string line, SupportedFileType type, bool fileIsGlyco, Dictionary<string, int> parsedHeader) where T : SpectrumMatchFromTsv
        {
            bool lineIsGlyco = fileIsGlyco && ResultIsGlyco(parsedHeader, line);

            T result = type switch
            {
                SupportedFileType.osmtsv => (T)(SpectrumMatchFromTsv)new OsmFromTsv(line, Split, parsedHeader),
                _ when lineIsGlyco => (T)(SpectrumMatchFromTsv)new GlycoPsmFromTsv(line, Split, parsedHeader),
                _ => (T)(SpectrumMatchFromTsv)new PsmFromTsv(line, Split, parsedHeader)
            };
            return result;
        }

        /// <summary>
        /// Legacy method for reading PsmFromTsv files, creates a generic SpectrumMatchFromTsv object for each line
        /// </summary>
        public static List<SpectrumMatchFromTsv> ReadTsv(string filePath, out List<string> warnings) =>
            ReadTsv<SpectrumMatchFromTsv>(filePath, out warnings);

        /// <summary>
        /// Reads a psmtsv file and returns PsmFromTsv objects
        /// </summary>
        public static List<PsmFromTsv> ReadPsmTsv(string filePath, out List<string> warnings) =>
            ReadTsv<PsmFromTsv>(filePath, out warnings);

        public static List<GlycoPsmFromTsv> ReadGlycoPsmTsv(string filePath, out List<string> warnings) =>
            ReadTsv<GlycoPsmFromTsv>(filePath, out warnings);

        /// <summary>
        /// Reads a osmtsv file and returns OsmFromTsv objects
        /// </summary>
        public static List<OsmFromTsv> ReadOsmTsv(string filePath, out List<string> warnings) =>
            ReadTsv<OsmFromTsv>(filePath, out warnings);

        public static Dictionary<string, int> ParseHeader(string header)
        {
            var parsedHeader = new Dictionary<string, int>();
            var spl = header.Split(Split);

            // Build a reverse lookup: column name -> index (single pass over header)
            var columnIndex = new Dictionary<string, int>(spl.Length);
            for (int i = 0; i < spl.Length; i++)
            {
                columnIndex[spl[i]] = i;
            }

            // Helper to get index or -1
            int IndexOf(string columnName) => columnIndex.TryGetValue(columnName, out int idx) ? idx : -1;

            parsedHeader.Add(SpectrumMatchFromTsvHeader.FullSequence, IndexOf(SpectrumMatchFromTsvHeader.FullSequence));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.Ms2ScanNumber, IndexOf(SpectrumMatchFromTsvHeader.Ms2ScanNumber));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.FileName, IndexOf(SpectrumMatchFromTsvHeader.FileName));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.TotalIonCurrent, IndexOf(SpectrumMatchFromTsvHeader.TotalIonCurrent));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.PrecursorScanNum, IndexOf(SpectrumMatchFromTsvHeader.PrecursorScanNum));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.PrecursorCharge, IndexOf(SpectrumMatchFromTsvHeader.PrecursorCharge));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.PrecursorIntensity, IndexOf(SpectrumMatchFromTsvHeader.PrecursorIntensity));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.PrecursorMz, IndexOf(SpectrumMatchFromTsvHeader.PrecursorMz));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.PrecursorMass, IndexOf(SpectrumMatchFromTsvHeader.PrecursorMass));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.OneOverK0, IndexOf(SpectrumMatchFromTsvHeader.OneOverK0));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.Score, IndexOf(SpectrumMatchFromTsvHeader.Score));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.DeltaScore, IndexOf(SpectrumMatchFromTsvHeader.DeltaScore));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.Notch, IndexOf(SpectrumMatchFromTsvHeader.Notch));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BaseSequence, IndexOf(SpectrumMatchFromTsvHeader.BaseSequence));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.EssentialSequence, IndexOf(SpectrumMatchFromTsvHeader.EssentialSequence));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.AmbiguityLevel, IndexOf(SpectrumMatchFromTsvHeader.AmbiguityLevel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.MissedCleavages, IndexOf(SpectrumMatchFromTsvHeader.MissedCleavages));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.MassDiffDa, IndexOf(SpectrumMatchFromTsvHeader.MassDiffDa));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.MassDiffPpm, IndexOf(SpectrumMatchFromTsvHeader.MassDiffPpm));

            //Handle legacy input
            if (columnIndex.ContainsKey(SpectrumMatchFromTsvHeader.Accession))
            {
                parsedHeader.Add(SpectrumMatchFromTsvHeader.SpectrumMatchCount, IndexOf(SpectrumMatchFromTsvHeader.SpectrumMatchCount));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.MonoisotopicMass, IndexOf(SpectrumMatchFromTsvHeader.MonoisotopicMass));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.Accession, IndexOf(SpectrumMatchFromTsvHeader.Accession));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.Name, IndexOf(SpectrumMatchFromTsvHeader.Name));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.Description, IndexOf(SpectrumMatchFromTsvHeader.Description));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.StartAndEndResiduesInFullSequence, IndexOf(SpectrumMatchFromTsvHeader.StartAndEndResiduesInFullSequence));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.NextResidue, IndexOf(SpectrumMatchFromTsvHeader.NextResidue));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.PreviousResidue, IndexOf(SpectrumMatchFromTsvHeader.PreviousResidue));
            }
            else
            {
                parsedHeader.Add(SpectrumMatchFromTsvHeader.SpectrumMatchCount, IndexOf(SpectrumMatchFromTsvHeader.PsmCount));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.MonoisotopicMass, IndexOf(SpectrumMatchFromTsvHeader.PeptideMonoMass));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.Accession, IndexOf(SpectrumMatchFromTsvHeader.ProteinAccession));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.Name, IndexOf(SpectrumMatchFromTsvHeader.ProteinName));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.Description, IndexOf(SpectrumMatchFromTsvHeader.PeptideDescription));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.StartAndEndResiduesInFullSequence, IndexOf(SpectrumMatchFromTsvHeader.StartAndEndResiduesInProtein));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.NextResidue, IndexOf(SpectrumMatchFromTsvHeader.NextAminoAcid));
                parsedHeader.Add(SpectrumMatchFromTsvHeader.PreviousResidue, IndexOf(SpectrumMatchFromTsvHeader.PreviousAminoAcid));
            }

            parsedHeader.Add(SpectrumMatchFromTsvHeader.FlankingResidues, IndexOf(SpectrumMatchFromTsvHeader.FlankingResidues));
            if (parsedHeader[SpectrumMatchFromTsvHeader.FlankingResidues] == -1) // try legacy name from previous versions
            {
                parsedHeader[SpectrumMatchFromTsvHeader.FlankingResidues] = IndexOf("FlankingResidues");
            }

            parsedHeader.Add(SpectrumMatchFromTsvHeader.NumberOfMods, IndexOf(SpectrumMatchFromTsvHeader.NumberOfMods));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.GeneName, IndexOf(SpectrumMatchFromTsvHeader.GeneName));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.OrganismName, IndexOf(SpectrumMatchFromTsvHeader.OrganismName));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.IntersectingSequenceVariations, IndexOf(SpectrumMatchFromTsvHeader.IntersectingSequenceVariations));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.IdentifiedSequenceVariations, IndexOf(SpectrumMatchFromTsvHeader.IdentifiedSequenceVariations));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.SpliceSites, IndexOf(SpectrumMatchFromTsvHeader.SpliceSites));

            parsedHeader.Add(SpectrumMatchFromTsvHeader.DecoyContaminantTarget, IndexOf(SpectrumMatchFromTsvHeader.DecoyContaminantTarget));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.MatchedIonMzRatios, IndexOf(SpectrumMatchFromTsvHeader.MatchedIonMzRatios));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.MatchedIonIntensities, IndexOf(SpectrumMatchFromTsvHeader.MatchedIonIntensities));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.MatchedIonMassDiffDa, IndexOf(SpectrumMatchFromTsvHeader.MatchedIonMassDiffDa));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.SpectralAngle, IndexOf(SpectrumMatchFromTsvHeader.SpectralAngle));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.QValue, IndexOf(SpectrumMatchFromTsvHeader.QValue));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.QValueNotch, IndexOf(SpectrumMatchFromTsvHeader.QValueNotch));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.PEP, IndexOf(SpectrumMatchFromTsvHeader.PEP));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.PEP_QValue, IndexOf(SpectrumMatchFromTsvHeader.PEP_QValue));

            parsedHeader.Add(SpectrumMatchFromTsvHeader.CrossTypeLabel, IndexOf(SpectrumMatchFromTsvHeader.CrossTypeLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.LinkResiduesLabel, IndexOf(SpectrumMatchFromTsvHeader.LinkResiduesLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.ProteinLinkSiteLabel, IndexOf(SpectrumMatchFromTsvHeader.ProteinLinkSiteLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.RankLabel, IndexOf(SpectrumMatchFromTsvHeader.RankLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideProteinAccessionLabel, IndexOf(SpectrumMatchFromTsvHeader.BetaPeptideProteinAccessionLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideProteinLinkSiteLabel, IndexOf(SpectrumMatchFromTsvHeader.BetaPeptideProteinLinkSiteLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideBaseSequenceLabel, IndexOf(SpectrumMatchFromTsvHeader.BetaPeptideBaseSequenceLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideFullSequenceLabel, IndexOf(SpectrumMatchFromTsvHeader.BetaPeptideFullSequenceLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideTheoreticalMassLabel, IndexOf(SpectrumMatchFromTsvHeader.BetaPeptideTheoreticalMassLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideScoreLabel, IndexOf(SpectrumMatchFromTsvHeader.BetaPeptideScoreLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideRankLabel, IndexOf(SpectrumMatchFromTsvHeader.BetaPeptideRankLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel, IndexOf(SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonIntensitiesLabel, IndexOf(SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonIntensitiesLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.XLTotalScoreLabel, IndexOf(SpectrumMatchFromTsvHeader.XLTotalScoreLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.ParentIonsLabel, IndexOf(SpectrumMatchFromTsvHeader.ParentIonsLabel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.Ms2ScanRetentionTime, IndexOf(SpectrumMatchFromTsvHeader.Ms2ScanRetentionTime));

            // Glyco
            parsedHeader.Add(SpectrumMatchFromTsvHeader.GlycanMass, IndexOf(SpectrumMatchFromTsvHeader.GlycanMass));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.GlycanStructure, IndexOf(SpectrumMatchFromTsvHeader.GlycanStructure));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.GlycanComposition, IndexOf(SpectrumMatchFromTsvHeader.GlycanComposition));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.LocalizedScores, IndexOf(SpectrumMatchFromTsvHeader.LocalizedScores));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.NumberOfGlycan, IndexOf(SpectrumMatchFromTsvHeader.NumberOfGlycan));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.TotalGlycanSite, IndexOf(SpectrumMatchFromTsvHeader.TotalGlycanSite));
            if (parsedHeader[SpectrumMatchFromTsvHeader.TotalGlycanSite] == -1) // try legacy name from previous versions
            {
                parsedHeader[SpectrumMatchFromTsvHeader.TotalGlycanSite] = IndexOf("Total Glycosylation sites");
            }
            parsedHeader.Add(SpectrumMatchFromTsvHeader.GlycanLocalizationLevel, IndexOf(SpectrumMatchFromTsvHeader.GlycanLocalizationLevel));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.LocalizedGlycanInPeptide, IndexOf(SpectrumMatchFromTsvHeader.LocalizedGlycanInPeptide));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.LocalizedGlycanInProtein, IndexOf(SpectrumMatchFromTsvHeader.LocalizedGlycanInProtein));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.YionScore, IndexOf(SpectrumMatchFromTsvHeader.YionScore));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.DiagonosticIonScore, IndexOf(SpectrumMatchFromTsvHeader.DiagonosticIonScore));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.NGlycanMotifCheck, IndexOf(SpectrumMatchFromTsvHeader.NGlycanMotifCheck));
            if (parsedHeader[SpectrumMatchFromTsvHeader.NGlycanMotifCheck] == -1)// try legacy name from previous versions
            {
                parsedHeader[SpectrumMatchFromTsvHeader.NGlycanMotifCheck] = IndexOf("N-Glycan motif Check");
            }
            parsedHeader.Add(SpectrumMatchFromTsvHeader.R138144, IndexOf(SpectrumMatchFromTsvHeader.R138144));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.AllPotentialGlycanLocalization, IndexOf(SpectrumMatchFromTsvHeader.AllPotentialGlycanLocalization));
            if (parsedHeader[SpectrumMatchFromTsvHeader.AllPotentialGlycanLocalization] == -1)// try legacy name from previous versions
            {
                parsedHeader[SpectrumMatchFromTsvHeader.AllPotentialGlycanLocalization] = IndexOf("All potential glycan localizations");
            }
            parsedHeader.Add(SpectrumMatchFromTsvHeader.AllSiteSpecificLocalizationProbability, IndexOf(SpectrumMatchFromTsvHeader.AllSiteSpecificLocalizationProbability));
            if (parsedHeader[SpectrumMatchFromTsvHeader.AllSiteSpecificLocalizationProbability] == -1)// try legacy name from previous versions
            {
                parsedHeader[SpectrumMatchFromTsvHeader.AllSiteSpecificLocalizationProbability] = IndexOf("AllSiteSpecificLocalizationProbability");
            }

            // Oligo
            parsedHeader.Add(SpectrumMatchFromTsvHeader.FivePrimeTerminus, IndexOf(SpectrumMatchFromTsvHeader.FivePrimeTerminus));
            parsedHeader.Add(SpectrumMatchFromTsvHeader.ThreePrimeTerminus, IndexOf(SpectrumMatchFromTsvHeader.ThreePrimeTerminus));

            return parsedHeader;
        }

        private static bool ResultIsGlyco(Dictionary<string, int> parsedHeader, string line)
        {
            if (!parsedHeader.ContainsKey(SpectrumMatchFromTsvHeader.GlycanMass))
                return false;
            int glycanMassIndex = parsedHeader[SpectrumMatchFromTsvHeader.GlycanMass];
            if (glycanMassIndex < 0)
                return false;

            // Count tabs to check if line has enough columns without splitting
            int tabCount = 0;
            foreach (char c in line)
            {
                if (c == '\t' && ++tabCount >= glycanMassIndex)
                    return true;
            }
            return false;
        }
    }
}
