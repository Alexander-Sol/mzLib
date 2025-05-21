using CsvHelper;
using CsvHelper.Configuration;
using Readers;

namespace Readers
{
    public class ProformaRecord
    {
        public static CsvConfiguration CsvContext => new CsvConfiguration(System.Globalization.CultureInfo.InvariantCulture)
        {
            Delimiter = "\t",
            HasHeaderRecord = true,
        };

        // For This repo

        public string Condition { get; set; }

        // Proforma
        public string FileName { get; set; }
        public double ScanNumber { get; set; }
        public double PrecursorCharge { get; set; }
        public string ProteinAccession { get; set; }
        public string BaseSequence { get; set; }
        public int ModificationMass { get; set; }
        public string FullSequence { get; set; }
    }

    public class ProformaFile : ResultFile<ProformaRecord>, IResultFile
    {

        public ProformaFile(string filepath) : base(filepath) { }

        public override void LoadResults()
        {
            using var csv = new CsvReader(new StreamReader(FilePath), ProformaRecord.CsvContext);
            Results = csv.GetRecords<ProformaRecord>().ToList();
        }

        public override void WriteResults(string outputPath)
        {
            using (var csv = new CsvWriter(new StreamWriter(File.Create(outputPath)), ProformaRecord.CsvContext))
            {
                csv.WriteHeader<ProformaRecord>();
                foreach (var result in Results)
                {
                    csv.NextRecord();
                    csv.WriteRecord(result);
                }
            }
        }

        public override SupportedFileType FileType { get; }
        public override Software Software { get; set; }
    }
}