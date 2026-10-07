using System;
using System.Collections.Generic;
using System.Linq;
using Readers;

namespace TopDownSimulator.Extraction;

public sealed record MmResultRecord(
    string FileNameWithoutExtension,
    int PrecursorScanNumber,
    int Ms2ScanNumber,
    int PrecursorCharge,
    double MonoisotopicMass,
    double RetentionTime,
    double Score,
    string FullSequence,
    string? Accession,
    string Identifier,
    double? PrecursorIntensity = null,
    double? PrecursorMass = null);

/// <summary>
/// Loads MetaMorpheus `.psmtsv` search results into a compact record shape the
/// simulator can use for anchoring extraction and Phase 3 co-eluter discovery.
/// </summary>
public sealed class MmResultLoader
{
    public IReadOnlyList<MmResultRecord> Load(string psmTsvPath)
    {
        var file = new PsmFromTsvFile(psmTsvPath, new SpectrumMatchParsingParameters
        {
            ParseMatchedFragmentIons = false,
        });
        file.LoadResults();

        return file.Results
            .Where(p => p is not null)
            .Where(p => p.MonoisotopicMass > 0 && p.RetentionTime >= 0)
            .Select(p => new MmResultRecord(
                FileNameWithoutExtension: p.FileNameWithoutExtension,
                PrecursorScanNumber: p.PrecursorScanNum,
                Ms2ScanNumber: p.Ms2ScanNumber,
                PrecursorCharge: p.PrecursorCharge,
                MonoisotopicMass: p.MonoisotopicMass,
                RetentionTime: p.RetentionTime,
                Score: p.Score,
                FullSequence: p.FullSequence,
                Accession: p.Accession,
                Identifier: BuildIdentifier(p),
                PrecursorIntensity: p.PrecursorIntensity,
                PrecursorMass: p.PrecursorMass))
            .OrderBy(p => p.FileNameWithoutExtension, StringComparer.OrdinalIgnoreCase)
            .ThenBy(p => p.RetentionTime)
            .ThenBy(p => p.MonoisotopicMass)
            .ToArray();
    }

    /// <summary>
    /// Identifications from one raw file at or below <paramref name="maxQValue"/>, best score first.
    /// </summary>
    /// <param name="fileNameWithoutExtension">The raw file the identifications must come from.</param>
    public IReadOnlyList<MmResultRecord> LoadQualified(string psmTsvPath, string fileNameWithoutExtension, double maxQValue)
    {
        var file = new PsmFromTsvFile(psmTsvPath, new SpectrumMatchParsingParameters
        {
            ParseMatchedFragmentIons = false,
        });
        file.LoadResults();

        return file.Results
            .Where(p => p is not null)
            .Where(p => string.Equals(p.FileNameWithoutExtension, fileNameWithoutExtension, StringComparison.OrdinalIgnoreCase))
            .Where(p => p.MonoisotopicMass > 0 && p.RetentionTime >= 0)
            .Where(p => !double.IsNaN(p.QValue) && p.QValue <= maxQValue)
            .Select(p => new MmResultRecord(
                FileNameWithoutExtension: p.FileNameWithoutExtension,
                PrecursorScanNumber: p.PrecursorScanNum,
                Ms2ScanNumber: p.Ms2ScanNumber,
                PrecursorCharge: p.PrecursorCharge,
                MonoisotopicMass: p.MonoisotopicMass,
                RetentionTime: p.RetentionTime,
                Score: p.Score,
                FullSequence: p.FullSequence,
                Accession: p.Accession,
                Identifier: BuildIdentifier(p),
                PrecursorIntensity: p.PrecursorIntensity,
                PrecursorMass: p.PrecursorMass))
            .OrderByDescending(r => r.Score)
            .ThenBy(r => r.RetentionTime)
            .ToArray();
    }

    private static string BuildIdentifier(PsmFromTsv psm)
    {
        if (!string.IsNullOrWhiteSpace(psm.Accession))
            return $"{psm.Accession}:{psm.FullSequence}:{psm.Ms2ScanNumber}";

        return $"{psm.FileNameWithoutExtension}:{psm.Ms2ScanNumber}";
    }
}
