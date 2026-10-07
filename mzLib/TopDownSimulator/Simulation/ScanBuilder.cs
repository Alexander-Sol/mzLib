using System;
using System.Linq;
using MassSpectrometry;
using MzLibUtil;
using Readers;

namespace TopDownSimulator.Simulation;

/// <summary>
/// Converts a rasterized simulation grid into mzLib scan/file objects.
/// </summary>
public sealed class ScanBuilder
{
    public MsDataScan[] BuildMs1Scans(
        RasterizedScanGrid grid,
        bool isCentroid = false,
        Polarity polarity = Polarity.Positive,
        MZAnalyzerType analyzer = MZAnalyzerType.Orbitrap,
        string scanFilter = "synthetic")
    {
        int nScans = grid.ScanTimes.Length;
        int nMz = grid.MzGrid.Length;
        var scans = new MsDataScan[nScans];
        var mzRange = new MzRange(grid.MzGrid[0], grid.MzGrid[^1]);

        for (int s = 0; s < nScans; s++)
        {
            var intensities = new double[nMz];
            double tic = 0;
            for (int b = 0; b < nMz; b++)
            {
                double intensity = grid.Intensities[s, b];
                intensities[b] = intensity;
                tic += intensity;
            }

            scans[s] = new MsDataScan(
                massSpectrum: new MzSpectrum((double[])grid.MzGrid.Clone(), intensities, false),
                oneBasedScanNumber: s + 1,
                msnOrder: 1,
                isCentroid: isCentroid,
                polarity: polarity,
                retentionTime: grid.ScanTimes[s],
                scanWindowRange: mzRange,
                scanFilter: scanFilter,
                mzAnalyzer: analyzer,
                totalIonCurrent: tic,
                injectionTime: 1.0,
                noiseData: null,
                nativeId: $"scan={s + 1}");
        }

        return scans;
    }

    /// <summary>
    /// Builds centroided MS1 scans from per-scan peak lists, each already in ascending m/z. Every
    /// scan carries the same scan window, the span of peaks over the whole run.
    /// </summary>
    public MsDataScan[] BuildCentroidedMs1Scans(
        double[] scanTimes,
        (double[] Mz, double[] Intensity)[] spectra,
        Polarity polarity = Polarity.Positive,
        MZAnalyzerType analyzer = MZAnalyzerType.Orbitrap,
        string scanFilter = "synthetic")
    {
        if (scanTimes.Length != spectra.Length)
            throw new ArgumentException("There must be one spectrum per scan time.", nameof(spectra));

        double min = double.PositiveInfinity, max = double.NegativeInfinity;
        foreach (var (mz, _) in spectra)
        {
            if (mz.Length == 0) continue;
            min = Math.Min(min, mz[0]);
            max = Math.Max(max, mz[^1]);
        }

        var mzRange = double.IsFinite(min) ? new MzRange(min, max) : new MzRange(0, 0);
        var scans = new MsDataScan[scanTimes.Length];
        for (int s = 0; s < scans.Length; s++)
        {
            var (mz, intensities) = spectra[s];
            scans[s] = new MsDataScan(
                massSpectrum: new MzSpectrum(mz, intensities, false),
                oneBasedScanNumber: s + 1,
                msnOrder: 1,
                isCentroid: true,
                polarity: polarity,
                retentionTime: scanTimes[s],
                scanWindowRange: mzRange,
                scanFilter: scanFilter,
                mzAnalyzer: analyzer,
                totalIonCurrent: intensities.Sum(),
                injectionTime: 1.0,
                noiseData: null,
                nativeId: $"scan={s + 1}");
        }

        return scans;
    }

    public GenericMsDataFile BuildMsDataFile(MsDataScan[] scans, SourceFile? sourceFile = null)
    {
        sourceFile ??= new SourceFile(
            nativeIdFormat: "scan number only nativeID format",
            massSpectrometerFileFormat: "mzML format",
            checkSum: null,
            fileChecksumType: null,
            id: "TopDownSimulator");

        return new GenericMsDataFile(scans, sourceFile);
    }
}
