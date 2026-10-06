#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;

namespace TopDownSimulator.Simulation;

/// <summary>
/// A single acquisition ordered by retention time and numbered contiguously from 1.
/// </summary>
/// <param name="Scans">Every scan, MS1 and MSn interleaved, in acquisition order.</param>
/// <param name="Ms1ScanNumbers">
/// The number each MS1 scan ended up with, parallel to the MS1 scans that were passed in. The
/// feature-level ground truth is written against these, so it has to be the post-merge numbering
/// rather than the numbering the MS1 simulation produced on its own.
/// </param>
/// <param name="Ms2ScansDroppedOutsideMs1TimeRange">
/// MSn scans discarded because no MS1 scan precedes them in the merged file. An MSn scan whose
/// <c>spectrumRef</c> points outside the file is not something a reader can resolve.
/// </param>
/// <param name="Ms2ScansDroppedWithEmptyPrecursor">
/// MSn scans discarded because the survey scan preceding them carries no peaks. See the remarks on
/// <see cref="ScanListMerger"/> for why writing these is worse than dropping them.
/// </param>
public sealed record MergedScanList(
    MsDataScan[] Scans,
    int[] Ms1ScanNumbers,
    int Ms2ScansDroppedOutsideMs1TimeRange,
    int Ms2ScansDroppedWithEmptyPrecursor);

/// <summary>
/// Interleaves simulated MS1 scans with MSn scans carried over from another run, renumbering both
/// into one contiguous acquisition.
/// </summary>
/// <remarks>
/// Renumbering is not cosmetic. <c>MsDataFile.GetOneBasedScan(n)</c> is <c>Scans[n-1]</c>, and the
/// mzML writer resolves an MSn scan's precursor reference through it, so scan numbers must equal
/// positions or every <c>spectrumRef</c> in the file points at the wrong spectrum.
/// <para>
/// An MSn scan whose preceding survey scan has no peaks is dropped rather than written. A clean
/// simulation carries peaks only where a proteoform is eluting, so most of its survey scans reduce
/// to nothing — and an empty precursor spectrum is fatal downstream, not merely unhelpful:
/// <see cref="MsDataScan.RefineSelectedMzAndIntensity"/> throws on one, and a consumer that catches
/// that per scan and moves on (MetaMorpheus does) can be left with an unassigned slot that fails
/// far from the cause, taking the whole file with it. Linking backwards to an earlier non-empty
/// survey scan would be worse: that scan's peaks belong to whatever was eluting *then*, so refining
/// the precursor against it fabricates a precursor rather than losing one.
/// </para>
/// </remarks>
public static class ScanListMerger
{
    /// <summary>
    /// Orders <paramref name="ms1Scans"/> and <paramref name="ms2Scans"/> by retention time,
    /// renumbers them from 1, and links each MSn scan to the MS1 scan that most recently preceded
    /// it.
    /// </summary>
    /// <remarks>
    /// The scans passed in are renumbered in place — their scan number, native id and precursor
    /// link are rewritten — so callers must own them. Both
    /// <see cref="Simulator"/> and <see cref="PrecursorMassShift.ApplyToMs2Scans"/> hand back fresh
    /// objects for exactly this reason.
    /// <para>
    /// MS1 wins ties in retention time, which puts an MSn scan after the survey scan it was
    /// triggered from when an instrument reports both at the same stamp.
    /// </para>
    /// </remarks>
    public static MergedScanList Merge(
        IReadOnlyList<MsDataScan> ms1Scans, IReadOnlyList<MsDataScan> ms2Scans)
    {
        if (ms1Scans is null) throw new ArgumentNullException(nameof(ms1Scans));
        if (ms2Scans is null) throw new ArgumentNullException(nameof(ms2Scans));
        if (ms1Scans.Count == 0)
            throw new ArgumentException("At least one MS1 scan is required to merge against.", nameof(ms1Scans));

        // Indices rather than the scans themselves, so the returned numbers can be reported against
        // the caller's ordering — the feature truth indexes MS1 scans by their position in the scan
        // time grid, not by retention time rank.
        var ms1Order = Enumerable.Range(0, ms1Scans.Count)
            .OrderBy(k => ms1Scans[k].RetentionTime)
            .ToArray();

        double firstMs1Rt = ms1Scans[ms1Order[0]].RetentionTime;
        double lastMs1Rt = ms1Scans[ms1Order[^1]].RetentionTime;

        var ms2 = ms2Scans
            .Where(s => s.RetentionTime >= firstMs1Rt && s.RetentionTime <= lastMs1Rt)
            .OrderBy(s => s.RetentionTime)
            .ToArray();
        int droppedOutsideRange = ms2Scans.Count - ms2.Length;

        var merged = new List<MsDataScan>(ms1Scans.Count + ms2.Length);
        var ms1ScanNumbers = new int[ms1Scans.Count];

        int i = 0, j = 0, nextNumber = 1, lastMs1Number = 0;
        bool lastMs1IsEmpty = true;
        int droppedEmptyPrecursor = 0;

        while (i < ms1Order.Length || j < ms2.Length)
        {
            bool takeMs1 = j >= ms2.Length
                || (i < ms1Order.Length
                    && ms1Scans[ms1Order[i]].RetentionTime <= ms2[j].RetentionTime);

            if (takeMs1)
            {
                var scan = ms1Scans[ms1Order[i]];
                scan.SetOneBasedScanNumber(nextNumber);
                scan.SetNativeID($"scan={nextNumber}");
                ms1ScanNumbers[ms1Order[i]] = nextNumber;
                lastMs1Number = nextNumber;
                lastMs1IsEmpty = scan.MassSpectrum.XArray.Length == 0;
                merged.Add(scan);
                i++;
            }
            else
            {
                var scan = ms2[j];
                j++;

                if (lastMs1IsEmpty)
                {
                    droppedEmptyPrecursor++;
                    continue;
                }

                scan.SetOneBasedScanNumber(nextNumber);
                scan.SetNativeID($"scan={nextNumber}");
                scan.SetOneBasedPrecursorScanNumber(lastMs1Number);
                merged.Add(scan);
            }

            nextNumber++;
        }

        return new MergedScanList(
            merged.ToArray(), ms1ScanNumbers, droppedOutsideRange, droppedEmptyPrecursor);
    }
}
