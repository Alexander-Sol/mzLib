#nullable enable
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;
using MassSpectrometry;
using TopDownSimulator.Model;

namespace TopDownSimulator.Simulation;

/// <summary>One proteoform's monoisotopic mass before and after a precursor shift.</summary>
public sealed record ShiftedProteoform(
    string Identifier,
    double OriginalMonoisotopicMass,
    double ShiftedMonoisotopicMass);

/// <summary>Counts from rewriting a set of MS2 scans onto shifted precursors.</summary>
/// <param name="Shifted">Scans whose precursor m/z fields were moved.</param>
/// <param name="DroppedWithoutChargeState">
/// Scans discarded because the shift is undefined for them: the neutral shift converts to an m/z
/// shift only through the precursor charge, so a scan with no charge state guess cannot be moved.
/// </param>
/// <param name="DroppedWithoutPrecursorMz">
/// Scans discarded because they carry neither a selected-ion m/z nor an isolation m/z, so there is
/// nothing to shift and nothing the mzML writer could emit.
/// </param>
/// <param name="DroppedProfileMode">
/// Profile-mode scans discarded on request. See the remarks on
/// <see cref="PrecursorMassShift.ApplyToMs2Scans"/> for why they cannot simply be written.
/// </param>
public sealed record Ms2ShiftSummary(
    int Shifted,
    int DroppedWithoutChargeState,
    int DroppedWithoutPrecursorMz,
    int DroppedProfileMode);

/// <summary>
/// Applies a single constant neutral mass offset to every simulated proteoform and to the MS2
/// scans that isolate them.
/// </summary>
/// <remarks>
/// A neutral shift of X daltons moves a charge-z precursor by X/z in m/z, so the MS1 envelope and
/// the MS2 isolation window move by different amounts and have to be shifted separately — which is
/// the whole reason this is not a one-line edit of the model masses.
/// <para>
/// Fragment m/z values are deliberately left alone. The shifted file is meant to read as an
/// unlocalized modification of mass X sitting on an otherwise unchanged proteoform, which is the
/// arrangement a decoy or open-search benchmark needs; shifting the fragments too would instead
/// describe a molecule that is uniformly heavier, and nothing would be testable against it.
/// </para>
/// </remarks>
public static class PrecursorMassShift
{
    /// <summary>The m/z offset a neutral shift of <paramref name="massShiftDa"/> produces at <paramref name="charge"/>.</summary>
    /// <remarks>
    /// The magnitude of the charge is used, so this is the offset in the direction the mass moved
    /// regardless of polarity.
    /// </remarks>
    public static double MzShift(double massShiftDa, int charge)
    {
        if (charge == 0)
            throw new ArgumentOutOfRangeException(nameof(charge), "A precursor shift is undefined at charge 0.");

        return massShiftDa / Math.Abs(charge);
    }

    /// <summary>
    /// Returns copies of <paramref name="proteoforms"/> with every monoisotopic mass moved by
    /// <paramref name="massShiftDa"/>. Everything else — abundance, elution, charge distribution —
    /// is carried through untouched.
    /// </summary>
    public static ProteoformModel[] ApplyToModels(
        IReadOnlyList<ProteoformModel> proteoforms, double massShiftDa)
    {
        if (proteoforms is null) throw new ArgumentNullException(nameof(proteoforms));
        if (!double.IsFinite(massShiftDa))
            throw new ArgumentOutOfRangeException(nameof(massShiftDa), massShiftDa, "The mass shift must be finite.");

        var shifted = new ProteoformModel[proteoforms.Count];
        for (int i = 0; i < proteoforms.Count; i++)
        {
            double mass = proteoforms[i].MonoisotopicMass + massShiftDa;
            if (mass <= 0)
                throw new ArgumentOutOfRangeException(nameof(massShiftDa), massShiftDa,
                    $"Shifting '{proteoforms[i].Identifier ?? $"proteoform{i}"}' " +
                    $"({proteoforms[i].MonoisotopicMass:F4} Da) by this much leaves a non-positive mass.");

            shifted[i] = proteoforms[i] with { MonoisotopicMass = mass };
        }

        return shifted;
    }

    /// <summary>
    /// The before/after mass of every proteoform, which is what makes a shifted file interpretable
    /// once it has been written: the ground-truth sidecars describe the file as it stands and so
    /// quote only the shifted masses.
    /// </summary>
    public static ShiftedProteoform[] Describe(
        IReadOnlyList<ProteoformModel> proteoforms, double massShiftDa)
    {
        if (proteoforms is null) throw new ArgumentNullException(nameof(proteoforms));

        var described = new ShiftedProteoform[proteoforms.Count];
        for (int i = 0; i < proteoforms.Count; i++)
        {
            described[i] = new ShiftedProteoform(
                proteoforms[i].Identifier ?? $"proteoform{i}",
                proteoforms[i].MonoisotopicMass,
                proteoforms[i].MonoisotopicMass + massShiftDa);
        }

        return described;
    }

    public static void WriteSidecar(
        IReadOnlyList<ShiftedProteoform> shifts, double massShiftDa, string outputPath)
    {
        if (shifts is null) throw new ArgumentNullException(nameof(shifts));
        if (string.IsNullOrWhiteSpace(outputPath))
            throw new ArgumentException("An output path is required.", nameof(outputPath));

        var sb = new StringBuilder();
        sb.AppendLine(string.Join('\t',
            "Identifier", "OriginalMonoisotopicMass", "ShiftedMonoisotopicMass", "MassShiftDa"));

        foreach (var shift in shifts)
        {
            sb.AppendLine(string.Join('\t',
                shift.Identifier,
                Num(shift.OriginalMonoisotopicMass),
                Num(shift.ShiftedMonoisotopicMass),
                Num(massShiftDa)));
        }

        File.WriteAllText(outputPath, sb.ToString());
    }

    /// <summary>
    /// Rebuilds <paramref name="ms2Scans"/> with their precursor m/z, monoisotopic guess and
    /// isolation window all moved by <paramref name="massShiftDa"/>/z. Peak lists, isolation width
    /// and every other scan property are carried through unchanged.
    /// </summary>
    /// <remarks>
    /// New scan objects are returned rather than the originals being mutated, so a scan list read
    /// from a real file can be shifted without disturbing the file it came from.
    /// <para>
    /// <c>SelectedIonIntensity</c> is carried across as measured. It describes a peak in the real
    /// run, not in the simulated MS1 the shifted scans are about to be written next to, but
    /// recomputing it would mean snapping the selected-ion m/z to the nearest simulated peak — which
    /// would quietly undo the shift wherever the simulation has no peak at the shifted position.
    /// </para>
    /// <para>
    /// Profile-mode scans are rejected. mzLib's mzML reader throws on the first profile spectrum it
    /// meets, so a single one carried across makes the entire written file unreadable — including
    /// the simulated MS1 next to it. Pass <paramref name="dropProfileScans"/> to exclude them
    /// instead, which is a real loss of data and so is not the default.
    /// </para>
    /// </remarks>
    /// <param name="dropProfileScans">
    /// Drop profile-mode scans rather than throwing. They are counted in the returned summary.
    /// </param>
    public static (MsDataScan[] Scans, Ms2ShiftSummary Summary) ApplyToMs2Scans(
        IReadOnlyList<MsDataScan> ms2Scans, double massShiftDa, bool dropProfileScans = false)
    {
        if (ms2Scans is null) throw new ArgumentNullException(nameof(ms2Scans));
        if (!double.IsFinite(massShiftDa))
            throw new ArgumentOutOfRangeException(nameof(massShiftDa), massShiftDa, "The mass shift must be finite.");

        var shifted = new List<MsDataScan>(ms2Scans.Count);
        int droppedWithoutCharge = 0;
        int droppedWithoutPrecursorMz = 0;
        int droppedProfile = 0;

        foreach (var scan in ms2Scans)
        {
            if (scan.MsnOrder < 2)
                throw new ArgumentException(
                    $"Scan {scan.OneBasedScanNumber} is MS{scan.MsnOrder}; only MSn scans carry a precursor to shift.",
                    nameof(ms2Scans));

            if (!scan.IsCentroid)
            {
                if (!dropProfileScans)
                    throw new ArgumentException(
                        $"Scan {scan.OneBasedScanNumber} is profile mode. mzLib's mzML reader rejects the whole " +
                        "file on the first profile spectrum, so carrying it across would make the written run " +
                        "unreadable. Centroid the source scans, or pass dropProfileScans to exclude them.",
                        nameof(ms2Scans));

                droppedProfile++;
                continue;
            }

            int? charge = scan.SelectedIonChargeStateGuess;
            if (charge is null || charge.Value == 0)
            {
                droppedWithoutCharge++;
                continue;
            }

            double? selectedIonMz = scan.SelectedIonMZ ?? scan.IsolationMz;
            if (selectedIonMz is null)
            {
                droppedWithoutPrecursorMz++;
                continue;
            }

            double mzShift = MzShift(massShiftDa, charge.Value);

            shifted.Add(new MsDataScan(
                massSpectrum: scan.MassSpectrum,
                oneBasedScanNumber: scan.OneBasedScanNumber,
                msnOrder: scan.MsnOrder,
                isCentroid: scan.IsCentroid,
                polarity: scan.Polarity,
                retentionTime: scan.RetentionTime,
                scanWindowRange: scan.ScanWindowRange,
                scanFilter: scan.ScanFilter,
                mzAnalyzer: scan.MzAnalyzer,
                totalIonCurrent: scan.TotalIonCurrent,
                injectionTime: scan.InjectionTime,
                noiseData: scan.NoiseData,
                nativeId: scan.NativeId,
                selectedIonMz: selectedIonMz.Value + mzShift,
                selectedIonChargeStateGuess: charge.Value,
                selectedIonIntensity: scan.SelectedIonIntensity,
                isolationMZ: (scan.IsolationMz ?? selectedIonMz.Value) + mzShift,
                isolationWidth: scan.IsolationWidth,
                // The writer dereferences this unconditionally, and Unknown is a term the mzML
                // dictionaries carry, so an unannotated activation survives the round trip.
                dissociationType: scan.DissociationType ?? DissociationType.Unknown,
                oneBasedPrecursorScanNumber: scan.OneBasedPrecursorScanNumber,
                selectedIonMonoisotopicGuessMz: scan.SelectedIonMonoisotopicGuessMz.HasValue
                    ? scan.SelectedIonMonoisotopicGuessMz.Value + mzShift
                    : null,
                hcdEnergy: scan.HcdEnergy,
                scanDescription: scan.ScanDescription,
                compensationVoltage: scan.CompensationVoltage));
        }

        return (shifted.ToArray(),
            new Ms2ShiftSummary(shifted.Count, droppedWithoutCharge, droppedWithoutPrecursorMz, droppedProfile));
    }

    private static string Num(double value) => value.ToString("R", CultureInfo.InvariantCulture);
}
