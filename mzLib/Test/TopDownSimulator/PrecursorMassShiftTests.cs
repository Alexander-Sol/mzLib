using System;
using System.Globalization;
using System.IO;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Readers;
using TopDownSimulator.Model;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

[TestFixture]
public class PrecursorMassShiftTests
{
    private const double Mass = 10000.0;
    private const int MinCharge = 6;
    private const int MaxCharge = 11;
    private const double SigmaMz = 0.012;
    private const double MassShift = 12.5;

    private string _outputDirectory;

    private static ProteoformModel BuildModel(
        double mass = Mass, double abundance = 1.5e6, double rtCenter = 20.0, string identifier = "P00000|TESTPROTEOFORM") =>
        new ProteoformModel(
            MonoisotopicMass: mass,
            Abundance: abundance,
            RtProfile: new EmgProfile(Mu: rtCenter, Sigma: 0.22, Tau: 0.08),
            ChargeDistribution: new GaussianChargeDistribution(MuZ: 8.3, SigmaZ: 1.15),
            Identifier: identifier);

    private static double[] ScanTimes(int count = 11, double start = 19.5, double step = 0.1) =>
        Enumerable.Range(0, count).Select(i => start + i * step).ToArray();

    /// <summary>
    /// An MS2 scan shaped like one a Thermo reader hands back: real fragment peaks, a charge state
    /// guess, and a populated isolation window.
    /// </summary>
    private static MsDataScan BuildMs2Scan(
        double retentionTime,
        int charge = 8,
        double precursorMass = Mass,
        int oneBasedScanNumber = 1,
        double isolationWidth = 4.0,
        int? chargeOverride = null,
        bool omitCharge = false)
    {
        double precursorMz = precursorMass.ToMz(charge);
        var mz = new[] { 500.25, 700.4, 900.55, 1100.7 };
        var intensities = new[] { 1000.0, 2000.0, 1500.0, 800.0 };

        return new MsDataScan(
            massSpectrum: new MzSpectrum(mz, intensities, false),
            oneBasedScanNumber: oneBasedScanNumber,
            msnOrder: 2,
            isCentroid: true,
            polarity: Polarity.Positive,
            retentionTime: retentionTime,
            scanWindowRange: new MzRange(200, 2000),
            scanFilter: "real ms2",
            mzAnalyzer: MZAnalyzerType.Orbitrap,
            totalIonCurrent: intensities.Sum(),
            injectionTime: 25.0,
            noiseData: null,
            nativeId: $"controllerType=0 controllerNumber=1 scan={oneBasedScanNumber}",
            selectedIonMz: precursorMz,
            selectedIonChargeStateGuess: omitCharge ? null : chargeOverride ?? charge,
            selectedIonIntensity: 5e5,
            isolationMZ: precursorMz,
            isolationWidth: isolationWidth,
            dissociationType: DissociationType.ETD,
            oneBasedPrecursorScanNumber: null,
            selectedIonMonoisotopicGuessMz: precursorMz);
    }

    private static MsDataScan BuildMs1Scan(double retentionTime, int oneBasedScanNumber) =>
        new MsDataScan(
            massSpectrum: new MzSpectrum(new[] { 600.0, 700.0 }, new[] { 10.0, 20.0 }, false),
            oneBasedScanNumber: oneBasedScanNumber,
            msnOrder: 1,
            isCentroid: true,
            polarity: Polarity.Positive,
            retentionTime: retentionTime,
            scanWindowRange: new MzRange(200, 2000),
            scanFilter: "synthetic",
            mzAnalyzer: MZAnalyzerType.Orbitrap,
            totalIonCurrent: 30.0,
            injectionTime: 1.0,
            noiseData: null,
            nativeId: $"scan={oneBasedScanNumber}");

    [SetUp]
    public void SetUp()
    {
        _outputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownSimulatorShiftOutput");
        Directory.CreateDirectory(_outputDirectory);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_outputDirectory))
            Directory.Delete(_outputDirectory, recursive: true);
    }

    [Test]
    public void MzShiftIsTheNeutralShiftDividedByChargeMagnitude()
    {
        Assert.That(PrecursorMassShift.MzShift(10.0, 5), Is.EqualTo(2.0).Within(1e-12));
        Assert.That(PrecursorMassShift.MzShift(-10.0, 4), Is.EqualTo(-2.5).Within(1e-12));

        // Polarity changes the sign of the charge, not the direction the mass moved.
        Assert.That(PrecursorMassShift.MzShift(10.0, -5), Is.EqualTo(2.0).Within(1e-12));

        Assert.Throws<ArgumentOutOfRangeException>(() => PrecursorMassShift.MzShift(10.0, 0));
    }

    [Test]
    public void ApplyToModelsMovesTheMassAndNothingElse()
    {
        var original = BuildModel();
        var shifted = PrecursorMassShift.ApplyToModels(new[] { original }, MassShift);

        Assert.That(shifted, Has.Length.EqualTo(1));
        Assert.That(shifted[0].MonoisotopicMass, Is.EqualTo(Mass + MassShift).Within(1e-9));
        Assert.That(shifted[0].Abundance, Is.EqualTo(original.Abundance));
        Assert.That(shifted[0].RtProfile, Is.EqualTo(original.RtProfile));
        Assert.That(shifted[0].ChargeDistribution, Is.SameAs(original.ChargeDistribution));
        Assert.That(shifted[0].Identifier, Is.EqualTo(original.Identifier));

        Assert.That(original.MonoisotopicMass, Is.EqualTo(Mass), "The input models must not be mutated.");
    }

    [Test]
    public void ApplyToModelsRefusesAShiftThatLeavesANonPositiveMass()
    {
        var ex = Assert.Throws<ArgumentOutOfRangeException>(
            () => PrecursorMassShift.ApplyToModels(new[] { BuildModel(mass: 5000.0) }, -5000.0));

        Assert.That(ex.Message, Does.Contain("non-positive"));
    }

    [Test]
    public void ApplyToMs2ScansMovesEveryPrecursorFieldButLeavesTheFragments()
    {
        const int charge = 8;
        var source = BuildMs2Scan(retentionTime: 20.0, charge: charge);
        double expectedShift = MassShift / charge;

        var (shifted, summary) = PrecursorMassShift.ApplyToMs2Scans(new[] { source }, MassShift);

        Assert.That(summary.Shifted, Is.EqualTo(1));
        Assert.That(summary.DroppedWithoutChargeState, Is.Zero);
        Assert.That(summary.DroppedWithoutPrecursorMz, Is.Zero);

        var scan = shifted[0];
        Assert.That(scan.SelectedIonMZ.Value,
            Is.EqualTo(source.SelectedIonMZ.Value + expectedShift).Within(1e-9));
        Assert.That(scan.IsolationMz.Value,
            Is.EqualTo(source.IsolationMz.Value + expectedShift).Within(1e-9));
        Assert.That(scan.SelectedIonMonoisotopicGuessMz.Value,
            Is.EqualTo(source.SelectedIonMonoisotopicGuessMz.Value + expectedShift).Within(1e-9));

        // The window moves; it does not widen.
        Assert.That(scan.IsolationWidth.Value, Is.EqualTo(source.IsolationWidth.Value).Within(1e-12));
        Assert.That(scan.IsolationRange.Width, Is.EqualTo(source.IsolationRange.Width).Within(1e-9));
        Assert.That(scan.IsolationRange.Mean,
            Is.EqualTo(source.IsolationRange.Mean + expectedShift).Within(1e-9));

        // The whole point of the arrangement: the fragments still describe the unshifted proteoform.
        Assert.That(scan.MassSpectrum.XArray, Is.EqualTo(source.MassSpectrum.XArray));
        Assert.That(scan.MassSpectrum.YArray, Is.EqualTo(source.MassSpectrum.YArray));

        Assert.That(scan.SelectedIonChargeStateGuess, Is.EqualTo(charge));
        Assert.That(scan.DissociationType, Is.EqualTo(DissociationType.ETD));
        Assert.That(scan.RetentionTime, Is.EqualTo(source.RetentionTime));
    }

    [Test]
    public void ApplyToMs2ScansDoesNotDisturbTheScansItWasGiven()
    {
        var source = BuildMs2Scan(retentionTime: 20.0);
        double originalSelectedIonMz = source.SelectedIonMZ.Value;

        PrecursorMassShift.ApplyToMs2Scans(new[] { source }, MassShift);

        Assert.That(source.SelectedIonMZ.Value, Is.EqualTo(originalSelectedIonMz).Within(1e-12));
    }

    [Test]
    public void ScansWithNoChargeStateAreDroppedBecauseTheShiftIsUndefinedForThem()
    {
        var scans = new[]
        {
            BuildMs2Scan(retentionTime: 20.0, oneBasedScanNumber: 1),
            BuildMs2Scan(retentionTime: 20.1, oneBasedScanNumber: 2, omitCharge: true),
        };

        var (shifted, summary) = PrecursorMassShift.ApplyToMs2Scans(scans, MassShift);

        Assert.That(shifted, Has.Length.EqualTo(1));
        Assert.That(summary.Shifted, Is.EqualTo(1));
        Assert.That(summary.DroppedWithoutChargeState, Is.EqualTo(1));
    }

    /// <summary>
    /// A single profile spectrum makes mzLib's reader throw on the whole file, so carrying one
    /// across would take the simulated MS1 down with it. Failing at the shift is the only place the
    /// problem is still legible.
    /// </summary>
    [Test]
    public void ProfileModeScansAreRefusedUnlessDroppingIsAskedFor()
    {
        var profile = new MsDataScan(
            massSpectrum: new MzSpectrum(new[] { 500.0, 500.01 }, new[] { 5.0, 9.0 }, false),
            oneBasedScanNumber: 7,
            msnOrder: 2,
            isCentroid: false,
            polarity: Polarity.Positive,
            retentionTime: 20.0,
            scanWindowRange: new MzRange(200, 2000),
            scanFilter: "real ms2",
            mzAnalyzer: MZAnalyzerType.Orbitrap,
            totalIonCurrent: 14.0,
            injectionTime: 25.0,
            noiseData: null,
            nativeId: "scan=7",
            selectedIonMz: 1251.0,
            selectedIonChargeStateGuess: 8,
            selectedIonIntensity: 5e5,
            isolationMZ: 1251.0,
            isolationWidth: 4.0,
            dissociationType: DissociationType.ETD);

        var ex = Assert.Throws<ArgumentException>(
            () => PrecursorMassShift.ApplyToMs2Scans(new[] { profile }, MassShift));
        Assert.That(ex.Message, Does.Contain("profile mode"));

        var (shifted, summary) = PrecursorMassShift.ApplyToMs2Scans(
            new[] { profile, BuildMs2Scan(retentionTime: 20.1, oneBasedScanNumber: 8) },
            MassShift,
            dropProfileScans: true);

        Assert.That(shifted, Has.Length.EqualTo(1));
        Assert.That(summary.DroppedProfileMode, Is.EqualTo(1));
    }

    [Test]
    public void ApplyToMs2ScansRejectsAnMs1()
    {
        Assert.Throws<ArgumentException>(
            () => PrecursorMassShift.ApplyToMs2Scans(new[] { BuildMs1Scan(20.0, 1) }, MassShift));
    }

    [Test]
    public void MergeInterleavesByRetentionTimeAndNumbersScansContiguously()
    {
        var ms1 = new[] { BuildMs1Scan(20.0, 1), BuildMs1Scan(20.2, 2), BuildMs1Scan(20.4, 3) };
        var ms2 = new[]
        {
            BuildMs2Scan(retentionTime: 20.15, oneBasedScanNumber: 900),
            BuildMs2Scan(retentionTime: 20.05, oneBasedScanNumber: 901),
            BuildMs2Scan(retentionTime: 20.25, oneBasedScanNumber: 902),
        };

        var merged = ScanListMerger.Merge(ms1, ms2);

        Assert.That(merged.Scans.Select(s => s.OneBasedScanNumber), Is.EqualTo(new[] { 1, 2, 3, 4, 5, 6 }));
        Assert.That(merged.Scans.Select(s => s.MsnOrder), Is.EqualTo(new[] { 1, 2, 2, 1, 2, 1 }));
        Assert.That(merged.Scans.Select(s => s.RetentionTime), Is.Ordered.Ascending);
        Assert.That(merged.Ms1ScanNumbers, Is.EqualTo(new[] { 1, 4, 6 }));
        Assert.That(merged.Ms2ScansDroppedOutsideMs1TimeRange, Is.Zero);

        // Scan number and native id have to agree, because the writer resolves precursor references
        // by native id but indexes scans by number.
        foreach (var scan in merged.Scans)
            Assert.That(scan.NativeId, Is.EqualTo($"scan={scan.OneBasedScanNumber}"));
    }

    [Test]
    public void EveryMergedMs2LinksToTheMostRecentPrecedingMs1()
    {
        var ms1 = new[] { BuildMs1Scan(20.0, 1), BuildMs1Scan(20.2, 2), BuildMs1Scan(20.4, 3) };
        var ms2 = new[]
        {
            BuildMs2Scan(retentionTime: 20.05, oneBasedScanNumber: 900),
            BuildMs2Scan(retentionTime: 20.25, oneBasedScanNumber: 901),
            BuildMs2Scan(retentionTime: 20.35, oneBasedScanNumber: 902),
        };

        var merged = ScanListMerger.Merge(ms1, ms2);

        foreach (var scan in merged.Scans.Where(s => s.MsnOrder > 1))
        {
            int precursorNumber = scan.OneBasedPrecursorScanNumber.Value;
            var precursor = merged.Scans[precursorNumber - 1];

            Assert.That(precursor.MsnOrder, Is.EqualTo(1));
            Assert.That(precursor.RetentionTime, Is.LessThanOrEqualTo(scan.RetentionTime));
            Assert.That(precursorNumber, Is.LessThan(scan.OneBasedScanNumber));

            // Nothing between them may be a survey scan, or this is not the *most recent* MS1.
            for (int n = precursorNumber + 1; n < scan.OneBasedScanNumber; n++)
                Assert.That(merged.Scans[n - 1].MsnOrder, Is.GreaterThan(1));
        }
    }

    /// <summary>
    /// A clean simulation carries peaks only where a proteoform elutes, so most survey scans reduce
    /// to nothing. An MS2 behind one of them cannot have its precursor refined —
    /// <see cref="MsDataScan.RefineSelectedMzAndIntensity"/> throws on an empty precursor spectrum —
    /// and a consumer that catches that per scan and continues can end up failing on the whole file.
    /// </summary>
    [Test]
    public void Ms2ScansBehindAnEmptyPrecursorScanAreDropped()
    {
        var empty = new MsDataScan(
            massSpectrum: new MzSpectrum(Array.Empty<double>(), Array.Empty<double>(), false),
            oneBasedScanNumber: 2,
            msnOrder: 1,
            isCentroid: true,
            polarity: Polarity.Positive,
            retentionTime: 20.2,
            scanWindowRange: new MzRange(200, 2000),
            scanFilter: "synthetic",
            mzAnalyzer: MZAnalyzerType.Orbitrap,
            totalIonCurrent: 0,
            injectionTime: 1.0,
            noiseData: null,
            nativeId: "scan=2");

        // The trailing survey scan keeps the last MS2 inside the MS1 time range, so the only reason
        // anything is dropped here is the empty precursor.
        var ms1 = new[] { BuildMs1Scan(20.0, 1), empty, BuildMs1Scan(20.4, 3), BuildMs1Scan(20.6, 4) };
        var ms2 = new[]
        {
            BuildMs2Scan(retentionTime: 20.1, oneBasedScanNumber: 900),  // behind a populated MS1
            BuildMs2Scan(retentionTime: 20.3, oneBasedScanNumber: 901),  // behind the empty MS1
            BuildMs2Scan(retentionTime: 20.5, oneBasedScanNumber: 902),  // behind a populated MS1
        };

        var merged = ScanListMerger.Merge(ms1, ms2);

        Assert.That(merged.Ms2ScansDroppedWithEmptyPrecursor, Is.EqualTo(1));
        Assert.That(merged.Ms2ScansDroppedOutsideMs1TimeRange, Is.Zero);
        Assert.That(merged.Scans.Count(s => s.MsnOrder > 1), Is.EqualTo(2));

        // Numbering must stay contiguous across the gap the dropped scan left behind.
        Assert.That(merged.Scans.Select(s => s.OneBasedScanNumber),
            Is.EqualTo(Enumerable.Range(1, merged.Scans.Length)));

        foreach (var scan in merged.Scans.Where(s => s.MsnOrder > 1))
        {
            var precursor = merged.Scans[scan.OneBasedPrecursorScanNumber.Value - 1];
            Assert.That(precursor.MassSpectrum.XArray, Is.Not.Empty,
                "No written MSn scan may point at an empty survey scan.");
        }

        // The empty survey scan itself still belongs in the file, and the feature truth still has to
        // be able to name it.
        Assert.That(merged.Scans.Count(s => s.MsnOrder == 1), Is.EqualTo(4));
        Assert.That(merged.Ms1ScanNumbers, Has.Length.EqualTo(4));
        Assert.That(merged.Ms1ScanNumbers, Is.Unique);
    }

    [Test]
    public void Ms2ScansOutsideTheMs1TimeRangeAreDropped()
    {
        var ms1 = new[] { BuildMs1Scan(20.0, 1), BuildMs1Scan(20.2, 2) };
        var ms2 = new[]
        {
            BuildMs2Scan(retentionTime: 19.5, oneBasedScanNumber: 900),  // before the first MS1
            BuildMs2Scan(retentionTime: 20.1, oneBasedScanNumber: 901),
            BuildMs2Scan(retentionTime: 25.0, oneBasedScanNumber: 902),  // after the last MS1
        };

        var merged = ScanListMerger.Merge(ms1, ms2);

        Assert.That(merged.Ms2ScansDroppedOutsideMs1TimeRange, Is.EqualTo(2));
        Assert.That(merged.Scans.Count(s => s.MsnOrder > 1), Is.EqualTo(1));
        Assert.That(merged.Scans.All(s => s.OneBasedPrecursorScanNumber is null or > 0), Is.True);
    }

    [Test]
    public void WrittenRunHasShiftedMs1EnvelopesAndShiftedIsolationWindows()
    {
        const int precursorCharge = 8;
        var model = BuildModel();
        double[] scanTimes = ScanTimes();
        string path = Path.Combine(_outputDirectory, "shifted.mzML");

        var sourceMs2 = new[]
        {
            BuildMs2Scan(retentionTime: scanTimes[3] + 0.01, charge: precursorCharge, oneBasedScanNumber: 500),
            BuildMs2Scan(retentionTime: scanTimes[6] + 0.01, charge: precursorCharge, oneBasedScanNumber: 501),
        };

        var export = new Simulator().WriteShiftedMzml(
            new[] { model }, MinCharge, MaxCharge, SigmaMz, scanTimes, sourceMs2, MassShift, path);

        Assert.That(export.Shift, Is.Not.Null);
        Assert.That(export.Shift.MassShiftDa, Is.EqualTo(MassShift));
        Assert.That(export.Shift.Ms2ScansWritten, Is.EqualTo(2));
        Assert.That(export.ScanCount, Is.EqualTo(scanTimes.Length + 2));

        var roundTripped = MsDataFileReader.GetDataFile(path);
        roundTripped.LoadAllStaticData();
        var scans = roundTripped.GetAllScansList();

        Assert.That(scans, Has.Count.EqualTo(scanTimes.Length + 2));
        Assert.That(scans.Select(s => s.OneBasedScanNumber), Is.Ordered.Ascending);
        Assert.That(scans.Count(s => s.MsnOrder == 1), Is.EqualTo(scanTimes.Length));

        // MS1: the envelope sits at the shifted mass and nothing is left at the original one.
        var apex = scans.Where(s => s.MsnOrder == 1).OrderByDescending(s => s.TotalIonCurrent).First();
        var tolerance = new PpmTolerance(5);
        double shiftedMonoMz = (Mass + MassShift).ToMz(precursorCharge);
        double originalMonoMz = Mass.ToMz(precursorCharge);

        Assert.That(apex.MassSpectrum.XArray.Any(mz => tolerance.Within(mz, shiftedMonoMz)), Is.True,
            $"No MS1 peak at the shifted monoisotopic m/z {shiftedMonoMz}.");
        Assert.That(apex.MassSpectrum.XArray.Any(mz => tolerance.Within(mz, originalMonoMz)), Is.False,
            $"An MS1 peak survived at the unshifted monoisotopic m/z {originalMonoMz}.");

        // MS2: the window moved by exactly mass/charge, and the fragments did not move at all.
        double expectedIsolationMz = originalMonoMz + MassShift / precursorCharge;
        foreach (var ms2 in scans.Where(s => s.MsnOrder > 1))
        {
            Assert.That(ms2.IsolationMz.Value, Is.EqualTo(expectedIsolationMz).Within(1e-4));
            Assert.That(ms2.SelectedIonMZ.Value, Is.EqualTo(expectedIsolationMz).Within(1e-4));
            Assert.That(ms2.IsolationWidth.Value, Is.EqualTo(4.0).Within(1e-6));
            Assert.That(ms2.MassSpectrum.XArray,
                Is.EqualTo(sourceMs2[0].MassSpectrum.XArray).Within(1e-6).AsCollection);

            int precursorNumber = ms2.OneBasedPrecursorScanNumber.Value;
            Assert.That(scans[precursorNumber - 1].MsnOrder, Is.EqualTo(1),
                "An MSn scan's precursor reference must resolve to a survey scan in the same file.");
        }
    }

    [Test]
    public void SidecarsDescribeTheShiftedFileAndRecordWhereItCameFrom()
    {
        var model = BuildModel();
        double[] scanTimes = ScanTimes();
        string path = Path.Combine(_outputDirectory, "shifted-sidecars.mzML");

        var export = new Simulator().WriteShiftedMzml(
            new[] { model }, MinCharge, MaxCharge, SigmaMz, scanTimes,
            new[] { BuildMs2Scan(retentionTime: scanTimes[5] + 0.01) }, MassShift, path);

        // The mass-shift sidecar is the only place the original mass survives.
        Assert.That(export.Shift.ShiftSidecarPath, Is.Not.Null);
        var shiftLines = File.ReadAllLines(export.Shift.ShiftSidecarPath);
        var shiftHeader = shiftLines[0].Split('\t');
        var shiftRow = shiftLines[1].Split('\t');

        Assert.That(shiftLines, Has.Length.EqualTo(2));
        Assert.That(shiftRow[Array.IndexOf(shiftHeader, "Identifier")], Is.EqualTo(model.Identifier));
        Assert.That(Parse(shiftRow[Array.IndexOf(shiftHeader, "OriginalMonoisotopicMass")]),
            Is.EqualTo(Mass).Within(1e-6));
        Assert.That(Parse(shiftRow[Array.IndexOf(shiftHeader, "ShiftedMonoisotopicMass")]),
            Is.EqualTo(Mass + MassShift).Within(1e-6));
        Assert.That(Parse(shiftRow[Array.IndexOf(shiftHeader, "MassShiftDa")]),
            Is.EqualTo(MassShift).Within(1e-9));

        // The feature truth describes the file as written, so it quotes shifted masses.
        var featureLines = File.ReadAllLines(export.FeatureGroundTruthPath);
        var featureHeader = featureLines[0].Split('\t');
        int massColumn = Array.IndexOf(featureHeader, "MonoisotopicMass");
        int chargeColumn = Array.IndexOf(featureHeader, "Charge");
        int monoMzColumn = Array.IndexOf(featureHeader, "MonoisotopicMz");

        Assert.That(featureLines, Has.Length.GreaterThan(1));
        foreach (var row in featureLines.Skip(1).Select(l => l.Split('\t')))
        {
            Assert.That(Parse(row[massColumn]), Is.EqualTo(Mass + MassShift).Within(1e-6));
            Assert.That(Parse(row[monoMzColumn]),
                Is.EqualTo((Mass + MassShift).ToMz(int.Parse(row[chargeColumn]))).Within(1e-6));
        }
    }

    /// <summary>
    /// The feature truth indexes scans by number, and merging renumbers every MS1 to make room for
    /// the MSn scans. If the truth were written against the pre-merge numbering it would point at
    /// MS2 scans.
    /// </summary>
    [Test]
    public void FeatureTruthScanNumbersPointAtMs1ScansInTheMergedFile()
    {
        var model = BuildModel();
        double[] scanTimes = ScanTimes();
        string path = Path.Combine(_outputDirectory, "shifted-scannumbers.mzML");

        // One MS2 after every MS1, so pre-merge and post-merge numbering cannot coincide.
        var sourceMs2 = scanTimes
            .Select((t, i) => BuildMs2Scan(retentionTime: t + 0.01, oneBasedScanNumber: 500 + i))
            .ToArray();

        var export = new Simulator().WriteShiftedMzml(
            new[] { model }, MinCharge, MaxCharge, SigmaMz, scanTimes, sourceMs2, MassShift, path);

        var roundTripped = MsDataFileReader.GetDataFile(path);
        roundTripped.LoadAllStaticData();
        var scans = roundTripped.GetAllScansList();

        var featureLines = File.ReadAllLines(export.FeatureGroundTruthPath);
        var header = featureLines[0].Split('\t');
        int apexScanColumn = Array.IndexOf(header, "ApexScanNumber");
        int firstScanColumn = Array.IndexOf(header, "FirstScanNumber");
        int lastScanColumn = Array.IndexOf(header, "LastScanNumber");
        int apexRtColumn = Array.IndexOf(header, "ApexRt");

        Assert.That(featureLines, Has.Length.GreaterThan(1));
        foreach (var row in featureLines.Skip(1).Select(l => l.Split('\t')))
        {
            foreach (int column in new[] { apexScanColumn, firstScanColumn, lastScanColumn })
            {
                int scanNumber = int.Parse(row[column]);
                Assert.That(scans[scanNumber - 1].MsnOrder, Is.EqualTo(1),
                    $"Scan {scanNumber} quoted by the feature truth is an MS{scans[scanNumber - 1].MsnOrder}.");
            }

            Assert.That(scans[int.Parse(row[apexScanColumn]) - 1].RetentionTime,
                Is.EqualTo(Parse(row[apexRtColumn])).Within(1e-6));
        }
    }

    private static double Parse(string value) =>
        double.Parse(value, NumberStyles.Float, CultureInfo.InvariantCulture);
}
