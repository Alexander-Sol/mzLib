using System;
using System.IO;
using System.Linq;
using NUnit.Framework;
using Readers;
using TopDownSimulator.Model;
using TopDownSimulator.Noise;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

/// <summary>
/// End-to-end checks that a simulation with an injected noise floor writes an mzML mzLib can read
/// back, and that the noise statistics survive the round trip.
/// </summary>
[TestFixture]
public class NoisyMzmlExportTests
{
    private const int MinCharge = 6;
    private const int MaxCharge = 11;

    /// <summary>Fitted to the Jurkat runs; see <see cref="OrbitrapPeakWidth"/>.</summary>
    private static readonly IPeakWidthModel Width = OrbitrapPeakWidth.FromSigmaAt(800, 0.011);

    private string _outputDirectory;

    private static ProteoformModel BuildModel(double abundance = 1e9) =>
        new(
            MonoisotopicMass: 6400.0,
            Abundance: abundance,
            RtProfile: new EmgProfile(Mu: 20.0, Sigma: 0.22, Tau: 0.08),
            ChargeDistribution: new GaussianChargeDistribution(MuZ: 8.3, SigmaZ: 1.15),
            Identifier: "P00000|TESTPROTEOFORM");

    private static double[] ScanTimes(int count = 12) =>
        Enumerable.Range(0, count).Select(i => 19.5 + i * 0.1).ToArray();

    [SetUp]
    public void SetUp()
    {
        _outputDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, "TopDownSimulatorNoisyOutput");
        Directory.CreateDirectory(_outputDirectory);
    }

    [TearDown]
    public void TearDown()
    {
        if (Directory.Exists(_outputDirectory))
            Directory.Delete(_outputDirectory, recursive: true);
    }

    private string OutputPath(string name) => Path.Combine(_outputDirectory, name);

    [Test]
    public void NoisyMzmlRoundTripsThroughMzLibWithTheExpectedPeakCount()
    {
        var noise = new NoiseFloorModel(densityScale: 0.02);
        double[] scanTimes = ScanTimes();
        string path = OutputPath("noisy.mzML");

        var export = new Simulator().WriteMzml(
            new[] { BuildModel() }, MinCharge, MaxCharge, Width, scanTimes, path, noise: noise);

        Assert.That(export.Noise, Is.Not.Null);
        Assert.That(export.Noise!.NoisePeaks, Is.GreaterThan(0));
        Assert.That(export.Noise.SignalPeaks, Is.GreaterThan(0));

        var reloaded = MsDataFileReader.GetDataFile(path);
        reloaded.LoadAllStaticData();
        var scans = reloaded.GetAllScansList();

        Assert.That(scans, Has.Count.EqualTo(scanTimes.Length));
        Assert.That(scans.Sum(s => s.MassSpectrum.XArray.Length), Is.EqualTo(export.PeakCount));

        // The modelled density governs the noise contribution specifically; the total also carries
        // signal, so it is only bounded below by it.
        double expectedPerScan = noise.ExpectedPeaksPerScan;
        double noisePerScan = export.Noise.NoisePeaks / (double)export.ScanCount;
        Assert.That(noisePerScan, Is.EqualTo(expectedPerScan).Within(0.15 * expectedPerScan));

        foreach (var scan in scans)
        {
            Assert.That(scan.MassSpectrum.XArray.Length,
                Is.GreaterThan(0.5 * expectedPerScan),
                $"scan {scan.OneBasedScanNumber} carries far less than the modelled noise density");
            Assert.That(scan.MassSpectrum.XArray, Is.Ordered.Ascending);
        }
    }

    /// <summary>
    /// Every scan must carry noise, including those far outside the proteoform's elution window.
    /// The whole point of the exercise is the post-elution part of a run.
    /// </summary>
    [Test]
    public void ScansOutsideTheElutionWindowStillCarryNoise()
    {
        // The proteoform elutes at 20 min; these scans run from 40 to 41.
        double[] scanTimes = Enumerable.Range(0, 6).Select(i => 40.0 + i * 0.2).ToArray();
        string path = OutputPath("empty-region.mzML");

        var export = new Simulator().WriteMzml(
            new[] { BuildModel() }, MinCharge, MaxCharge, Width, scanTimes, path,
            noise: new NoiseFloorModel(densityScale: 0.02));

        Assert.That(export.Noise!.SignalPeaks, Is.Zero, "the proteoform should contribute nothing here");
        Assert.That(export.Noise.NoisePeaks, Is.GreaterThan(0));

        var reloaded = MsDataFileReader.GetDataFile(path);
        reloaded.LoadAllStaticData();
        foreach (var scan in reloaded.GetAllScansList())
            Assert.That(scan.MassSpectrum.XArray.Length, Is.GreaterThan(0), "an empty scan means no noise was added");
    }

    /// <summary>
    /// Supplying a noise model must switch reduction onto the noise-relative floor. Left on the
    /// default fraction-of-the-brightest-peak rule, the file would contain noise peaks fainter than
    /// any surviving signal peak, which no instrument can produce.
    /// </summary>
    [Test]
    public void NoiseModelSwitchesReductionOntoTheNoiseRelativeFloor()
    {
        var noise = new NoiseFloorModel(densityScale: 0.02);
        double[] scanTimes = ScanTimes();

        var clean = new Simulator().WriteMzml(
            new[] { BuildModel() }, MinCharge, MaxCharge, Width, scanTimes, OutputPath("clean.mzML"));
        var noisy = new Simulator().WriteMzml(
            new[] { BuildModel() }, MinCharge, MaxCharge, Width, scanTimes, OutputPath("noisy2.mzML"), noise: noise);

        // The noise-relative floor sits far below 1e-4 of the brightest peak, so more of the
        // proteoform's own faint signal survives.
        Assert.That(noisy.Noise!.SignalPeaks, Is.GreaterThan(clean.PeakCount));

        // And the faintest thing in the file is of order the noise amplitude, not of order the
        // brightest peak.
        var reloaded = MsDataFileReader.GetDataFile(OutputPath("noisy2.mzML"));
        reloaded.LoadAllStaticData();
        var all = reloaded.GetAllScansList().SelectMany(s => s.MassSpectrum.YArray).ToArray();
        double faintest = all.Min();
        double brightest = all.Max();

        Assert.That(faintest, Is.LessThan(1e-4 * brightest),
            "reduction still appears to be using the global relative floor");
        Assert.That(faintest, Is.GreaterThan(0.1 * noise.NoiseLevelAt(noise.MinMz)),
            "peaks far below the noise amplitude should not be in the file");
    }

    [Test]
    public void FeatureGroundTruthStillDescribesOnlyPeaksThatAreInTheFile()
    {
        var noise = new NoiseFloorModel(densityScale: 0.02);
        double[] scanTimes = ScanTimes();
        string path = OutputPath("truth.mzML");

        var export = new Simulator().WriteMzml(
            new[] { BuildModel() }, MinCharge, MaxCharge, Width, scanTimes, path, noise: noise);

        Assert.That(export.FeatureGroundTruthPath, Is.Not.Null);
        Assert.That(File.Exists(export.FeatureGroundTruthPath!), Is.True);
        Assert.That(export.FeatureCount, Is.GreaterThan(0));

        var reloaded = MsDataFileReader.GetDataFile(path);
        reloaded.LoadAllStaticData();
        var scansByNumber = reloaded.GetAllScansList().ToDictionary(s => s.OneBasedScanNumber);

        var header = File.ReadAllLines(export.FeatureGroundTruthPath!);
        var columns = header[0].Split('\t');
        int apexScanColumn = Array.IndexOf(columns, "ApexScanNumber");
        int apexMzColumn = Array.IndexOf(columns, "ApexMz");
        Assert.That(apexScanColumn, Is.GreaterThanOrEqualTo(0), "ApexScanNumber column is missing");
        Assert.That(apexMzColumn, Is.GreaterThanOrEqualTo(0), "ApexPeakMz column is missing");

        foreach (string line in header.Skip(1))
        {
            var fields = line.Split('\t');
            int apexScan = int.Parse(fields[apexScanColumn]);
            double apexMz = double.Parse(fields[apexMzColumn]);

            Assert.That(scansByNumber.ContainsKey(apexScan), Is.True, $"scan {apexScan} is not in the file");
            var x = scansByNumber[apexScan].MassSpectrum.XArray;
            double tolerance = 2 * Width.SigmaAt(apexMz);
            Assert.That(x.Any(mz => Math.Abs(mz - apexMz) <= tolerance), Is.True,
                $"ground truth claims a peak at {apexMz} in scan {apexScan} that is not in the file");
        }
    }

    [Test]
    public void OmittingTheNoiseModelReproducesTheCleanSimulationExactly()
    {
        double[] scanTimes = ScanTimes();

        var a = new Simulator().WriteMzml(
            new[] { BuildModel() }, MinCharge, MaxCharge, Width, scanTimes, OutputPath("a.mzML"));
        var b = new Simulator().WriteMzml(
            new[] { BuildModel() }, MinCharge, MaxCharge, Width, scanTimes, OutputPath("b.mzML"), noise: null);

        Assert.That(a.Noise, Is.Null);
        Assert.That(b.PeakCount, Is.EqualTo(a.PeakCount));
        Assert.That(b.FeatureCount, Is.EqualTo(a.FeatureCount));
    }
}
