using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.Json;
using MassSpectrometry;
using NUnit.Framework;
using TopDownSimulator.Model;
using TopDownSimulator.Simulation;

namespace Test.TopDownSimulator;

/// <summary>
/// Generates the data-free SimSpec v1 golden fixture that the Rust sim-core port parity-tests
/// against (see ProteoStar agent_info/TopDown-SimSpec-Contract-v1.md §5). Signal only — no noise,
/// no reduction — so the fixture is reproducible offline and isolates the forward model +
/// centroider. [Explicit]: run deliberately to (re)write the committed fixture, not on every build.
/// </summary>
[TestFixture]
[Explicit("Writes committed golden fixtures; run on demand, not in CI.")]
public class SimSpecGoldenGenerator
{
    // One small proteoform (monoisotopic peak present) and one above ~31 kDa (monoisotopic peak
    // dropped by averagine's 1e-8 cutoff) so the Rust port must exercise the mass-anchored shift.
    private static ProteoformModel[] GoldenModels() => new[]
    {
        new ProteoformModel(5000.25, 1.0e9,
            new EmgProfile(20.0, 0.22, 0.08),
            new GaussianChargeDistribution(7.5, 1.2), "P1|SMALL"),
        new ProteoformModel(35000.75, 4.2e8,
            new EmgProfile(20.6, 0.30, 0.12),
            new GaussianChargeDistribution(30.0, 3.5), "P2|LARGE"),
    };

    private const int MinCharge = 3;
    private const int MaxCharge = 40;
    private static readonly OrbitrapPeakWidth Width = new(1.2e-7);

    private static double[] GoldenScanTimes() =>
        Enumerable.Range(0, 20).Select(i => 19.5 + i * 0.1).ToArray();

    [Test]
    public void WriteGoldenV1()
    {
        var models = GoldenModels();
        double[] scanTimes = GoldenScanTimes();

        string dir = FixtureDir();
        Directory.CreateDirectory(dir);

        // 1. The spec (signal only: noise/reduction null).
        var spec = SimSpec.FromModels(
            models, MinCharge, MaxCharge, Width, scanTimes,
            generator: SimGenerator.Now(),
            sourceRun: new SimSourceRun(Description: "hand-built data-free golden; no raw source"));
        spec.Write(Path.Combine(dir, "golden-v1.simspec.json"));

        // 2. The expected centroided output of exactly what the spec describes.
        var inputs = spec.ToSimulationInputs();
        var simulation = new Simulator().SimulateCentroided(
            inputs.Proteoforms, inputs.MinCharge, inputs.MaxCharge, inputs.WidthModel, inputs.ScanTimes);
        WriteScansTsv(Path.Combine(dir, "golden-v1.scans.tsv"), simulation);

        // 3. The averagine table the envelope kernel reads, so Rust can match envelope shapes.
        WriteAveragine(Path.Combine(dir, "averagine.v1.json"));

        long peakCount = simulation.Scans.Sum(s => (long)s.MassSpectrum.XArray.Length);
        Assert.That(simulation.Scans, Has.Length.EqualTo(scanTimes.Length));
        Assert.That(peakCount, Is.GreaterThan(0), "Golden produced no peaks.");
        TestContext.Out.WriteLine($"Wrote golden-v1 to {dir}: {simulation.Scans.Length} scans, {peakCount} peaks.");
    }

    private static void WriteScansTsv(string path, SimulationResult simulation)
    {
        var sb = new StringBuilder();
        sb.AppendLine(string.Join('\t', "scanIndex", "oneBasedScanNumber", "retentionTime", "mz", "intensity"));
        for (int s = 0; s < simulation.Scans.Length; s++)
        {
            var scan = simulation.Scans[s];
            double[] mz = scan.MassSpectrum.XArray;
            double[] intensity = scan.MassSpectrum.YArray;
            for (int i = 0; i < mz.Length; i++)
                sb.AppendLine(string.Join('\t',
                    s.ToString(CultureInfo.InvariantCulture),
                    scan.OneBasedScanNumber.ToString(CultureInfo.InvariantCulture),
                    R(scan.RetentionTime), R(mz[i]), R(intensity[i])));
        }
        File.WriteAllText(path, sb.ToString());
    }

    private static void WriteAveragine(string path)
    {
        var table = new
        {
            allMasses = Averagine.AllMasses,
            allIntensities = Averagine.AllIntensities,
            mostIntenseMasses = Averagine.MostIntenseMasses,
            diffToMonoisotopic = Averagine.DiffToMonoisotopic,
        };
        File.WriteAllText(path, JsonSerializer.Serialize(table, new JsonSerializerOptions { WriteIndented = false }));
    }

    private static string R(double v) => v.ToString("R", CultureInfo.InvariantCulture);

    /// <summary>Walks up from the test output dir to the Test project, then into SimSpecFixtures.</summary>
    private static string FixtureDir()
    {
        var d = new DirectoryInfo(TestContext.CurrentContext.TestDirectory);
        while (d is not null && !File.Exists(Path.Combine(d.FullName, "Test.csproj")))
            d = d.Parent;
        if (d is null)
            throw new DirectoryNotFoundException("Could not locate the Test project directory (no Test.csproj found above the test output dir).");
        return Path.Combine(d.FullName, "TopDownSimulator", "SimSpecFixtures");
    }
}
