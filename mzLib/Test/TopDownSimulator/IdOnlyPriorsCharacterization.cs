using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using NUnit.Framework;
using TopDownSimulator.Extraction;

namespace Test.TopDownSimulator;

/// <summary>
/// How well the fields an identification carries predict the parameters a fit to the raw data
/// produces. That decides what an ID-only simulation can know about each species and what it has to
/// draw from a population prior.
/// </summary>
[TestFixture]
public class IdOnlyPriorsCharacterization
{
    private const string PsmPath = @"D:\JurkatTopdown\Frac7_GPTMD_Search\Task2-TopDownSearch\Individual File Results\02-18-20_jurkat_td_rep2_fract7_Proteoforms.psmtsv";
    private const string Stem = "02-18-20_jurkat_td_rep2_fract7";
    private const string FittedPath = @"D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.full.v2.noisy.simulated.groundtruth.tsv";

    private sealed record Joined(IdentifiedSpecies Species, double Mass, double Abundance,
        double RtMu, double RtSigma, double RtTau, double ChargeMu, double ChargeSigma);

    [Test]
    [Explicit("Joins fitted rep2 fract7 models to their identifications and reports predictive power")]
    public static void CharacterizePredictors()
    {
        var species = SpeciesGrouper.Group(new MmResultLoader().LoadQualified(PsmPath, Stem, 0.01));
        var byAnchor = species.ToDictionary(s => s.Anchor.Identifier);

        var lines = File.ReadAllLines(FittedPath);
        var h = lines[0].Split('\t');
        int C(string n) => Array.IndexOf(h, n);
        var joined = new List<Joined>();
        foreach (var line in lines.Skip(1))
        {
            var f = line.Split('\t');
            if (!byAnchor.TryGetValue(f[C("Identifier")], out var s)) continue;
            double D(string n) => double.Parse(f[C(n)], CultureInfo.InvariantCulture);
            joined.Add(new Joined(s, D("MonoisotopicMass"), D("Abundance"), D("RtMu"), D("RtSigma"), D("RtTau"),
                D("ChargeMu"), D("ChargeSigma")));
        }

        Console.WriteLine($"species {species.Length}, fitted {lines.Length - 1}, joined {joined.Count}");
        var fit = joined.Where(j => j.Abundance > 0).ToArray();
        Console.WriteLine($"with positive abundance after refit: {fit.Length}");

        // --- Retention time apex ---
        Console.WriteLine();
        Console.WriteLine("RT: fitted mu minus ...");
        Report("anchor MS2 RT", fit.Select(j => j.RtMu - j.Species.Anchor.RetentionTime));
        Report("median member RT", fit.Select(j => j.RtMu - Median(j.Species.Members.Select(m => m.RetentionTime))));
        Report("intensity-weighted member RT", fit.Select(j => j.RtMu - Weighted(j.Species.Members)));

        // --- Abundance ---
        Console.WriteLine();
        Console.WriteLine("log10 abundance regressed on:");
        var logA = fit.Select(j => Math.Log10(j.Abundance)).ToArray();
        Regress("log10 anchor precursor intensity", fit.Select(j => Log(j.Species.Anchor.PrecursorIntensity)), logA);
        Regress("log10 max member precursor intensity", fit.Select(j => Log(j.Species.Members.Max(m => m.PrecursorIntensity ?? 0))), logA);
        Regress("log10 summed member precursor intensity", fit.Select(j => Log(j.Species.Members.Sum(m => m.PrecursorIntensity ?? 0))), logA);
        Regress("log10 member count", fit.Select(j => Math.Log10(j.Species.Members.Length)), logA);
        Report("log10 abundance (prior)", logA);

        // --- Charge ---
        Console.WriteLine();
        Console.WriteLine("charge mu regressed on:");
        var muZ = fit.Select(j => j.ChargeMu).ToArray();
        Regress("mass (kDa)", fit.Select(j => j.Mass / 1000), muZ);
        Regress("sqrt mass", fit.Select(j => Math.Sqrt(j.Mass)), muZ);
        Regress("mean member precursor charge", fit.Select(j => j.Species.Members.Average(m => (double)m.PrecursorCharge)), muZ);
        Regress("anchor precursor charge", fit.Select(j => (double)j.Species.Anchor.PrecursorCharge), muZ);
        Regress("basic residues K+R", fit.Select(j => (double)Basic(j, "KR")), muZ);
        Regress("basic residues K+R+H", fit.Select(j => (double)Basic(j, "KRH")), muZ);
        Regress("basic residues K+R+H + N-terminus", fit.Select(j => Basic(j, "KRH") + 1.0), muZ);
        Regress2("mean member charge + K+R",
            fit.Select(j => j.Species.Members.Average(m => (double)m.PrecursorCharge)).ToArray(),
            fit.Select(j => (double)Basic(j, "KR")).ToArray(), muZ);
        Regress2("mass (kDa) + K+R",
            fit.Select(j => j.Mass / 1000).ToArray(), fit.Select(j => (double)Basic(j, "KR")).ToArray(), muZ);
        Regress("charge sigma ~ K+R", fit.Select(j => (double)Basic(j, "KR")), fit.Select(j => j.ChargeSigma).ToArray());

        // --- Shape priors ---
        Console.WriteLine();
        Report("rt sigma", fit.Select(j => j.RtSigma));
        Report("rt tau", fit.Select(j => j.RtTau));
        Console.WriteLine($"  tau == 0 fraction: {fit.Count(j => j.RtTau <= 1e-9) / (double)fit.Length:P1}");
        Report("charge sigma", fit.Select(j => j.ChargeSigma));
        Regress("rt sigma ~ rt mu", fit.Select(j => j.RtMu), fit.Select(j => j.RtSigma).ToArray());
        Regress("charge sigma ~ mass (kDa)", fit.Select(j => j.Mass / 1000), fit.Select(j => j.ChargeSigma).ToArray());
        Report("member count", fit.Select(j => (double)j.Species.Members.Length));
    }

    /// <summary>Residues of <paramref name="residues"/> in the anchor's sequence, modifications stripped.</summary>
    private static int Basic(Joined j, string residues) =>
        System.Text.RegularExpressions.Regex.Replace(j.Species.Anchor.FullSequence, @"\[[^\]]*\]", "")
            .Count(residues.Contains);

    /// <summary>Two-predictor least squares, solved through the 2x2 normal equations on centred data.</summary>
    private static void Regress2(string label, double[] x1, double[] x2, double[] y)
    {
        int n = y.Length;
        double m1 = x1.Average(), m2 = x2.Average(), my = y.Average();
        double s11 = 0, s22 = 0, s12 = 0, s1y = 0, s2y = 0, syy = 0;
        for (int i = 0; i < n; i++)
        {
            double a = x1[i] - m1, b = x2[i] - m2, c = y[i] - my;
            s11 += a * a; s22 += b * b; s12 += a * b; s1y += a * c; s2y += b * c; syy += c * c;
        }

        double det = s11 * s22 - s12 * s12;
        double b1 = (s22 * s1y - s12 * s2y) / det, b2 = (s11 * s2y - s12 * s1y) / det;
        double a0 = my - b1 * m1 - b2 * m2;
        double rss = 0;
        for (int i = 0; i < n; i++) rss += Math.Pow(y[i] - (a0 + b1 * x1[i] + b2 * x2[i]), 2);
        Console.WriteLine($"  {label,-42} n={n,4}  y = {a0,8:F3} + {b1,7:F4} x1 + {b2,7:F4} x2   R2 {1 - rss / syy,6:F3}   resid sd {Math.Sqrt(rss / (n - 3)),7:F3}");
    }

    private static double Log(double? v) => v is > 0 ? Math.Log10(v.Value) : double.NaN;

    private static double Weighted(MmResultRecord[] members)
    {
        double w = members.Sum(m => m.PrecursorIntensity ?? 0);
        return w > 0 ? members.Sum(m => (m.PrecursorIntensity ?? 0) * m.RetentionTime) / w : members[0].RetentionTime;
    }

    private static void Regress(string label, IEnumerable<double> xs, double[] ys)
    {
        var pairs = xs.Zip(ys).Where(p => double.IsFinite(p.First) && double.IsFinite(p.Second)).ToArray();
        double mx = pairs.Average(p => p.First), my = pairs.Average(p => p.Second);
        double sxx = pairs.Sum(p => (p.First - mx) * (p.First - mx));
        double sxy = pairs.Sum(p => (p.First - mx) * (p.Second - my));
        double syy = pairs.Sum(p => (p.Second - my) * (p.Second - my));
        double b = sxy / sxx, a = my - b * mx;
        double rss = pairs.Sum(p => Math.Pow(p.Second - (a + b * p.First), 2));
        Console.WriteLine($"  {label,-42} n={pairs.Length,4}  y = {a,8:F3} + {b,7:F4} x   R2 {1 - rss / syy,6:F3}   resid sd {Math.Sqrt(rss / (pairs.Length - 2)),7:F3}   (y sd {Math.Sqrt(syy / (pairs.Length - 1)):F3})");
    }

    private static void Report(string label, IEnumerable<double> values)
    {
        var v = values.Where(double.IsFinite).OrderBy(x => x).ToArray();
        double Q(double q) => v[(int)Math.Min(v.Length - 1, Math.Floor(q * v.Length))];
        Console.WriteLine($"  {label,-32} n={v.Length,4}  p5 {Q(0.05),8:F3}  p25 {Q(0.25),8:F3}  p50 {Q(0.5),8:F3}  p75 {Q(0.75),8:F3}  p95 {Q(0.95),8:F3}  mean {v.Average(),8:F3}");
    }

    private static double Median(IEnumerable<double> values)
    {
        var v = values.OrderBy(x => x).ToArray();
        return v.Length % 2 == 1 ? v[v.Length / 2] : 0.5 * (v[v.Length / 2 - 1] + v[v.Length / 2]);
    }
}
