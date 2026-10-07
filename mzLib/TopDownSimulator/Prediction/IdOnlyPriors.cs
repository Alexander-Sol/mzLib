#nullable enable
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.Json;
using TopDownSimulator.Extraction;
using TopDownSimulator.Model;

namespace TopDownSimulator.Prediction;

/// <summary>
/// One fitted species' shape parameters, kept together so that draws preserve their correlation.
/// </summary>
public sealed record ShapeSample(double RtSigma, double RtTau, double ChargeSigma);

/// <summary>
/// Everything needed to turn identifications into <see cref="ProteoformModel"/>s without the raw
/// data they came from, learned from a run where models were fitted to the raw data.
/// </summary>
/// <remarks>
/// <para>
/// What an identification can predict was measured on Jurkat rep2 fract7, 451 species fitted from
/// the raw data:
/// </para>
/// <list type="bullet">
/// <item><b>Elution apex:</b> the anchor's MS2 time, median offset −0.014 min, interquartile range
/// 0.12 min, about one elution σ.</item>
/// <item><b>Abundance:</b> log10 of the summed precursor intensity of the species' members,
/// R² 0.74, residual 0.34 dex.</item>
/// <item><b>Charge-distribution centre:</b> the mean member precursor charge plus the anchor's count
/// of lysines and arginines, R² 0.946, residual 0.66 charges. The precursor charges carry almost all
/// of it (0.939 alone); mass alone manages 0.48 and K+R alone 0.30.</item>
/// <item><b>Elution width, tail and charge-distribution width:</b> no measurable dependence on
/// retention time or mass, so they are drawn jointly from the fitted population.</item>
/// </list>
/// <para>
/// Predictions are regressions, not draws, wherever there is a predictor. The regression residual
/// is the error of not having the raw data, not real variation between species, so adding it back
/// would make individual species worse without making the population any more realistic.
/// </para>
/// </remarks>
public sealed record IdOnlyPriors(
    double RtApexOffset,
    double LogAbundanceIntercept,
    double LogAbundanceSlope,
    double FallbackLogAbundance,
    double ChargeMuIntercept,
    double ChargeMuChargeSlope,
    double ChargeMuBasicSlope,
    ShapeSample[] Shapes,
    string Source)
{
    /// <summary>
    /// Learns the priors from species whose models were fitted to the raw data.
    /// </summary>
    /// <param name="fitted">Fitted models, matched to <paramref name="species"/> by the anchor's identifier.</param>
    /// <param name="source">A description of the training run, carried into the saved file.</param>
    public static IdOnlyPriors Fit(
        IReadOnlyList<IdentifiedSpecies> species,
        IReadOnlyList<ProteoformModel> fitted,
        string source)
    {
        if (species is null) throw new ArgumentNullException(nameof(species));
        if (fitted is null) throw new ArgumentNullException(nameof(fitted));

        var byAnchor = species.ToDictionary(s => s.Anchor.Identifier);
        var pairs = fitted
            .Where(m => m.Identifier is not null && m.Abundance > 0 && byAnchor.ContainsKey(m.Identifier))
            .Select(m => (Species: byAnchor[m.Identifier!], Model: m))
            .ToArray();
        if (pairs.Length < 10)
            throw new InvalidOperationException(
                $"Only {pairs.Length} fitted models could be matched to species; at least 10 are needed.");

        double rtOffset = Median(pairs.Select(p => p.Model.RtProfile.Mu - p.Species.Anchor.RetentionTime));

        var abundance = pairs
            .Select(p => (X: LogSummedPrecursorIntensity(p.Species), Y: Math.Log10(p.Model.Abundance)))
            .Where(p => double.IsFinite(p.X))
            .ToArray();
        var (aIntercept, aSlope) = LeastSquares(abundance);
        double fallback = Median(pairs.Select(p => Math.Log10(p.Model.Abundance)));

        var charge = pairs
            .Where(p => p.Model.ChargeDistribution is GaussianChargeDistribution)
            .Select(p => (X1: MeanPrecursorCharge(p.Species), X2: (double)BasicResidueCount(p.Species),
                Y: ((GaussianChargeDistribution)p.Model.ChargeDistribution).MuZ))
            .ToArray();
        var (zIntercept, zChargeSlope, zBasicSlope) = LeastSquares(charge);

        var shapes = pairs
            .Where(p => p.Model.ChargeDistribution is GaussianChargeDistribution)
            .Select(p => new ShapeSample(
                p.Model.RtProfile.Sigma,
                p.Model.RtProfile.Tau,
                ((GaussianChargeDistribution)p.Model.ChargeDistribution).SigmaZ))
            .ToArray();

        return new IdOnlyPriors(rtOffset, aIntercept, aSlope, fallback, zIntercept, zChargeSlope, zBasicSlope, shapes, source);
    }

    /// <summary>
    /// A model for every species. Shape parameters are drawn per species from a stream seeded by
    /// <paramref name="seed"/> and the species' anchor identifier, so a species gets the same draw
    /// whatever else is in the set.
    /// </summary>
    public ProteoformModel[] Predict(IReadOnlyList<IdentifiedSpecies> species, int seed = 0)
    {
        if (species is null) throw new ArgumentNullException(nameof(species));
        if (Shapes.Length == 0) throw new InvalidOperationException("The priors carry no shape samples.");

        return species.Select(s =>
        {
            double logIntensity = LogSummedPrecursorIntensity(s);
            double logAbundance = double.IsFinite(logIntensity)
                ? LogAbundanceIntercept + LogAbundanceSlope * logIntensity
                : FallbackLogAbundance;

            var shape = Shapes[StableIndex(seed, s.Anchor.Identifier, Shapes.Length)];
            double muZ = ChargeMuIntercept + ChargeMuChargeSlope * MeanPrecursorCharge(s) + ChargeMuBasicSlope * BasicResidueCount(s);

            return new ProteoformModel(
                s.Anchor.MonoisotopicMass,
                Math.Pow(10, logAbundance),
                new EmgProfile(s.Anchor.RetentionTime + RtApexOffset, shape.RtSigma, shape.RtTau),
                new GaussianChargeDistribution(muZ, shape.ChargeSigma),
                s.Anchor.Identifier);
        }).ToArray();
    }

    public void Save(string path) =>
        File.WriteAllText(path, JsonSerializer.Serialize(this, new JsonSerializerOptions { WriteIndented = true }));

    public static IdOnlyPriors Load(string path) =>
        JsonSerializer.Deserialize<IdOnlyPriors>(File.ReadAllText(path))
        ?? throw new InvalidDataException($"{path} does not hold ID-only priors.");

    /// <summary>log10 of the members' summed precursor intensity, or NaN when none reported one.</summary>
    public static double LogSummedPrecursorIntensity(IdentifiedSpecies species)
    {
        double sum = species.Members.Sum(m => m.PrecursorIntensity is > 0 ? m.PrecursorIntensity.Value : 0);
        return sum > 0 ? Math.Log10(sum) : double.NaN;
    }

    public static double MeanPrecursorCharge(IdentifiedSpecies species) =>
        species.Members.Average(m => (double)m.PrecursorCharge);

    /// <summary>
    /// Lysines and arginines in the anchor's sequence, modifications stripped. With the observed
    /// precursor charges it explains the fitted charge centre slightly better than either alone
    /// (R² 0.946 against 0.939 on rep2 fract7). On its own it is weaker than mass (0.30 against 0.48).
    /// </summary>
    public static int BasicResidueCount(IdentifiedSpecies species)
    {
        int count = 0, depth = 0;
        foreach (char c in species.Anchor.FullSequence)
        {
            if (c == '[') depth++;
            else if (c == ']') depth--;
            else if (depth == 0 && (c == 'K' || c == 'R')) count++;
        }

        return count;
    }

    private static int StableIndex(int seed, string identifier, int count)
    {
        // FNV-1a, so the draw does not depend on string.GetHashCode's per-process randomisation.
        unchecked
        {
            uint hash = 2166136261u ^ (uint)seed;
            foreach (char c in identifier)
                hash = (hash ^ c) * 16777619u;
            return (int)(hash % (uint)count);
        }
    }

    /// <summary>Two-predictor least squares through the normal equations on centred data.</summary>
    private static (double Intercept, double Slope1, double Slope2) LeastSquares(IReadOnlyList<(double X1, double X2, double Y)> points)
    {
        if (points.Count < 3) throw new InvalidOperationException("At least three points are needed for a two-predictor regression.");
        double m1 = points.Average(p => p.X1), m2 = points.Average(p => p.X2), my = points.Average(p => p.Y);
        double s11 = 0, s22 = 0, s12 = 0, s1y = 0, s2y = 0;
        foreach (var p in points)
        {
            double a = p.X1 - m1, b = p.X2 - m2, c = p.Y - my;
            s11 += a * a; s22 += b * b; s12 += a * b; s1y += a * c; s2y += b * c;
        }

        double det = s11 * s22 - s12 * s12;
        if (!(Math.Abs(det) > 1e-12 * Math.Max(1, s11 * s22)))
        {
            var (i, s) = LeastSquares(points.Select(p => (p.X1, p.Y)).ToArray());
            return (i, s, 0);
        }

        double b1 = (s22 * s1y - s12 * s2y) / det, b2 = (s11 * s2y - s12 * s1y) / det;
        return (my - b1 * m1 - b2 * m2, b1, b2);
    }

    private static (double Intercept, double Slope) LeastSquares(IReadOnlyList<(double X, double Y)> points)
    {
        if (points.Count < 2) throw new InvalidOperationException("At least two points are needed for a regression.");
        double mx = points.Average(p => p.X), my = points.Average(p => p.Y);
        double sxx = points.Sum(p => (p.X - mx) * (p.X - mx));
        double sxy = points.Sum(p => (p.X - mx) * (p.Y - my));
        if (!(sxx > 0)) return (my, 0);
        double slope = sxy / sxx;
        return (my - slope * mx, slope);
    }

    private static double Median(IEnumerable<double> values)
    {
        var v = values.Where(double.IsFinite).OrderBy(x => x).ToArray();
        if (v.Length == 0) throw new InvalidOperationException("No finite values to take a median of.");
        return v.Length % 2 == 1 ? v[v.Length / 2] : 0.5 * (v[v.Length / 2 - 1] + v[v.Length / 2]);
    }
}
