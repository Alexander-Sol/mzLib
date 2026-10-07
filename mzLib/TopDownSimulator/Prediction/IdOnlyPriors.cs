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

/// <summary>Where <see cref="IdOnlyPriors.Predict"/> takes the charge-distribution centre from.</summary>
public enum ChargePredictor
{
    /// <summary>
    /// Mean precursor charge of the species' identifications, plus K+R. R² 0.45 on rep2 fract7 fits
    /// over the whole charge range; biased low, because MS2 selects charges below the apex.
    /// </summary>
    ObservedCharges,

    /// <summary>√mass and net basic residues (K+R+H − D−E). R² 0.70.</summary>
    SequenceOnly,

    /// <summary>
    /// √mass and a separate weight for each charged residue (<see cref="ResidueChargeModel"/>). R² 0.72,
    /// the best found. Falls back to <see cref="SequenceOnly"/> when the priors carry no residue model.
    /// </summary>
    ResidueWeights,
}

/// <summary>
/// μ_z = intercept + slope·√mass + Σ weight_r·count_r over the residues in <see cref="Residues"/>,
/// counted in the anchor's sequence with modifications stripped.
/// </summary>
/// <param name="ResidualSd">Residual standard deviation on the training species, in charges.</param>
public sealed record ResidueChargeModel(
    double Intercept,
    double SqrtMassSlope,
    string Residues,
    double[] Weights,
    double ResidualSd)
{
    /// <summary>The residues the default model weighs: the basic K, R and H and the acidic D and E.</summary>
    public const string ChargedResidues = "KRHDE";

    public double Predict(IdentifiedSpecies species)
    {
        double mu = Intercept + SqrtMassSlope * Math.Sqrt(species.Anchor.MonoisotopicMass);
        for (int i = 0; i < Residues.Length; i++)
            mu += Weights[i] * IdOnlyPriors.CountResidues(species, Residues[i].ToString());
        return mu;
    }

    /// <summary>Ordinary least squares of <paramref name="targets"/> on √mass and each residue count.</summary>
    public static ResidueChargeModel Fit(
        IReadOnlyList<IdentifiedSpecies> species, IReadOnlyList<double> targets, string residues = ChargedResidues)
    {
        if (species.Count != targets.Count) throw new ArgumentException("One target per species is needed.");
        int p = residues.Length + 2;
        if (species.Count <= p) throw new InvalidOperationException($"At least {p + 1} species are needed.");

        var x = species.Select(s =>
        {
            var row = new double[p];
            row[0] = 1;
            row[1] = Math.Sqrt(s.Anchor.MonoisotopicMass);
            for (int i = 0; i < residues.Length; i++)
                row[i + 2] = IdOnlyPriors.CountResidues(s, residues[i].ToString());
            return row;
        }).ToArray();

        double[] beta = LinearLeastSquares.Solve(x, targets);
        double rss = 0;
        for (int n = 0; n < x.Length; n++)
        {
            double fit = 0;
            for (int j = 0; j < p; j++) fit += beta[j] * x[n][j];
            rss += (targets[n] - fit) * (targets[n] - fit);
        }

        return new ResidueChargeModel(beta[0], beta[1], residues, beta.Skip(2).ToArray(), Math.Sqrt(rss / (x.Length - p)));
    }
}

/// <summary>Least squares through the normal equations, solved by Gaussian elimination with partial pivoting.</summary>
public static class LinearLeastSquares
{
    public static double[] Solve(IReadOnlyList<double[]> x, IReadOnlyList<double> y)
    {
        int p = x[0].Length;
        var a = new double[p, p + 1];
        for (int n = 0; n < x.Count; n++)
        for (int i = 0; i < p; i++)
        {
            for (int j = 0; j < p; j++) a[i, j] += x[n][i] * x[n][j];
            a[i, p] += x[n][i] * y[n];
        }

        double scale = 0;
        for (int i = 0; i < p; i++) scale = Math.Max(scale, Math.Abs(a[i, i]));

        for (int col = 0; col < p; col++)
        {
            int pivot = col;
            for (int r = col + 1; r < p; r++)
                if (Math.Abs(a[r, col]) > Math.Abs(a[pivot, col])) pivot = r;
            if (!(Math.Abs(a[pivot, col]) > 1e-10 * scale))
                throw new InvalidOperationException("The predictors are collinear.");
            for (int j = col; j <= p; j++) (a[col, j], a[pivot, j]) = (a[pivot, j], a[col, j]);
            for (int r = 0; r < p; r++)
            {
                if (r == col) continue;
                double f = a[r, col] / a[col, col];
                for (int j = col; j <= p; j++) a[r, j] -= f * a[col, j];
            }
        }

        var beta = new double[p];
        for (int i = 0; i < p; i++) beta[i] = a[i, p] / a[i, i];
        return beta;
    }
}

/// <summary>
/// Everything needed to turn identifications into <see cref="ProteoformModel"/>s without the raw
/// data they came from, learned from a run where models were fitted to the raw data.
/// </summary>
/// <remarks>
/// <para>
/// What an identification can predict was measured on Jurkat rep2 fract7, 436 species fitted from
/// the raw data over each species' whole observable charge range (the <c>.v3</c> fits):
/// </para>
/// <list type="bullet">
/// <item><b>Elution apex:</b> the anchor's MS2 time, median offset −0.015 min, interquartile range
/// 0.12 min, about one elution σ.</item>
/// <item><b>Abundance:</b> log10 of the brightest member's precursor intensity, R² 0.59, residual
/// 0.44 dex. Priors fitted before this carry <see cref="AbundanceFromMaxIntensity"/> false and use
/// the members' summed intensity.</item>
/// <item><b>Charge-distribution centre:</b> √mass and a weight per charged residue
/// (<see cref="ResidueChargeModel"/>), R² 0.72, residual 1.57 charges. The precursor charges add
/// nothing to it and alone manage 0.36: MS2 selects charges below the envelope's apex. (Fits over
/// the members' charges ± 2 only made them look like a near-perfect predictor.)</item>
/// <item><b>Charge-distribution width:</b> σ_z grows with μ_z (R² 0.42), so it is predicted from
/// the predicted μ_z.</item>
/// <item><b>Elution width and tail:</b> no measurable dependence on retention time or mass, so
/// they are drawn jointly from the fitted population.</item>
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
    string Source,
    double SequenceChargeIntercept = 0,
    double SequenceChargeSqrtMassSlope = 0,
    double SequenceChargeNetBasicSlope = 0,
    ResidueChargeModel? ResidueCharge = null,
    double? ChargeSigmaIntercept = null,
    double ChargeSigmaMuSlope = 0,
    bool AbundanceFromMaxIntensity = false)
{
    /// <summary>The narrowest charge distribution <see cref="Predict"/> will produce from the σ_z regression.</summary>
    public const double MinimumPredictedChargeSigma = 0.5;
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

        // The brightest member's precursor intensity, not the members' sum: R² 0.59 against 0.53 for
        // models fitted over the scan's whole charge range on rep2 fract7.
        var abundance = pairs
            .Select(p => (X: LogMaxPrecursorIntensity(p.Species), Y: Math.Log10(p.Model.Abundance)))
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

        var sequenceCharge = pairs
            .Where(p => p.Model.ChargeDistribution is GaussianChargeDistribution)
            .Select(p => (X1: Math.Sqrt(p.Species.Anchor.MonoisotopicMass), X2: (double)NetBasicResidueCount(p.Species),
                Y: ((GaussianChargeDistribution)p.Model.ChargeDistribution).MuZ))
            .ToArray();
        var (sIntercept, sMassSlope, sNetBasicSlope) = LeastSquares(sequenceCharge);

        // Small or uniform training sets cannot separate the residues; the net-basic model then stands alone.
        var gaussian = pairs.Where(p => p.Model.ChargeDistribution is GaussianChargeDistribution).ToArray();
        ResidueChargeModel? residueCharge;
        try
        {
            residueCharge = ResidueChargeModel.Fit(
                gaussian.Select(p => p.Species).ToArray(),
                gaussian.Select(p => ((GaussianChargeDistribution)p.Model.ChargeDistribution).MuZ).ToArray());
        }
        catch (InvalidOperationException)
        {
            residueCharge = null;
        }

        // Broader envelopes at higher charge: σ_z on μ_z, R² 0.42 on rep2 fract7 fitted over the scan's
        // whole charge range.
        var (sigmaIntercept, sigmaSlope) = LeastSquares(gaussian
            .Select(p => (GaussianChargeDistribution)p.Model.ChargeDistribution)
            .Select(c => (X: c.MuZ, Y: c.SigmaZ))
            .ToArray());

        var shapes = pairs
            .Where(p => p.Model.ChargeDistribution is GaussianChargeDistribution)
            .Select(p => new ShapeSample(
                p.Model.RtProfile.Sigma,
                p.Model.RtProfile.Tau,
                ((GaussianChargeDistribution)p.Model.ChargeDistribution).SigmaZ))
            .ToArray();

        return new IdOnlyPriors(rtOffset, aIntercept, aSlope, fallback, zIntercept, zChargeSlope, zBasicSlope, shapes, source,
            sIntercept, sMassSlope, sNetBasicSlope, residueCharge, sigmaIntercept, sigmaSlope, AbundanceFromMaxIntensity: true);
    }

    /// <summary>
    /// σ_z for a distribution centred on <paramref name="muZ"/>: the regression on μ_z when the priors
    /// carry one, otherwise <paramref name="drawn"/>, the bootstrap draw.
    /// </summary>
    public double ChargeSigmaFor(double muZ, double drawn) =>
        ChargeSigmaIntercept is { } a ? Math.Max(MinimumPredictedChargeSigma, a + ChargeSigmaMuSlope * muZ) : drawn;

    /// <summary>
    /// The charge-distribution centre predicted from mass and sequence alone: the residue-weighted
    /// model when the priors carry one, otherwise √mass and net basic residues.
    /// </summary>
    public double SequenceChargeMu(IdentifiedSpecies species, bool useResidueWeights = true) =>
        useResidueWeights && ResidueCharge is not null
            ? ResidueCharge.Predict(species)
            : SequenceChargeIntercept + SequenceChargeSqrtMassSlope * Math.Sqrt(species.Anchor.MonoisotopicMass)
              + SequenceChargeNetBasicSlope * NetBasicResidueCount(species);

    /// <summary>
    /// A charge range that should hold the species' envelope without knowing its precursor charges:
    /// the sequence-predicted centre ± <paramref name="sigmas"/> times the population's median σ_z,
    /// widened in quadrature by the predictor's residual when the residue model is present.
    /// </summary>
    public (int Min, int Max) SequenceChargeRange(IdentifiedSpecies species, double sigmas = 3, int minCharge = 2, int maxCharge = 80)
    {
        double sigmaZ = Median(Shapes.Select(s => s.ChargeSigma));
        double residual = ResidueCharge?.ResidualSd ?? 0;
        double halfWidth = sigmas * Math.Sqrt(sigmaZ * sigmaZ + residual * residual);
        double mu = SequenceChargeMu(species);
        return (Math.Max(minCharge, (int)Math.Floor(mu - halfWidth)), Math.Min(maxCharge, (int)Math.Ceiling(mu + halfWidth)));
    }

    /// <summary>
    /// A model for every species. Shape parameters are drawn per species from a stream seeded by
    /// <paramref name="seed"/> and the species' anchor identifier, so a species gets the same draw
    /// whatever else is in the set.
    /// </summary>
    /// <param name="chargePredictor">
    /// Whether the charge centre comes from the identifications' precursor charges or, for input that
    /// does not carry them, from mass and sequence alone.
    /// </param>
    public ProteoformModel[] Predict(
        IReadOnlyList<IdentifiedSpecies> species, int seed = 0, ChargePredictor chargePredictor = ChargePredictor.ObservedCharges)
    {
        if (species is null) throw new ArgumentNullException(nameof(species));
        if (Shapes.Length == 0) throw new InvalidOperationException("The priors carry no shape samples.");

        return species.Select(s =>
        {
            double logIntensity = AbundanceFromMaxIntensity ? LogMaxPrecursorIntensity(s) : LogSummedPrecursorIntensity(s);
            double logAbundance = double.IsFinite(logIntensity)
                ? LogAbundanceIntercept + LogAbundanceSlope * logIntensity
                : FallbackLogAbundance;

            var shape = Shapes[StableIndex(seed, s.Anchor.Identifier, Shapes.Length)];
            double muZ = chargePredictor switch
            {
                ChargePredictor.ObservedCharges =>
                    ChargeMuIntercept + ChargeMuChargeSlope * MeanPrecursorCharge(s) + ChargeMuBasicSlope * BasicResidueCount(s),
                ChargePredictor.ResidueWeights => SequenceChargeMu(s),
                _ => SequenceChargeMu(s, useResidueWeights: false),
            };

            return new ProteoformModel(
                s.Anchor.MonoisotopicMass,
                Math.Pow(10, logAbundance),
                new EmgProfile(s.Anchor.RetentionTime + RtApexOffset, shape.RtSigma, shape.RtTau),
                new GaussianChargeDistribution(muZ, ChargeSigmaFor(muZ, shape.ChargeSigma)),
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

    /// <summary>log10 of the brightest member's precursor intensity, or NaN when none reported one.</summary>
    public static double LogMaxPrecursorIntensity(IdentifiedSpecies species)
    {
        double max = species.Members.Max(m => m.PrecursorIntensity is > 0 ? m.PrecursorIntensity.Value : 0);
        return max > 0 ? Math.Log10(max) : double.NaN;
    }

    public static double MeanPrecursorCharge(IdentifiedSpecies species) =>
        species.Members.Average(m => (double)m.PrecursorCharge);

    /// <summary>
    /// Lysines and arginines in the anchor's sequence, modifications stripped. With the observed
    /// precursor charges it explains the fitted charge centre slightly better than either alone
    /// (R² 0.946 against 0.939 on rep2 fract7). On its own it is weaker than mass (0.30 against 0.48).
    /// </summary>
    public static int BasicResidueCount(IdentifiedSpecies species) => CountResidues(species, "KR");

    /// <summary>
    /// K+R+H minus D+E in the anchor's sequence, modifications stripped. With √mass it is the best
    /// charge-centre predictor that needs no precursor charge: R² 0.62 on rep2 fract7, against 0.50
    /// for √mass alone and 0.59 for √mass with K+R.
    /// </summary>
    public static int NetBasicResidueCount(IdentifiedSpecies species) =>
        CountResidues(species, "KRH") - CountResidues(species, "DE");

    internal static int CountResidues(IdentifiedSpecies species, string residues)
    {
        int count = 0, depth = 0;
        foreach (char c in species.Anchor.FullSequence)
        {
            if (c == '[') depth++;
            else if (c == ']') depth--;
            else if (depth == 0 && residues.IndexOf(c) >= 0) count++;
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
