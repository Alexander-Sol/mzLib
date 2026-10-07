#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;

namespace TopDownSimulator.Extraction;

/// <summary>
/// The identifications that describe one MS1 species, which is simulated as a single model.
/// </summary>
/// <param name="Anchor">The best-scoring member. Its mass and retention time anchor the species.</param>
/// <param name="Members">Every identification grouped into the species, the anchor included.</param>
/// <param name="MinCharge">Lowest member precursor charge minus 2, floored at 2.</param>
/// <param name="MaxCharge">Highest member precursor charge plus 2, capped at 80.</param>
public sealed record IdentifiedSpecies(
    MmResultRecord Anchor,
    MmResultRecord[] Members,
    int MinCharge,
    int MaxCharge);

/// <summary>
/// Groups identifications that describe the same MS1 signal, so each is modelled once.
/// </summary>
/// <remarks>
/// <para>
/// Two identifications are one species when they elute within <see cref="DefaultRtTolerance"/>
/// minutes and their masses differ by a whole number of isotopologue spacings (0 to ±3) within
/// <see cref="DefaultMassTolerance"/> Da. That covers:
/// </para>
/// <list type="bullet">
/// <item>the same proteoform identified at several precursor charges;</item>
/// <item>identical sequences under different accessions;</item>
/// <item>isobaric localization variants;</item>
/// <item>off-by-one-dalton monoisotopic assignments;</item>
/// <item>near-isobaric pairs such as deamidation (+0.984 Da), which top-down MS1 cannot separate
/// from a +1 isotopologue shift.</item>
/// </list>
/// <para>
/// Modelled separately, each of those claims the whole observed envelope; at H2B in Jurkat rep2
/// fract7 seven such models each predicted the full height of the same peak. Grouping is greedy in
/// descending score, so the best-scoring identification anchors each species.
/// </para>
/// </remarks>
public static class SpeciesGrouper
{
    public const double DefaultRtTolerance = 0.5;
    public const double DefaultMassTolerance = 0.03;

    /// <summary>Isotopologue spacing of averagine, in daltons.</summary>
    public const double AveragineIsotopeSpacing = 1.00235;

    public static IdentifiedSpecies[] Group(
        IEnumerable<MmResultRecord> records,
        double rtTolerance = DefaultRtTolerance,
        double massTolerance = DefaultMassTolerance)
    {
        if (records is null) throw new ArgumentNullException(nameof(records));

        var groups = new List<List<MmResultRecord>>();
        foreach (var record in records.OrderByDescending(r => r.Score).ThenBy(r => r.RetentionTime))
        {
            var group = groups.FirstOrDefault(g =>
                Math.Abs(g[0].RetentionTime - record.RetentionTime) <= rtTolerance
                && WithinIsotopeSpacings(g[0].MonoisotopicMass, record.MonoisotopicMass, massTolerance));

            if (group is null)
                groups.Add(new List<MmResultRecord> { record });
            else
                group.Add(record);
        }

        return groups.Select(g => new IdentifiedSpecies(
            g[0],
            g.ToArray(),
            Math.Max(2, g.Min(r => r.PrecursorCharge) - 2),
            Math.Min(80, g.Max(r => r.PrecursorCharge) + 2))).ToArray();
    }

    /// <summary>Each identification as its own species, with no grouping.</summary>
    public static IdentifiedSpecies[] Ungrouped(IEnumerable<MmResultRecord> records) =>
        records.Select(r => new IdentifiedSpecies(
            r, new[] { r }, Math.Max(2, r.PrecursorCharge - 2), Math.Min(80, r.PrecursorCharge + 2))).ToArray();

    public static bool WithinIsotopeSpacings(double a, double b, double tolerance)
    {
        double delta = b - a;
        int n = (int)Math.Round(delta / AveragineIsotopeSpacing);
        return Math.Abs(n) <= 3 && Math.Abs(delta - n * AveragineIsotopeSpacing) <= tolerance;
    }
}
