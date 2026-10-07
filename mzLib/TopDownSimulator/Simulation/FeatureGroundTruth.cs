#nullable enable
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Text;
using Chemistry;
using TopDownSimulator.Model;

namespace TopDownSimulator.Simulation;

/// <summary>
/// One simulated feature: a single (proteoform, charge) that survived scan reduction, described in
/// the terms a feature finder emits rather than in the terms that generated it.
/// </summary>
/// <remarks>
/// This is the artifact a benchmark scores against. The parameter sidecar written by
/// <see cref="Simulator.WriteGroundTruth"/> reproduces a run; it cannot be scored without
/// reimplementing the simulator, because deciding which charges are present, where each
/// isotopologue lands, and which scans clear the reduction threshold is exactly what the simulator
/// does.
/// </remarks>
public sealed record SimulatedFeature(
    string Identifier,
    double MonoisotopicMass,
    int Charge,
    double MonoisotopicMz,
    double ApexMz,
    double ApexRt,
    double RtStart,
    double RtEnd,
    int ApexScanNumber,
    int FirstScanNumber,
    int LastScanNumber,
    double SummedIntensity,
    double ApexIntensity,
    int NumIsotopologues);

/// <summary>
/// Derives the feature-level ground truth for a simulation and writes it as TSV.
/// </summary>
public static class FeatureGroundTruth
{
    /// <summary>
    /// Enumerates every (proteoform, charge) whose own simulated signal clears
    /// <paramref name="intensityFloor"/> in at least one scan.
    /// </summary>
    /// <remarks>
    /// Presence is decided on the feature's <b>own</b> contribution, not on the total intensity in
    /// the file. The total is the sum over all proteoforms and so is never smaller, which makes
    /// this the conservative rule: everything listed here is genuinely above the floor in the mzML.
    /// A charge state that reduction removed entirely is absent, so a benchmark does not start with
    /// unfixable false negatives.
    /// <para>
    /// The apex m/z is snapped to the nearest signal centroid in the apex scan, so the truth quotes
    /// a position that is actually in the file. Where an envelope's isotopologues are not resolved
    /// from each other or from a neighbour's, that is the merged maximum the instrument would
    /// report, not the theoretical isotopologue position.
    /// </para>
    /// </remarks>
    /// <param name="oneBasedScanNumbers">Scan numbers parallel to <paramref name="scanTimes"/>.</param>
    /// <param name="scanMzs">Per scan, the ascending m/z of the signal centroids written to it.</param>
    public static SimulatedFeature[] Build(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        double[] scanTimes,
        int[] oneBasedScanNumbers,
        IReadOnlyList<double[]> scanMzs,
        IIntensityFloor intensityFloor)
    {
        if (intensityFloor is null) throw new ArgumentNullException(nameof(intensityFloor));
        if (scanTimes is null) throw new ArgumentNullException(nameof(scanTimes));
        var floors = new IIntensityFloor[scanTimes.Length];
        Array.Fill(floors, intensityFloor);
        return Build(proteoforms, minCharge, maxCharge, widthModel, scanTimes, oneBasedScanNumbers, scanMzs, floors);
    }

    /// <summary>
    /// As above, with the detection limit of scan s given by <paramref name="scanFloors"/>[s], for
    /// simulations whose noise floor changes from scan to scan.
    /// </summary>
    public static SimulatedFeature[] Build(
        IReadOnlyList<ProteoformModel> proteoforms,
        int minCharge,
        int maxCharge,
        IPeakWidthModel widthModel,
        double[] scanTimes,
        int[] oneBasedScanNumbers,
        IReadOnlyList<double[]> scanMzs,
        IReadOnlyList<IIntensityFloor> scanFloors)
    {
        if (scanFloors is null) throw new ArgumentNullException(nameof(scanFloors));
        if (proteoforms is null) throw new ArgumentNullException(nameof(proteoforms));
        if (widthModel is null) throw new ArgumentNullException(nameof(widthModel));
        if (scanTimes is null) throw new ArgumentNullException(nameof(scanTimes));
        if (oneBasedScanNumbers is null) throw new ArgumentNullException(nameof(oneBasedScanNumbers));
        if (scanMzs is null) throw new ArgumentNullException(nameof(scanMzs));
        if (oneBasedScanNumbers.Length != scanTimes.Length)
            throw new ArgumentException("Scan numbers must be parallel to scan times.", nameof(oneBasedScanNumbers));
        if (scanMzs.Count != scanTimes.Length)
            throw new ArgumentException("There must be one peak list per scan.", nameof(scanMzs));
        if (scanFloors.Count != scanTimes.Length)
            throw new ArgumentException("There must be one floor per scan.", nameof(scanFloors));
        if (minCharge < 1 || maxCharge < minCharge)
            throw new ArgumentException("Charge range must satisfy 1 ≤ minCharge ≤ maxCharge.");

        var features = new List<SimulatedFeature>();
        int nScans = scanTimes.Length;
        var rtValues = new double[nScans];

        for (int p = 0; p < proteoforms.Count; p++)
        {
            var model = proteoforms[p];
            var kernel = new IsotopeEnvelopeKernel(model.MonoisotopicMass);

            double maxRt = 0;
            for (int s = 0; s < nScans; s++)
            {
                double rt = model.RtProfile.Evaluate(scanTimes[s]);
                rtValues[s] = rt > 0 ? rt : 0;
                if (rtValues[s] > maxRt)
                    maxRt = rtValues[s];
            }

            if (model.Abundance <= 0 || maxRt <= 0)
                continue;

            string identifier = model.Identifier ?? $"proteoform{p}";

            for (int z = minCharge; z <= maxCharge; z++)
            {
                double fz = model.ChargeDistribution.Evaluate(z);
                if (fz <= 0)
                    continue;

                double[] centroids = kernel.CentroidMzs(z);
                int nIso = centroids.Length;
                if (nIso == 0)
                    continue;

                // Per-isotopologue unit-abundance peak height, including the overlap contributed by
                // neighbouring isotopologues of the same charge — which is what makes this agree
                // with the rendered spectrum once the envelope starts to merge at high m/z.
                var peakShape = new double[nIso];
                double maxShape = 0;
                for (int i = 0; i < nIso; i++)
                {
                    peakShape[i] = kernel.Evaluate(centroids[i], z, widthModel);
                    if (peakShape[i] > maxShape)
                        maxShape = peakShape[i];
                }

                // The brightest this (proteoform, charge) can ever be, against the lowest floor
                // anywhere under its envelope in any scan. Charge states far out on the Gaussian
                // never clear the floor, and skipping them here is what keeps this cheap over a
                // wide global charge range. Taking the minimum floor keeps the early-out
                // conservative: a charge state dropped here could not have cleared the floor at
                // any m/z it occupies in any scan.
                double lowestFloor = double.PositiveInfinity;
                for (int s = 0; s < nScans; s++)
                    lowestFloor = Math.Min(lowestFloor, scanFloors[s].MinOver(centroids[0], centroids[^1]));
                if (model.Abundance * maxRt * fz * maxShape < lowestFloor)
                    continue;

                int firstScan = -1, lastScan = -1, apexScan = -1;
                double apexChargeIntensity = 0, summedIntensity = 0;
                double apexPeakIntensity = 0, apexPeakMz = 0;
                int apexIsotopologueCount = 0;

                for (int s = 0; s < nScans; s++)
                {
                    var intensityFloor = scanFloors[s];
                    double scale = model.Abundance * rtValues[s] * fz;
                    if (scale * maxShape < intensityFloor.MinOver(centroids[0], centroids[^1]))
                        continue;

                    double scanSum = 0, scanBest = 0, scanBestMz = 0;
                    int scanIsotopologues = 0;

                    for (int i = 0; i < nIso; i++)
                    {
                        double value = scale * peakShape[i];
                        if (value < intensityFloor.At(centroids[i]))
                            continue;

                        scanSum += value;
                        scanIsotopologues++;
                        if (value > scanBest)
                        {
                            scanBest = value;
                            scanBestMz = centroids[i];
                        }
                    }

                    if (scanIsotopologues == 0)
                        continue;

                    if (firstScan < 0)
                        firstScan = s;
                    lastScan = s;
                    summedIntensity += scanSum;

                    if (scanSum > apexChargeIntensity)
                    {
                        apexChargeIntensity = scanSum;
                        apexScan = s;
                        apexPeakIntensity = scanBest;
                        apexPeakMz = scanBestMz;
                        apexIsotopologueCount = scanIsotopologues;
                    }
                }

                if (apexScan < 0)
                    continue;

                features.Add(new SimulatedFeature(
                    Identifier: identifier,
                    MonoisotopicMass: model.MonoisotopicMass,
                    Charge: z,
                    // Computed from the mass rather than read off centroids[0]. The averagine tables
                    // drop isotopologues below an absolute probability of 1e-8, and the all-light
                    // composition falls under that at roughly 31 kDa — above which the kernel's
                    // lowest-mass entry is not the monoisotopic peak at all. This column is an
                    // identifier a finder's reported monoisotopic m/z is scored against, so it must
                    // be exact whether or not that peak survived into the envelope.
                    MonoisotopicMz: model.MonoisotopicMass.ToMz(z),
                    ApexMz: SnapToAxis(scanMzs[apexScan], apexPeakMz),
                    ApexRt: scanTimes[apexScan],
                    RtStart: scanTimes[firstScan],
                    RtEnd: scanTimes[lastScan],
                    ApexScanNumber: oneBasedScanNumbers[apexScan],
                    FirstScanNumber: oneBasedScanNumbers[firstScan],
                    LastScanNumber: oneBasedScanNumbers[lastScan],
                    SummedIntensity: summedIntensity,
                    ApexIntensity: apexPeakIntensity,
                    NumIsotopologues: apexIsotopologueCount));
            }
        }

        return features.ToArray();
    }

    public static void Write(IReadOnlyList<SimulatedFeature> features, string outputPath)
    {
        if (features is null) throw new ArgumentNullException(nameof(features));
        if (string.IsNullOrWhiteSpace(outputPath))
            throw new ArgumentException("An output path is required.", nameof(outputPath));

        var sb = new StringBuilder();
        sb.AppendLine(string.Join('\t',
            "Identifier", "MonoisotopicMass", "Charge", "MonoisotopicMz", "ApexMz",
            "ApexRt", "RtStart", "RtEnd",
            "ApexScanNumber", "FirstScanNumber", "LastScanNumber",
            "SummedIntensity", "ApexIntensity", "NumIsotopologues"));

        foreach (var f in features)
        {
            sb.AppendLine(string.Join('\t',
                f.Identifier,
                Num(f.MonoisotopicMass),
                Int(f.Charge),
                Num(f.MonoisotopicMz),
                Num(f.ApexMz),
                Num(f.ApexRt),
                Num(f.RtStart),
                Num(f.RtEnd),
                Int(f.ApexScanNumber),
                Int(f.FirstScanNumber),
                Int(f.LastScanNumber),
                Num(f.SummedIntensity),
                Num(f.ApexIntensity),
                Int(f.NumIsotopologues)));
        }

        File.WriteAllText(outputPath, sb.ToString());
    }

    private static double SnapToAxis(double[] axis, double mz)
    {
        if (axis.Length == 0)
            return mz;

        int lo = 0, hi = axis.Length;
        while (lo < hi)
        {
            int mid = lo + ((hi - lo) >> 1);
            if (axis[mid] < mz)
                lo = mid + 1;
            else
                hi = mid;
        }

        if (lo == 0)
            return axis[0];
        if (lo == axis.Length)
            return axis[^1];

        return mz - axis[lo - 1] <= axis[lo] - mz ? axis[lo - 1] : axis[lo];
    }

    private static string Num(double value) => value.ToString("R", CultureInfo.InvariantCulture);
    private static string Int(int value) => value.ToString(CultureInfo.InvariantCulture);
}
