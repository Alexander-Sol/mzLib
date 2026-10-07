# Top-Down Simulator — Scan Realism

**Written:** 2026-10-07, on `TopDownSimClean`.
**Question:** why simulated scans show many peaks between the isotopologues of a strong, resolved envelope
where the real scan shows almost none. The reference case is rep2 fract7 real scan 2534, which is simulated
scan 1291 (RT 41.13), over m/z 860.4–864.4. That is H2B at z=16.

## The harness

`mzLib/Test/TopDownSimulator/ScanRealismComparison.cs`, all `[Explicit]`:

- `CompareSimulatedScansToReal` rebuilds the models from the ground-truth sidecar of the full noisy export
  (`*.full.noisy.simulated.groundtruth.tsv`). It re-simulates four target scans, one per regime (busy at
  IT 0.3 ms, busy at IT 2.8 ms, post-elution, pre-elution), and scores each variant against the real scan.
  It takes about 30 s and writes `D:\JurkatTopdown\scan-realism\comparison.tsv`, plus each spectrum as
  `scan<N>.<variant>.tsv`.
  - To try an idea, add a `Variant`.
  - `file` (the mzML on disk) and `current` (the re-simulation) should agree. If they don't, the sidecar no
    longer describes the file.
- `AttributeSimulatedPeak` lists which models put intensity on one simulated peak. It also shows what the
  real centroid stream and mzLib's reader report there; the two agree.
- `NoiseLevelVsInjectionTime` reports per-scan injection time (IT) against the noise level Thermo reports
  for each centroid.

The metrics are in `TopDownSimulator/Comparison/CentroidSpectrumComparison.cs`:

- peak counts;
- the matched fraction each way (10 ppm, one-to-one);
- cosine, and cosine on √intensity;
- the ratio of summed intensity (TIC);
- the KS distance between the two distributions of log(I/Imax);
- the fraction of peaks under 1 % of the maximum.

**Low `sim_matched` is the signature of spurious peaks.**

## Findings

### 1. The peaks between isotopologues are simulated signal, not noise (main cause)

`GridRasterizer.RasterizeAtCentroids` evaluates the *summed* profile at the union of every model's
isotopologue m/z values. In a crowded region, one real peak is therefore sampled at many nearby positions.

- The samples on its shoulders survive as peaks between the isotopologues.
- `NoiseInjector.Collapse` then *sums* the samples that fall within 1σ, which counts the same profile height
  several times.

A real instrument reports one centroid per local maximum. Prototype variant `picked-signal` renders the
profile on a dense grid and keeps the local maxima, with a log-parabola refinement. For 860.4–864.4:

| variant | peaks (real 72) | sim matched | √cos | TIC ratio |
|---|---|---|---|---|
| current | 256 | 0.28 | 0.77 | 153 |
| signal-only (no noise peaks) | 238 | 0.30 | 0.80 | 153 |
| picked-signal | 64 | 0.98 | 0.97 | 5.8 |

The noise floor contributes only about 18 of the 256 peaks.

### 2. The fitted abundances over-count crowded families

Even after picking, the summed intensity in that region is about 6× the real value.

`AttributeSimulatedPeak` at 861.979 finds 43 models contributing. Seven of them each predict about the
whole real peak: H2B masses 13765–13768 Da, at off-by-one-Da monoisotopic errors and +0.98 Da deamidations.

The causes:

- The global refit was **skipped** for the full run, because 1087 models exceeds
  `GLOBAL_REFIT_MAX_MODELS` (default 200).
- `DeduplicateByProteoform` keys on accession, full sequence and *precursor charge*. That keeps one model
  per charge, and identical sequences under different accessions survive as separate models.

Deduplicating by mass (0.05 Da, 0.5 min) roughly halves the excess. The real fix is in fitting.

### 3. The noise amplitude scales as 1/IT, so a single level cannot fit every scan

Wherever AGC limits the fill, noise × IT ≈ 2×10⁴ at m/z 650. In empty scans at the 50 ms maximum IT it is
≈ 9×10³.

Across the run, the noise level at 650 ranges from 180 to 49,000. The export used a constant 663.

In scan 2534 (IT 0.31 ms) the real noise is about 73k. The `noise-it` variant (2e4/IT) gives about 65k.

### 4. Noise-peak density is not stationary

Measured peaks per scan:

| regime | peaks per scan |
|---|---|
| pre-elution | ~400 |
| busy | 8–12k, which includes signal |
| post-elution | ~17k |

The model emits 17k in every scan. After peak picking, this is what still clutters the reference region:
177 simulated peaks against 72 real. Pre-elution scans are about 50× too dense. The calibration window
(50–55 min) is the densest part of the run, so it is not representative.

## Open next steps

1. Move peak picking into the library centroid path, replacing union sampling plus summing collapse. Check
   what `FeatureGroundTruth` assumes about the centroid axis first.
2. Make the noise level per scan from IT. The real IT is available because the simulation reuses the real
   scan grid; an AGC model driven by the simulated TIC is the forward-model alternative.
3. Make the noise density per scan. It needs a measurement first, e.g. a per-scan count of the noise-like
   peaks (charge 0, low S/N), and against what it scales.
4. Fitting: refit in RT blocks so 1087 models fit under the cap, and deduplicate by mass rather than by
   (accession, sequence, charge).
