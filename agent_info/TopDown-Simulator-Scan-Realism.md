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

## What was done (2026-10-07)

| Commit | Change |
|---|---|
| `eefc593` | `ProfileCentroider`: one centroid per local maximum of the summed profile, Newton-polished so the height is the forward model at the maximum. Feature truth snaps apex m/z to the apex scan's own peaks. |
| `0719ddc3` | The refit basis is built through an exact m/z index of (model, charge) envelopes, so the 200-model cap is gone (default 10000). Records are grouped into species (0.5 min, 0 to ±3 isotope spacings within 0.03 Da) and fitted once, over every member's charge ± 2. |
| `4e265786` | `ScanNoiseConditions.FromSourceScans`: one noise model per scan, with amplitude `JurkatNoiseTimesInjectionTime` / IT and per-bin density counted from the source scan's peaks under S/N 10. |

**Density is conditioned on the source, not predicted.** `NoiseDensityAcrossRun` shows noise-like peak
counts follow neither IT nor TIC: about 350 per scan pre-elution and 17k in the wash, both at the 50 ms
maximum, and the m/z shape changes too. Most of it is presumably faint unidentified analyte and background.

The full rep2 export with all of this ran as `MZLIB_TOPDOWN_SIM_OUTPUT_TAG=.v2`:

- 1089 records grouped into 474 species, and the refit over 473 models took 25 s;
- the unexplained energy fell from 0.72 to 0.45;
- the whole export took 72 s and is 427 MB;
- it averages 7058 peaks per scan, against 7210 real.

Harness, `file` rows, before (`full`) → after (`full.v2`):

| scan / region | peaks (real) | sim matched | √cos | TIC ratio | KS |
|---|---|---|---|---|---|
| 2534, 860.4–864.4 | 251 → 88 (72) | 0.28 → 0.72 | 0.76 → 0.94 | 152 → 1.26 | 0.42 → 0.20 |
| 2534, whole scan | 20576 → 8137 (8222) | 0.22 → 0.33 | 0.61 → 0.61 | 60 → 1.14 | 0.89 → 0.15 |
| 1840, whole scan | 17228 → 11417 (11565) | 0.27 → 0.35 | 0.25 → 0.29 | 0.42 → 1.08 | 0.99 → 0.18 |
| 1200, pre-elution | 16969 → 359 (344) | – | – | 503 → 5.0 | – |
| 3680, wash | 16941 → 16985 (17134) | – | – | 1.75 → 1.07 | – |

## Still open

- **Busy scans are mostly unidentified analyte.** In scan 1840 only 25 simulated signal peaks clear the
  floor; the real scan's structure (cos 0.19) is analyte the simulation renders as unstructured noise.
  Only identified proteoforms are simulated.
- **The level's m/z shape is the wash shape.** Busy scans rise only about 2× from m/z 650 to 1150, against
  about 3.8× in the model, so the floor above m/z 1000 sits too high there.
- **Pre-elution TIC is about 5× real.** Empty scans have amplitude × IT ≈ 9e3, not 2e4. The S/N 10
  counting window makes the sampled noise a little bright there.
- **The refit still reports `converged: false`** (see Known-Bugs).
