# Top-Down Simulator — Noise Model

**Added:** 2026-07-26
**Code:** `mzLib/TopDownSimulator/Noise/`
**Measurement tool:** `mzLib/Test/TopDownSimulator/NoiseCharacterization.cs` (`[Explicit]`)
**Tests:** `NoiseFloorModelTests`, `NoiseInjectorTests`, `NoisyMzmlExportTests`

This closes the Phase 4 gap recorded in `TopDown-Simulator-Status.md` ("the biggest gap: there is no
noise model").

## Do we need to simulate transients and FFT them?

**No.** It was considered and rejected on three grounds:

1. **The output is centroided.** mzLib's mzML reader rejects profile-mode files, so the simulator
   writes centroids only. Everything an FFT gives you beyond a peak's (position, intensity, width) is
   discarded at the centroiding step. What survives is exactly the four things measured below.
2. **Cost.** Covering m/z 600–2000 at R≈100 000 needs a ~10⁶-point transform per scan, ~2900 scans
   per run. That is minutes to hours, against microseconds for direct centroid draws.
3. **It would have been wrong anyway.** The one thing a transient model buys for free is the
   analytic density law — noise uniform in frequency, and since f ∝ (m/z)^−1/2 the density in m/z
   goes as (m/z)^−1.5. Measured against real data that prediction is off by two orders of magnitude
   (see below). The reported peaks are local maxima above a detection threshold, not occupied
   frequency bins, and no amount of transient fidelity fixes that.

A transient model would be the right tool for one thing only: profile-mode peak coalescence at high
ion load. That is moot while the output is centroided.

## What the real noise actually is

Measured on the 50–55 min post-elution window of four Jurkat runs — rep1 fract6, rep1 fract7,
rep2 fract5, rep2 fract7 — read straight off the Thermo `CentroidStream`, whose `Noises`,
`Baselines` and `Resolutions` columns mzLib fetches and then discards
(`ThermoRawFileReader.cs:417`, `noiseData: null`).

### The quiet window has *more* peaks than the busy one

| Window | peaks/scan (median) |
|---|---|
| 31–35 min (busy) | 8 419 – 9 867 |
| 50–55 min (quiet) | 16 434 – 17 554 |

Roughly twice as many. Almost certainly the AGC: with little analyte the injection time runs long,
the transient fills, and thousands of noise maxima clear the reporting threshold.

For scale, a clean simulation of rep2 fract7's 1087 q≤0.01 proteoforms over 45–56 min emits
~10 600 peaks/scan against the real file's ~15 500 in that window — and every one of those simulated
peaks is a real isotopologue of a real proteoform, whereas most of the real file's are noise. A
smaller model set makes the gap far worse: the 25-proteoform run recorded in
`TopDown-Simulator-Status.md` produced ~32 peaks/scan.

### It is unstructured — not chemical background

Three independent measurements agree, and they are why the model has no persistent-contaminant
component:

- The instrument assigned **no charge to 98.4 %** of quiet-window peaks (88 % in the busy window).
- The fraction of peaks with a +1.00335/z partner was ~45 % **identically for every z from 1 to 7**.
  A flat response across implausible charges is chance collision at this peak density, not isotopic
  structure.
- The **brightest 10 %** of quiet-window peaks persisted to the next scan no more often than random
  ones (42–48 % either way, against a ~25 % chance-collision floor). In the busy window that number
  jumps to **77 %**, which is what real analyte looks like.

### Density vs m/z — the analytic law fails

Peaks per scan per 100-Th bin, mean of the four files:

| m/z | 600 | 700 | 800 | 900 | 1000 | 1100 | 1200 | 1300 | 1400 | 1500 | 1600 | 1700 | 1800 | 1900 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| n/scan | 3527 | 4303 | 3504 | 2519 | 1337 | 745 | 383 | 235 | 141 | 87.5 | 50.8 | 29.0 | 16.4 | 8.2 |

Rises to a maximum at 700–800 and then falls by a factor of ~500. A (m/z)^−1.5 law predicts a
monotonic fall of about 4× over the same span. Hence: **measure the curve, do not derive it.**

### Noise amplitude vs m/z — same shape everywhere, amplitude varies 3×

Relative to the 600–700 bin, the four files agreed to within a few percent:

| m/z | 600 | 700 | 800 | 900 | 1000 | 1100 | 1200 | 1300 | 1400 | 1500 | 1600 | 1700 | 1800 | 1900 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| relative | 1.00 | 1.47 | 2.04 | 2.68 | 3.26 | 3.77 | 4.04 | 4.47 | 4.63 | 4.57 | 4.55 | 4.07 | 3.82 | 3.32 |

Absolute amplitude at m/z 650 was 663 (rep2 f7), 715 (rep1 f7), 1640 (rep2 f5), 2050 (rep1 f6) —
a factor of 3. So amplitude is the model's one free parameter; shape is fixed.

### S/N is universal

Quantiles of intensity ÷ the instrument's own per-peak noise estimate, quiet window:

| | p1 | p25 | p50 | p75 | p95 | p99 |
|---|---|---|---|---|---|---|
| rep1 f7 | 1.37 | 2.08 | 2.75 | 3.90 | 6.68 | 9.96 |
| rep2 f5 | 1.37 | 2.10 | 2.79 | 4.02 | 7.28 | 11.65 |
| rep1 f6 | 1.37 | 2.10 | 2.79 | 4.00 | 7.00 | 10.77 |
| rep2 f7 | 1.38 | 2.09 | 2.76 | 3.91 | 6.69 | 9.95 |

Identical across files **and** across both RT windows, diverging only at p99 in the busy window
where real analyte lives. A shifted exponential `S/N = 1.565 + Exp(mean 1.707)` fits the body to
under 2 %.

### Within-m/z amplitude dispersion

Pooling peaks in 5-Th slices and dividing each by its own slice median gives a right-skewed
multiplicative spread: p25 ≈ 0.72, p50 = 1, p75 ≈ 1.58, p95 ≈ 2.8, i.e. σ_log ≈ 0.6. Feeding 0.6 in
verbatim over-widens the simulated intensities, because within a slice the local noise estimate and
the S/N a peak achieves are negatively correlated and the model treats them as independent. The
shipped `LevelDispersion = 0.42` is fitted instead to reproduce the measured **intensity** quantile
ratios, which is what a feature finder sees.

## The model

```
intensity = NoiseLevelAt(m/z) · LogNormal(0, 0.42) · (1.565 + Exp(mean 1.707))
```

with m/z drawn from the measured density curve, counts Poisson per bin, and every scan independent.

Verified at the sampler: simulated intensity quantile ratios against rep1 fract7's measured ones come
out at p25/p50 = 0.974×, p75/p50 = 1.025×, p95/p50 = 1.002× of measured
(`SampledIntensityDistributionHasTheMeasuredShape`).

Verified end to end against the real file, 50–55 min window of rep2 fract7, comparing the written
mzML to the raw (`NoiseCharacterization.CompareSimulatedMzmlToRealFile`): **16 923 simulated
peaks/scan against 16 370 real**, with per-100-Th density ratios between 0.95 and 1.09 across the
whole range, and median intensity ratios between 0.80 and 1.25.

The intensity ratio drifts monotonically — 0.80 at m/z 650, 1.25 at m/z 1500 — because the shipped
relative-level curve is the four-file average while the amplitude supplied was rep2 fract7's own.
Using a per-file shape would remove that drift; whether it is worth the extra parameter is untested.

### Threshold coherence

Injecting noise also switches reduction from `SimulatedScanReducer`'s
fraction-of-the-brightest-peak floor onto `NoiseRelativeFloor`, an m/z-dependent multiple of the
local noise — the rule an instrument actually applies. Without this the file would contain noise
peaks fainter than any surviving signal peak, which no instrument can produce. This is why
`IIntensityFloor` replaced the scalar floor in `SimulatedScanReducer.ComputeFloor` and
`FeatureGroundTruth.Build`.

Note that not every noise peak clears the nominal floor, and that is correct: the floor uses the
median amplitude curve while individual peaks scatter about it by `LevelDispersion`.

### Two merge rules that are easy to get wrong

`NoiseInjector.Collapse` collapses peaks the instrument could not have reported separately. Two
things about it are load-bearing at this peak density:

- **Groups are bounded by their own first member, not grown transitively.** Single-linkage chaining
  lets a run of peaks each just inside the tolerance collapse into one centroid spanning many σ. At
  the measured density the mean noise spacing near m/z 700 is about 2σ, so chaining merged 5.1 M of
  9.0 M peaks in the first real run — over half the spectrum.
- **Groups made purely of noise are left alone.** The density curve was measured from the centroids
  a real instrument *reported*, so it already accounts for whatever that instrument's peak detection
  merged. Merging noise against noise applies the same reduction twice. Groups containing a signal
  peak are collapsed, because there the coincidence is genuine.

### Cost

At `DensityScale = 1.0` the output carries ~16 900 extra peaks per scan. `DensityScale` trades
fidelity for size linearly, and `MZLIB_TOPDOWN_SIM_NOISE_DENSITY_SCALE` sets it for the export tests.

Two `[Explicit]` exporters in `Test/FileReadingTests/AnalysisExample.cs`, both on rep2 fract7:

| Test | Window | Output |
|---|---|---|
| `ExportRep2NoisySimulationForInspection` | 45–56 min, 334 scans | clean 73 MB + noisy 119 MB |
| `ExportRep2FullNoisySimulation` | whole run, ~2900 scans | noisy only |

The windowed one writes a matched clean/noisy pair for eyeballing against the real data; the full-run
one is what a benchmark should actually consume, since a slice cannot exercise anything depending on
elution order, the full dynamic range, or reaching the post-elution region from a busy one.

## Scan-to-scan jitter on real peaks

**Code:** `Noise/PeakJitterModel.cs` · **Measurement:** `Test/TopDownSimulator/SignalJitterCharacterization.cs`

Without this, simulated XICs are the analytical EMG evaluated at each scan time — perfectly smooth
and perfectly reproducible, which is nothing like real data and is what XIC smoothing and apex
detection actually contend with.

Measured on rep2 fract7 against the top 400 MetaMorpheus identifications (26 730 isotopologue
observations), using estimators built to **cancel the elution profile rather than model it**:

- m/z error is each observation's deviation from its own isotopologue's median, split into a
  per-scan common mode and a per-peak remainder.
- Intensity error is the scan-to-scan scatter of the **ratio between adjacent isotopologues of one
  species**. That ratio is fixed by chemistry, so the elution profile divides out exactly and what
  remains is measurement noise.

### m/z

Per-scan common-mode offset (calibration drift): **σ = 1.89 ppm**. Per-peak remainder: **σ = 3.98 ppm**
overall, so the per-peak term dominates but neither is negligible. They stress a feature finder
differently — a common-mode shift moves a whole spectrum together and can be calibrated out;
independent per-peak error cannot.

| S/N | 2–4 | 4–8 | 8–16 | 16–32 | 32–64 | 64+ |
|---|---|---|---|---|---|---|
| measured σ (ppm) | 5.447 | 3.842 | 2.766 | 1.780 | 1.579 | 1.137 |

`σ × (S/N)` is not constant but `σ × √(S/N)` nearly is, so the law is
**σ_ppm = √(9.06²/(S/N) + 0.81²)**, which reproduces all six bins to within 15 %.

### Intensity

| S/N | 2–4 | 4–8 | 8–16 | 16–32 | 32–64 | 64–256 |
|---|---|---|---|---|---|---|
| measured σ of log intensity | 0.6318 | 0.4622 | 0.3644 | 0.2883 | 0.2368 | 0.1881 |

**σ_log = √(1.026²/(S/N) + 0.165²)**, reproducing all six bins to within 6 %.

**The 0.165 floor is the headline number**: even at S/N in the hundreds a real peak still scatters
~17 % scan to scan, where the simulator produced exactly 0 %.

The multiplicative factor is drawn **mean-preserving**, `exp(N(−σ²/2, σ))`. At the σ ≈ 0.84 that a
peak near the detection floor carries, a median-preserving draw would inflate expected intensity by
43 %, quietly brightening exactly the faint peaks whose detectability the simulation exists to test.

### XIC dropout falls out of it

A peak just above the floor has S/N ≈ 1.6 and hence σ_log ≈ 0.84, so it drops below the detection
limit in a good fraction of scans — measured at **61 %** for a peak sitting 5 % over the floor, and
**0 %** for a bright one. That gappy-XIC behaviour is real and is a genuine stressor, and it emerges
from the model rather than being added to it. Consequence: a feature's first and last scan numbers in
`FeatureGroundTruth` become approximate under jitter, though its identity, charge, mass and apex do
not.

### Two things that are easy to get wrong

- **Noise peaks are not jittered.** They are drawn from a distribution measured on centroids a real
  instrument already reported, so their scatter is baked in; jittering again double-counts it.
- **Jitter and noise draw from independent RNG streams** (`NoiseInjector.JitterStream` /
  `NoiseStream`). Sharing one stream means enabling jitter consumes draws before the noise sampler
  runs, so the same seed yields a different noise floor and a clean A/B comparison becomes
  impossible.

### Unmeasured

`CommonModeLogIntensitySigma` defaults to **0 because it is unmeasured, not because it is believed
absent**. The isotopologue-ratio estimator divides out anything affecting both peaks equally, so it
is blind to precisely this term, and a real AGC fluctuation would produce one.

Nor was the m/z jitter binned by m/z. The width law says σ_m ∝ (m/z)^1.5, so jitter expressed in ppm
should grow as √(m/z); the pooled figure above is dominated by wherever these proteoforms sit
(mostly m/z 600–1200).

## Known gaps / follow-ups

1. **Noise density does not vary with retention time — now the largest remaining error.** The model
   emits one density everywhere, calibrated on the 50–55 min window. The real profile, measured by
   `NoiseCharacterization.CharacterizeDensityAcrossRun`, spans a factor of **40**:

   | RT (min) | 0–10 | 10–20 | 20–30 | 30–40 | 40–50 | 50–65 | 65–75 | 75–90 |
   |---|---|---|---|---|---|---|---|---|
   | rep2 f7 peaks/scan | ~400 | 2000–3400 | 400–1400 | 7400–11400 | 11200–14300 | 14700–16500 | 9300–12200 | 2500–3900 |

   It rises with the gradient, peaks at 50–65 min, and collapses at the wash — so it tracks
   chromatographic background, not simply the absence of analyte. Consequence: over a full run the
   real file averages **7 210 peaks/scan while the simulation puts 17 154 everywhere**, roughly 2.4×
   too noisy overall and ~40× too noisy in the first ten minutes. Within the calibrated window the
   agreement is excellent (16 923 vs 16 370). Fixing this means an RT-dependent density scale plus a
   fitter that reads the profile off a target file; the measurement already exists.
2. **No harmonics or satellite artifacts.** These stress charge assignment specifically. Not
   measured yet; the charge-call and isotope-partner statistics above would not have detected them.
3. **The noise amplitude is set by hand.** There is no fitter that reads it off a target raw file —
   run `NoiseCharacterization.CharacterizeNoise` and read the median noise in the 600–700 bin.
4. **`ThermoRawFileReader` still discards `CentroidStream.Noises`/`.Baselines`.** Wiring them into
   `MsDataScan.NoiseData` (which already has the 3×N layout the mzML writer emits) would let the
   amplitude be fitted automatically, and would let simulated files carry a noise array too.
5. **Common-mode intensity jitter and the m/z dependence of m/z jitter are unmeasured** — see the
   jitter section above.
