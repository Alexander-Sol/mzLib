# Top-Down Simulator — Agent Context

**Written:** 2026-10-06
**Branch:** `TopDownSimClean`, the simulator alone squashed onto `master` on 2026-10-06. The full
history is on `TopDownSimPolished`.
**Sources:** Claude memory files for this repo, plus the one Claude session that ran on this branch.
Read this first, then the topic docs listed under "Where the detail lives".

## What this work is

A parametric forward model for top-down MS1 (and basic MS2) data, in `mzLib/TopDownSimulator/`.
It fits per-proteoform models (EMG elution profile, charge-state distribution, isotope envelope,
abundance) to a real run using MetaMorpheus IDs, then rasterizes synthetic centroided mzML with
feature-level ground truth. The intended use is **benchmarking feature finders and deconvolution**
(recall/precision against known truth), not precise quantification.

The source run for everything so far is Jurkat top-down rep2 fract7:

- Raw: `D:\JurkatTopdown\02-18-20_jurkat_td_rep2_fract7.raw`
- IDs: `D:\JurkatTopdown\Frac7_GPTMD_Search\Task2-TopDownSearch\Individual File Results\02-18-20_jurkat_td_rep2_fract7_Proteoforms.psmtsv`
  (some entry points use `...\Task2-TopDownSearch\AllProteoforms.psmtsv`)
- An existing benchmark report: `D:\JurkatTopdown\deconvolution-benchmark.pdf`, which was generated
  from **centroided** simulated data.

## Branch lineage

| Branch | State | Notes |
|---|---|---|
| `DataSimulator` | 2024-06 → 2025-07, on `origin` | Predecessor. `MassSpectrometry/SimulatedData/` (EMG, empirically transformed Gaussian, Nelder-Mead peak fitting), grant figure. No surviving Claude sessions. |
| `SimulatorExperimental` | 2026-04-16, local | Early simulator plus export services. |
| `TopDownSimulator` | 2026-07-26, on `upstream` | Main line. Fitting, refit, peak-width model, feature ground truth, MS2. |
| `TopDownSimPolished` | 2026-08-03, local | `TopDownSimulator` plus the noise model, then a merge of `master`. |
| `TopDownSimClean` | 2026-10-06 | **Work here.** One squash commit on `master`. |

The two commits that only exist on `TopDownSimPolished`:

- `c6f5c0e27` "idk" (2026-07-31): the noise model (`Noise/`), `PrecursorMassShift`,
  `ScanListMerger`, the noise and jitter characterization tests, and the noisy and mass-shifted export
  entry points. It also commits two profiler traces, about 160 MB and 95 MB
  (`mzLib/analysisexample_*.nettrace.etlx`). Both are over GitHub's 100 MB limit, so no branch
  containing this commit can be pushed.
- `5bc54cb32` "lot of changes" (2026-08-03): a **merge of `master` as of 2026-07-22** (`220f67f2`).
  It did not change any simulator files.

`TopDownSimClean` keeps only the simulator code, its tests, `AnalysisExample.cs` and the
`TopDown-Simulator-*` docs. It leaves out the following, which remain on `TopDownSimPolished`:

- TopDownEngine (the feature finder), with its tests and docs
- the ProSightPD reader and the IMSP export services
- the deletion of `mzPlot`, a leftover from the macOS work
- the `Development/` plotting helpers
- the FLASHDeconv and TauriCS notes
- the agent and editor config files, and the traces

The simulator was retargeted to `net10.0` to match master.

## Running it

Build and test from `mzLib/`, not the repo root:

```
cd mzLib
dotnet build ./Test/Test.csproj
dotnet test  ./Test/Test.csproj --filter "FullyQualifiedName~TopDownSimulator"
```

The real-data entry points are `[Explicit]` tests in `mzLib/Test/FileReadingTests/AnalysisExample.cs`:

- `SimulateJurkatRep2AndWriteMzml`: fits models and writes centroided simulated mzML.
- `ExportRep2SliceAndQValueSimulations`: covers the 31–35 min slice and the full q ≤ 0.01 run.
- `ExportRep2NoisySimulationForInspection` and `ExportRep2FullNoisySimulation`: clean and noisy
  output. The full run is about 1 GB and takes several minutes.
- `ExportRep2MassShiftedSimulation` and `ExportRep2MassShiftedNoisySimulation`, plus
  `VerifyRep2MassShiftedSimulation`: MS1 simulated with shifted precursor masses, MS2 carried over
  from the raw file.

Behaviour is controlled with `MZLIB_TOPDOWN_SIM_*` environment variables, all read in
`AnalysisExample.cs`: `MODE`, `NO_DEDUP`, `NO_GLOBAL_ABUNDANCE_REFIT`, `GLOBAL_REFIT_MAX_MODELS`,
`CONSTANT_PEAK_WIDTH`, `MIN_SAMPLES_PER_SIGMA`, `PEAK_WIDTH_K`, `MASS_SHIFT_DA` and
`NOISE_DENSITY_SCALE`.

## Hard-won facts — do not re-derive

### The source raw is centroided, so σ_m cannot be fitted from it
Any peak-width estimator run on this file measures its own extraction window. That window is
`0.5/z`, which is proportional to m/z, so the fit confidently returns σ ∝ (m/z)^0.94 ± 0.08 instead
of the physical 1.5. In the measurement, 4084 of 4097 windows were too coarsely sampled to be peak
shapes. Guards now refuse this case. To reproduce the artifact, set `MIN_SAMPLES_PER_SIGMA=0`.
Fitting a width law needs a profile-mode acquisition. Otherwise set k directly: R = 60,000 at m/z 400
gives `PEAK_WIDTH_K ≈ 3.5e-7` for σ = k·(m/z)^1.5.

### Averagine drops the monoisotopic peak above ~31 kDa
`AverageResidue.MinRes = 1e-8` is an **absolute** probability cutoff. Above about 31 kDa the
all-light isotopologue falls below it and is removed. As a result, `IsotopeEnvelopeKernel.CentroidMzs(z)[0]`
is the lowest *surviving* isotopologue: +1.003 Da at 32 kDa, +2.007 at 40 kDa, +5.015 at 50 kDa.
Compute the monoisotopic m/z from the mass (`mass.ToMz(z)`). This bug once mislabeled the feature
ground truth for exactly the large proteoforms that matter most.

### The noise model is empirical; the analytic FT law is wrong by about 100×
The argument "noise is uniform in frequency, so density ∝ (m/z)^−1.5" predicts a gentle 4× falloff.
The measured density instead *rises* to a peak at m/z 700–800 and then falls about 500×, because the
instrument reports local maxima above a threshold rather than occupied frequency bins. Simulating
transients and running an FFT was considered and rejected: the simulator writes centroids, and mzLib
cannot read profile mzML. **Do not re-propose the FFT route.** To extend the model, measure first
with the `[Explicit]` tests in `Test/TopDownSimulator/NoiseCharacterization.cs`. Those tests read
`CentroidStream.Noises`/`.Baselines` directly because `ThermoRawFileReader` discards them.

### The post-elution noise is unstructured, not chemical background
After elution (50–55 min) there are about 16–17.6k centroids per scan, roughly twice the busy region,
because AGC extends injection time. Three measurements show it is random:

- 98.4 % of those peaks have no charge assigned.
- About 45 % have a +1.00335/z partner, the same for every z from 1 to 7, which is chance at this density.
- The brightest peaks persist scan-to-scan no more often than random ones (77 % in the busy region).

So the model uses independent per-scan singlets with no persistent contaminants. Do not add a
polymer or background generator.

### Peak jitter laws (`Noise/PeakJitterModel.cs`)
These were measured on rep2 fract7 against the top 400 IDs (`SignalJitterCharacterization.cs`):

- **m/z:** a per-scan common-mode offset with σ = 1.89 ppm, plus a per-peak term
  σ = √(9.06²/(S/N) + 0.81²) ppm.
- **Intensity:** σ of log intensity = √(1.026²/(S/N) + 0.165²). The 0.165 floor means even
  high-S/N peaks scatter about 17 % from scan to scan.
- **Reusable trick:** estimate intensity noise from the scan-to-scan scatter of the *ratio of adjacent
  isotopologues* of one species. The elution profile cancels exactly. The trick cannot see common-mode
  (AGC) fluctuation, which is why `CommonModeLogIntensitySigma` defaults to 0 as "unmeasured", not as
  "absent".

### Clean simulated files crash MetaMorpheus
Without injected noise, 2377 of 2906 survey scans in a simulated file are empty. MetaMorpheus
`GetMs2Scans` then throws `ArgumentNullException('source')` from LINQ, far from the cause. The
underlying `MzLibException` ("precursor scan contains no peaks") is swallowed in `_GetMs2Scans` and
leaves a null slot. To diagnose, count MS1 scans with `XArray.Length == 0`; the precursor references
are usually fine. To avoid it, use the noisy exports.

## Known defects (verified 2026-07-26)

Full list: `agent_info/TopDown-Simulator-Known-Issues.md`. The ones most likely to cause trouble:

1. `GlobalAbundanceRefitter.BuildTranspose` allocates outside `MaxBasisCacheBytes`, which roughly
   doubles peak memory, and there is no fallback. The uncached fallback branch has no test coverage.
2. The Gauss-Seidel update minimizes over each proteoform's own sample set, but the reported residual
   sums over all sets. Monotone descent has only been observed on the test fixture; it is not
   guaranteed. `ResidualFractionFallsMonotonicallyAndConvergesQuickly` is not a proof. Making the
   update globally exact was deliberately deprioritized.
3. `OverpredictedFraction` is not deduplicated, so it grows with crowding. Use it only as a within-run
   diagnostic.
4. `IsotopeEnvelopeKernel` is not thread-safe (its lazy caches are plain dictionaries), and
   `ForwardModel` holds kernels as instance state. Parallelizing `Rasterize` would corrupt them.

## Session history

- **Session `c8a948d9`** (2026-08-26 → 2026-09-04) was the only Claude session recorded on this
  branch. You asked for a centroided simulation of rep2 fract7 so the deconvolution benchmark report
  could be redone on it. You then withdrew the request: "Nevermind / The old one was centroided." That
  means the existing `deconvolution-benchmark.pdf` already uses centroided data. Nothing else was
  done. The transcript has since been deleted by the 30-day cleanup; only the prompts survive in
  `~/.claude/history.jsonl`.
- The substantive simulator work happened in earlier sessions (`e46ad809`, `6c92fbb1`, `154e71d2`,
  July 2026). Their transcripts are gone, but their findings are captured above and in the topic docs.

## Where the detail lives (`agent_info/`)

- `TopDown-Simulator-Plan.md`: the forward model, the phased plan and the locked decisions.
- `TopDown-Simulator-Status.md`: component status as of 2026-04-11 (outdated).
- `TopDown-Simulator-Forward-Model-Fidelity.md`: the σ_m-vs-m/z model, feature ground truth, residual
  dedup and the outcome of that work.
- `TopDown-Simulator-Known-Bugs.md`: the isotopologue ordering and stale-total refit bugs, both fixed.
  Also lists context "do not redo".
- `TopDown-Simulator-Known-Issues.md`: the deferred defects and the test-suite baseline.
- `TopDown-Simulator-Noise-Model.md`: noise measurements and model design. This file is untracked.
- `TopDown-Simulator-Performance-Incident.md`: the apparent hang in `AnalysisExample` and its fixes.
- `mzLib/TopDownSimulator/GlobalAbundanceRefitPlan.md`: the refit design.
