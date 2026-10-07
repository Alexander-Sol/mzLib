# Top-Down Simulator — Simulating From Identifications Only

**Written:** 2026-10-07, on `TopDownSimClean`.
**Question:** given only a run's identifications, from the same instrument and method, how close is
the simulated run to the real one?

## Setup

- **Training run:** Jurkat rep2 fract7, models fitted to its raw data (`full.v2` export).
- **Held-out run:** rep1 fract7. Same fraction, same search
  (`Frac7_GPTMD_Search\Task2-TopDownSearch`), q ≤ 0.01.
  - 1113 identifications, grouped into 478 species.
  - Its raw file is read only to score the result.
- **Upper bound:** rep1 fitted to its own raw file (`AnalysisExample.ExportRep1FullNoisySimulation`).

The code:

- `TopDownSimulator/Extraction/SpeciesGrouper.cs`: species grouping, moved out of `AnalysisExample`.
- `TopDownSimulator/Prediction/IdOnlyPriors.cs`: learns the priors, saves and loads them as JSON, and
  predicts models.
- `Simulator.ReadGroundTruth`: reads a sidecar back into models.
- `MmResultLoader.LoadQualified`: loads one file's identifications at a q-value cutoff, carrying
  `PrecursorIntensity`.
- `Test/TopDownSimulator/IdOnlySimulation.cs`, run in order: `TrainPriors`, then
  `SimulateHeldOutFromIdsOnly` (with `MZLIB_TOPDOWN_SIM_IDONLY_TEMPLATE=train|self`), then
  `CompareHeldOut`.
- `Test/TopDownSimulator/IdOnlyPriorsCharacterization.cs`: the predictor measurements below.

## What identifications predict (rep2 fract7, 451 fitted species)

| Fitted parameter | Predictor | Fit |
|---|---|---|
| RT apex μ | anchor MS2 time | offset −0.014 min, IQR 0.12 min |
| log₁₀ abundance | log₁₀ Σ member precursor intensity | R² 0.74, residual 0.34 dex |
| charge μ_z | mean member precursor charge + K+R count | R² 0.946, residual 0.66 |
| RT σ, τ, charge σ_z | none found | joint bootstrap from the fitted population |

Charge predictors compared:

| Predictor | R² |
|---|---|
| mass | 0.48 |
| √mass | 0.50 |
| K+R | 0.30 |
| K+R+H | 0.34 |
| mass + K+R | 0.58 |
| mean member charge | 0.939 |
| mean member charge + K+R | 0.946 |

Basic residues add information beyond mass, but the observed precursor charges carry almost all of it.
Charge σ_z depends on neither mass (R² 0.015) nor K+R (0.03). Half the fitted τ are 0. RT σ does not
depend on RT.

Regressions are used without adding the residual back: the residual is the error of not having raw
data, not real variation between species.

## The scan grid and noise template

Both come from a template run on the same method. The scan grid is its MS1 times; the noise comes from
`ScanNoiseConditions` on its scans.

- `train` uses rep2, which is the honest ID-only setting.
- `self` uses rep1's own scans. It leaks the held-out noise, so it separates errors in the signal model
  from errors in the template.

## Results on rep1 fract7

Predicted against rep1's own fitted parameters (453 species):

- RT μ error: mean |0.067| min, which is about 0.6 σ;
- log₁₀ abundance error: mean |0.24| dex, r 0.83;
- μ_z error: mean |0.46| charges.

Whole run, simulation interpolated onto the real scan times:

| Simulation | log TIC r | peaks/scan r | TIC ratio | peak ratio |
|---|---|---|---|---|
| fitted to rep1 raw (upper bound) | 0.966 | 1.000 | 0.97 | 0.98 |
| IDs only, `self` template | 0.966 | 1.000 | 0.94 | 0.98 |
| IDs only, `train` template | 0.971 | 0.990 | 0.94 | 0.97 |

Per scan, m/z 600–2000:

| RT | Regime | Simulation | Peaks (real) | cos | √cos | TIC | KS |
|---|---|---|---|---|---|---|---|
| 41.1 | busy, histones | fitted | 7366 (7277) | 0.881 | 0.623 | 1.11 | 0.14 |
| 41.1 | busy, histones | IDs, `self` | 7212 | 0.678 | 0.488 | 0.89 | 0.13 |
| 41.1 | busy, histones | IDs, `train` | 8016 | 0.668 | 0.488 | 0.88 | 0.20 |
| 35.1 | busy | fitted / IDs `self` / IDs `train` | 11749 / 11749 / 11408 (11950) | 0.197 / 0.197 / 0.195 | ~0.30 | 1.03–1.10 | 0.13–0.15 |
| 53.8 | wash | all three | ~16.4k–17.0k (16534) | ~0.21 | ~0.30 | 1.04–1.10 | 0.51–0.65 |

## Reading it

- **The signal model is where ID-only loses.** At the signal-dominated scan (41.1) cosine drops from
  0.88 to 0.67. `self` and `train` agree, so the template is not the cause.
- **The template transfers well within a fraction.** Rep2's grid and noise reproduce rep1's
  peaks-per-scan trace (r 0.99) and TIC trace (r 0.97). This is untested across fractions.
- **Whole-scan metrics in noise-dominated scans do not discriminate.** At 35.1 and 53.8 the fitted and
  ID-only simulations score identically, because the conditioned noise dominates. A signal-focused
  metric, such as peaks above S/N 10 or per-species apex envelopes, is needed to judge the signal model
  there.

## Without precursor charges

`ChargePredictor.SequenceOnly` predicts μ_z from mass and sequence alone. To use it in
`SimulateHeldOutFromIdsOnly`, set `MZLIB_TOPDOWN_SIM_IDONLY_CHARGE=sequence`. The simulated charge range
then comes from the predicted μ_z ± 3σ_z.

On rep2 fract7, the best sequence-only predictor is √mass plus net basic residues (K+R+H − D+E):
R² 0.62, residual 1.74. The alternatives scored lower:

| Predictor | R² |
|---|---|
| √mass + K+R | 0.59 |
| √mass + D+E | 0.59 |
| mass + K+R | 0.58 |
| √mass alone | 0.50 |
| K+R alone | 0.30 |

On rep1, the error against rep1's fitted μ_z is mean |1.37| charges, against 0.46 with precursor charges.
The whole-scan metrics do not change, because they cannot see the charge distribution.

`CompareChargeProfiles` compares each species' observed per-charge profile, summed over the anchor
RT ± 0.25 min across charges 5–40 of rep1's raw data, with each model's f(z):

| Model | Median cos | p25 cos | Most intense charge within 1 | Offset p10 / median / p90 |
|---|---|---|---|---|
| fitted to rep1 raw | 0.812 | 0.699 | 50 % | −3 / −1 / +2 |
| IDs, observed charges | 0.765 | 0.665 | 44 % | −4 / −1 / +3 |
| IDs, sequence only | 0.842 | 0.756 | 58 % | −2 / 0 / +1 |
| (mean member precursor charge) | – | – | – | −4 / −1 / +3 |

**The precursor charges chosen for MS2 sit about one below the envelope's most intense charge.** The
fitter only extracts members' charges ± 2, so the fitted μ_z inherits that selection bias, and so does
every predictor trained on the fits. The sequence-only predictor is the least biased against the
observed profiles.

Read the comparison with care. The observed profile includes random peaks and overlapping species
inside the extraction windows, the same for every model, so it is a noisy reference. Still, it suggests
the fitter's charge window is too narrow.

## Next

0. Widen the fitter's extraction charge range, for example to the sequence-predicted μ_z ± 3σ_z combined
   with members ± 2. Then refit, retrain the priors and rerun `CompareChargeProfiles`. This removes the
   MS2-selection bias from the fitted charge centres.

1. A signal-focused comparison: per species, real versus simulated envelope intensity around the
   apex, and peaks above S/N 10 only.
2. A held-out test across fractions (e.g. rep2 fract5 or fract6 IDs with a fract7 template), to test
   how far the template transfers. The MM114 searches under `Rep2_Raw` have fract5 and fract6
   identifications.
3. An unidentified-analyte component and an AGC model in place of the template (step 4 of the plan).
