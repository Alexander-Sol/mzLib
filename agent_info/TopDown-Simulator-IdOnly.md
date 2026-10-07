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

## The `.v3` fits: the whole charge range

Everything above was measured on the `.v2` fits. Those extracted each species over its members'
precursor charges ± 2 and fitted μ_z, σ_z by moments over that window.

**Why the window was wrong.** Over charges 5–40, rep1's per-charge apex profiles show that the
envelopes are broad, with σ_z ≈ 3 for 14–17 kDa species. They end where the species leaves the scan
range (m/z < 600). The ± 2 window truncated them. Moments over a truncated window are pulled
towards the window centre, which sits on the MS2-selected charges, about one charge below the apex.

**What changed** (`ExportRep{1,2}FullNoisySimulation` with `OUTPUT_TAG=.v3`):

- **Window:** members ± 2 plus every charge at which the species lands inside the scan window
  (`GroundTruthExtractor.ObservableCharges`, `CHARGE_WINDOW=scan`). A sequence-predicted window
  (`CHARGE_WINDOW=sequence`) was tried first. With `.v2`'s narrow σ_z it was still too narrow, and its
  edges came from the predictor itself, which is circular.
- **Charge fit:** `ChargeDistributionFitter(trimFraction: 0.05)`. It keeps the contiguous run of
  charges around the most intense one that stays above 5 % of it and does not climb back into a
  neighbouring envelope. It then fits a parabola to log apex weighted by apex² (Caruana), which is
  exact for a truncated Gaussian.
- **Unexplained energy after the global refit:** rep2 0.45 → 0.359, rep1 0.278.

**Charge profiles against rep1's observed profiles** (`CompareChargeProfiles`; `.v2` numbers are
from the table above):

| Model | Median cos `.v2` → `.v3` | Most intense charge within 1 `.v2` → `.v3` | Offset median `.v3` |
|---|---|---|---|
| fitted to rep1 raw | 0.812 → **0.961** | 50 % → **83 %** | 0 |
| IDs, observed charges | 0.765 → 0.948 | 44 % → 50 % | −1 |
| IDs, sequence only (√mass + KRH − DE) | 0.842 → 0.970 | 58 % → 60 % | 0 |
| IDs, residue weights | – → 0.970 | – → 61 % | 0 |

The fitted model now beats every predictor on apex placement. Its cosine is about the same as the
sequence predictors' (0.961 against 0.970). The ID-only rows also gain from σ_z predicted from μ_z,
below; with bootstrap σ_z the sequence predictor scored 0.958.

**What identifications predict, re-measured on the `.v3` fits** (rep2, 436 species,
`IdOnlyPriorsCharacterization` with `FIT_TAG=.v3`):

| μ_z predictor | R² `.v2` → `.v3` | Residual sd `.v3` |
|---|---|---|
| mean member precursor charge | 0.939 → 0.36 | 2.35 |
| √mass | 0.50 → 0.64 | 1.77 |
| √mass + net basic (KRH − DE) | 0.62 → 0.70 | 1.61 |
| √mass + a weight per K, R, H, D, E | – → **0.72** | 1.57 |
| the same + mean member charge | – → 0.72 | 1.57 |
| the same + sequence length | – → 0.72 | 1.56 |

- **The precursor charge predicts little.** Its R² of 0.94 was an artefact of the window: the fitted μ_z
  was anchored to the members' charges. With the residue model it adds nothing.
- **Residue weights** (`ResidueChargeModel`, fitted on rep2):
  μ_z = −13.05 + 0.228·√M + 0.146 K + 0.218 R + 0.189 H − 0.048 D − 0.019 E.
  Arginine and histidine count more than lysine. The acidic residues count little.
- **σ_z grows with μ_z:** σ_z = −0.66 + 0.21·μ_z (R² 0.42; against mass, 0.19). `IdOnlyPriors` now
  predicts σ_z from the predicted μ_z instead of drawing it.
- **Abundance:** fitted over the whole envelope, it relates less closely to one precursor's
  intensity. Summed member intensity: R² 0.74 → 0.53. The brightest member's intensity scores 0.59,
  and `IdOnlyPriors` now uses it. Dividing by the predicted f(z) at the precursor's charge (0.60),
  adding mass, predicted σ_z or member count (≤ 0.59) do not help.
- **Retention time and shapes are unchanged:** apex offset −0.015 min; σ and τ as before.

**Held out on rep1, `.v3`** (`CompareHeldOut` with `FIT_TAG=.v3`, residue-weight charges):

| | `.v2` | `.v3` |
|---|---|---|
| μ_z error vs rep1's fit, residue weights / sequence / observed | – / 1.37 / 0.46 | 0.84 / 0.92 / 1.51 |
| log₁₀ abundance error, mean \|err\| (r) | 0.24 (0.83) | 0.34 (0.69) |
| whole run log TIC r, `train` template | 0.971 | 0.971 |
| cos at 41.1 min, fitted | 0.881 | 0.902 |
| cos at 41.1 min, IDs `self` / `train` | 0.678 / 0.668 | 0.580 / 0.557 |

At 41.1 min the bright H2A proteoforms (13.75–13.80 kDa) get μ_z within 0.2 and RT within 0.1–0.2 min.
Their abundances come out 0.2–0.7 dex low, and in the wrong order among near-isobaric proteoforms. A
regression with slope below one compresses the top of the range, and the `.v2` abundance R² was
inflated by the same window artefact as the charge. The 35.1 and 53.8 min scans are unchanged
(noise-dominated).

## Next

1. A signal-focused comparison: per species, real versus simulated envelope intensity around the
   apex, and peaks above S/N 10 only.
2. A held-out test across fractions (e.g. rep2 fract5 or fract6 IDs with a fract7 template), to test
   how far the template transfers. The MM114 searches under `Rep2_Raw` have fract5 and fract6
   identifications.
3. An unidentified-analyte component and an AGC model in place of the template (step 4 of the plan).
