# FlashLFQ → Rust Rewrite: Feasibility Analysis & Plan

**Date:** 2026-06-25
**Author:** Claude (Opus 4.8) + Alex
**Status:** Feasibility approved; scope decided; no code written yet.

## Goal

Rewrite the core peak-tracing functionality of mzLib's **FlashLFQ** and **Readers**
projects into **Rust**, exposed to **Python** via bindings. The pipeline:

1. Read identification results (starting with `.psmtsv`).
2. Read raw data files (`.mzML` first; `.raw` later).
3. Detect/trace isotopic-envelope peak traces (XICs) in the data that correspond to
   identified peptides, and integrate them to intensities.

Most logic stays behind the scenes in Rust. Python is for **driving** the pipeline and
**visualizing** detected peak traces.

## Verdict

**Clearly feasible.** The hardest, most tedious part — reading vendor data files — is
already solved by a mature Rust mass-spec ecosystem, and the core peak-tracing algorithm
is conceptually simple. With the scope below, nothing is a research-grade unknown. This is
a favorable first serious Rust project: small, algorithmic, with the gnarly parts (vendor
I/O, ML) delegated to mature crates and to Python.

Main residual risks:
1. Isotope-distribution **numeric parity** (Rust crate vs mzLib's `IsotopicDistribution`).
2. Discipline on building a **golden-file parity harness** early (mzLib vs Rust, peptide-by-peptide).

## Scope decisions (agreed)

| Area | Decision |
|---|---|
| **Core XIC quant** | In scope, first. psmtsv + mzML → isotope-envelope peak tracing + integration. |
| **MBR (match-between-runs)** | Later. RT alignment + acceptor-peak search in Rust; **ML/PEP model in Python**. |
| **IsoTracker** | Out of scope (for now). |
| **Bayesian protein quant** | Out of scope (for now). |
| **Thermo `.raw`** | Deferred. mzML-only MVP; add `.raw` later via the `mzdata` crate (reader swap only — IndexingEngine is format-agnostic). |
| **Serialization** | Move off custom NetSerializer blob → **Parquet** (or SQLite only if random-access query pattern emerges). Worth doing in mzLib C# too, independently. |

### Why ML-in-Python for MBR

This sidesteps the weakest part of the Rust ecosystem. MBR splits cleanly along a
Rust/Python boundary:

- **Rust half (heavy I/O + search):** RT alignment between runs (`GetRtCalSpline`),
  donor-file selection, acceptor-peak search (`FindAllAcceptorPeaks`) — all reuse the
  IndexingEngine + XIC machinery already built for core quant. Produces a **feature table**
  (one row per candidate transferred peak).
- **Python half (the model):** the PEP step is currently a FastTree gradient-boosted model
  in `Microsoft.ML`. In Python that's a few lines of `lightgbm`/`xgboost`/`scikit-learn`
  on the feature table, returning scores Rust folds back into FDR.

MBR runs **once at the end** on an already-built table, so the Rust→Python→Rust hand-off
is negligible. You get a mature, debuggable, retrainable model and avoid Rust GBM immaturity.

## Architecture

```
              ┌────────────────────────── Python ──────────────────────────┐
              │  driving + visualization + MBR model (lightgbm/sklearn)     │
              └───────────────▲───────────────────────────▲─────────────────┘
                              │ PyO3 / maturin             │ Arrow / NumPy
              ┌───────────────┴───────────────────────────┴─────────────────┐
              │                          Rust                                │
              │  psmtsv parse → isotope dist → binned IndexingEngine →       │
              │  XIC tracing/integration (GetIsotopicEnvelopes + CutPeak) →  │
              │  FDR → Parquet output                                        │
              │  (later) MBR: RT align + acceptor search → feature table     │
              └───────────────▲──────────────────────────────────────────────┘
                              │ mzdata crate
                       .mzML  │  (.raw later)
```

### Python binding surfaces (the public API)

1. **`quant(results_path, raw_paths, params) -> peaks`**
   Returns either an in-memory peak table or a **path to a Parquet/SQLite store**.
   For full runs with many files, return a **Parquet path** (cheap to hand back,
   lazy-loadable in polars/pandas). Keep an option to materialize small result sets in memory.

2. **IndexingEngine calls — the visualization workhorse.**
   Expose: build/load an index from a data file, and query it. The key methods:
   - `GetIndexedPeak(mz, scanIndex, tolerance)` (point query)
   - an **XIC extractor**: "give me (RT, intensity) across scans for this m/z + charge +
     tolerance" → return as **two `f64` NumPy arrays** (via the `numpy` PyO3 crate), the
     natural shape for matplotlib/plotly.
   A detected peak isn't just an intensity number — it's the underlying XIC array plus the
   integration bounds (`CutPeak` gives apex + boundaries), which you'll want to overlay when plotting.

3. **MBR** — its own call (or a flag on `quant`), with the Python-side model hook described above.

### Data interchange

Use **Arrow as the interchange spine**: Rust `arrow`/`parquet` ↔ Python `pyarrow`/`polars`
is zero-copy and avoids hand-rolling PyO3 conversions per struct.
- Big result tables → **Parquet** on disk.
- Small interactive things (one peptide's XIC for a plot) → **NumPy** arrays.

## Rust ecosystem leverage (don't start from zero)

- **`mzdata`** — mature reader for mzML, mzMLb, MGF, and **Thermo `.raw`** (bridges the
  Thermo .NET RawFileReader / ThermoRawFileParser). Eliminates the ~4,800-line mzML port
  *and* gives the eventual `.raw` path. Biggest de-risking factor.
- **`sage`** (Michael Lazear's proteomics search engine) — production Rust codebase that
  already does XIC-based label-free quant with the **same** binned-index + isotope-envelope
  approach. Strong reference implementation.
- **`rustyms`** — peptide/proteoform chemistry: chemical formulas from sequence+mods and
  isotopic distribution calculation. Can replace much of the `Chemistry`/`Proteomics` port.
- **`timsrust`** — native Bruker timsTOF `.d` reader, if ever needed.
- **PyO3 + maturin** — standard, mature Python-bindings path. Low risk.

## What the Rust port actually owns (with this scope)

**Port faithfully from mzLib:**
- The binned m/z `IndexingEngine<T>` — `mzLib/MassSpectrometry/PeakIndexing/IndexingEngine.cs`
  (~340 lines). Jagged array bucketed by `m/z × BinsPerDalton`, binary search over scan
  index within each bin. Trivial to port.
- `FlashLfqEngine.GetIsotopicEnvelopes` + `CheckIsotopicEnvelopeCorrelation` + `CutPeak`
  — the XIC tracing / envelope-correlation / peak-integration core. **The algorithmic heart.**
- `.psmtsv` parsing — `Readers/InternalResults/IndividualResultRecords/PsmFromTsv.cs`
  (~144 lines) + `MzLibExtensions.MakeIdentifications` (~66 lines). Trivial with `csv` + `serde`.
- Charge-state aggregation, FDR.

**Borrow from crates (and validate numeric agreement vs mzLib):**
- mzML reading → `mzdata`.
- Isotope distributions + peptide→formula → `rustyms`.
  (mzLib reference: `Chemistry.IsotopicDistribution.GetDistribution`,
  `Proteomics.AminoAcidPolymer.Peptide.GetChemicalFormula`, `PeriodicTable`, `ChemicalFormula`.)

**Defer:**
- MBR Rust-side feature generation until core quant has parity.
- The MBR model lives in Python from the start.

## Key mzLib reference points (for the port)

| Concept | C# location |
|---|---|
| Engine entry / orchestration | `FlashLFQ/FlashLfqEngine.cs` (`Run`, ~1966 lines total) |
| Theoretical isotope dists | `FlashLfqEngine.CalculateTheoreticalIsotopeDistributions` (~110 lines) |
| XIC tracing core | `FlashLfqEngine.GetIsotopicEnvelopes`, `CheckIsotopicEnvelopeCorrelation`, `CutPeak` |
| Binned index + query | `MassSpectrometry/PeakIndexing/IndexingEngine.cs` (`GetIndexedPeak`, `IndexPeaks`) |
| FlashLFQ index wrapper | `FlashLFQ/PeakIndexingEngine/PeakIndexingEngine.cs` (note: `Serialize/DeserializeIndex` use NetSerializer → replace with Parquet) |
| Identification model | `FlashLFQ/Identification.cs` (BaseSequence, ModifiedSequence, MonoisotopicMass, charge, RT, file, PeakfindingMass) |
| psmtsv → Identification | `FlashLFQ/ResultsReading/MzLibExtensions.cs` (`MakeIdentifications`, `MakeSpectraFileDict`) |
| psmtsv record / reader | `Readers/InternalResults/IndividualResultRecords/PsmFromTsv.cs`, `.../ResultFiles/PsmFromTsvFile.cs` |
| Data file dispatch | `Readers/MsDataFileReader.cs` (`GetDataFile`) → mzML reader in `Readers/MzML/*` |
| MBR (later) | `FlashLfqEngine.GetRtCalSpline`, `FindPeptideDonorFiles`, `FindAllAcceptorPeaks`, `FindIndividualAcceptorPeak`; PEP in `FlashLFQ/MBR/PEP/PepAnalysisEngine.cs` (Microsoft.ML FastTree) |

### External C# dependencies (context for what NOT to port)
- FlashLFQ.csproj: CsvHelper, MathNet.Numerics, **Microsoft.ML + Microsoft.ML.FastTree**
  (MBR/PEP only → Python), NetSerializer (→ Parquet), SharpLearning.Optimization.
- Readers.csproj: CsvHelper, OpenMcdf, System.Data.SQLite, ZstdSharp, **Thermo vendor DLLs**
  (closed-source .NET — the reason `.raw` is deferred to `mzdata`), Bruker native DLLs.

## Serialization change (Parquet)

Two separate things currently use NetSerializer:
- **The peak index** (`PeakIndexingEngine.SerializeIndex`, the `.ind` blob) — columnar
  (m/z, intensity, scan, RT) compresses extremely well → **Parquet** (Rust `arrow` + `parquet`).
- **The results** (peaks per peptide, intensities per file) — natural tabular Parquet output,
  directly readable from Python pandas/polars.

Use **SQLite only** if a *random-access query* pattern emerges ("all peaks for peptide X")
rather than bulk load. FlashLFQ's "load all → process → write all" pattern favors Parquet.
This change is worth making in mzLib's C# too, independent of the Rust rewrite.

## Phased plan

- **Phase 0 — spike (days):** Rust crate + maturin/PyO3 skeleton. Read an `.mzML` via
  `mzdata`, expose MS1 scan counts to Python. Validates the whole toolchain cheaply.
- **Phase 1 — MVP (core ask):** psmtsv parse → isotope dists (`rustyms` or port) → binned
  index (port) → XIC tracing/integration (port `GetIsotopicEnvelopes`/`CutPeak`) → Parquet
  output. **mzML only.** Build the golden-file parity harness here. Expose the three Python
  surfaces (Quant, IndexingEngine/XIC, stub MBR).
- **Phase 2:** Thermo `.raw` via `mzdata`; charge aggregation, FDR hardening.
- **Phase 3:** MBR — Rust RT alignment + acceptor search → feature table → Python ML model → FDR.

## Build the parity harness early

The thing that consumes time isn't writing Rust — it's verifying numeric parity. Run mzLib
and the Rust impl on identical inputs and diff peptide-by-peptide. Subtle off-by-one
isotope-index or tolerance differences are the usual failure mode. Make this a first-class
test from Phase 1.
