# ThermoRawFileReader.LoadAllStaticData — loading performance

Benchmarked on a 978 MB / 48,017-scan Thermo .raw file (`HFX_MB_14751_19.raw`),
Xeon 4214R (12c/24t), .NET 8, Release, warm OS file cache, best of 3 reps.

## TL;DR

The dominant cost of `LoadAllStaticData` was **not** the parallel scan loop — it was
`GetSourceFile()`, which computes a **SHA-1 over the entire file** (~3.0 s, ~55% of the
load) and previously ran serially *after* the scan loop. The fix is to run that hash on
a background task so it overlaps the scan reading. This makes every thread-count
configuration faster (biggest win at low thread counts) and, as a side effect, makes
sequential reading optimal.

## Correcting an earlier wrong conclusion

An earlier version of this note claimed the parallel loop was "the worst of every
option" and that per-chunk Thermo accessor creation ("accessor churn") cost ~3 s. That
was a measurement error: it compared a loop-only harness against the full
`LoadAllStaticData` (which also hashes the file) and attributed the difference to the
loop. Isolating the components disproves it — the old and new loops are within noise of
each other, and the ~3 s gap was the SHA-1, not accessor creation.

## Component breakdown (the corrected measurements)

| Component | @1 thread | @4 threads |
| --- | --: | --: |
| SHA-1 over the whole file (`GetSourceFile`) | 3060 ms | 3060 ms (fixed, does not parallelize) |
| Scan loop — old code (Partitioner, per-chunk accessor), loop only | 3841 | 2371 |
| Scan loop — contiguous blocks, per-worker accessor, loop only | 3797 | 2254 |
| Full `LoadAllStaticData` (loop + SHA-1, serial) | 6908 | 5387 |

The old Partitioner loop parallelized fine (~3.8 → 2.4 s, ~1.6×, plateauing at ~4
threads — the Thermo reader serializes internally beyond that). Rewriting the loop to
contiguous per-worker blocks changed nothing measurable (2371 vs 2254 ms), so that
rewrite was dropped.

## The fix that shipped: overlap the SHA-1 with the scan loop

`GetSourceFile` does not touch the Thermo reader (it opens its own `FileStream`), so it
can run concurrently with the scan-reading loop on a background task:

```csharp
var sourceFileTask = Task.Run(GetSourceFile);   // hashes the file while scans are read
// ... existing Parallel.ForEach scan loop ...
SourceFile = sourceFileTask.GetAwaiter().GetResult();
```

Full `LoadAllStaticData`, before vs after the overlap:

| maxThreads | before | after | note |
| --- | --: | --: | --- |
| 1 (sequential) | 6908 | 4299 | ~38% faster |
| 2 | — | 4250 | fastest |
| 4 | 5387 | 4370 | ~19% faster |
| 8 | — | 4871 | slower than sequential |
| 24 | ~5600 | 5142 | slower than sequential |

The overlap is a strict improvement at every thread count and needs no behavioural
change (the `SourceFile` field ends up identical).

## Scan reading is now sequential

Once the SHA-1 is overlapped, **sequential reading is the fastest configuration**, and
higher thread counts only get slower because the loop workers contend with the hash task
and with each other (the reader serializes internally and already caps out near 4
threads). So the scan loop was changed to a plain sequential read over a single accessor:
one thread reads the scans while one thread hashes the file. `maxThreads` is retained on
the signature for API compatibility but no longer affects reading.

Full `LoadAllStaticData` is now flat at ~4.2 s across all `maxThreads` values (identical
peak count), versus the original ~6.9 s single-threaded / ~5.4 s at four threads.

## Validation across file sizes and acquisition types

Old (parallel loop + serial SHA-1) vs new (sequential read + overlapped SHA-1), full
`LoadAllStaticData`, best of 3, warm cache. "new" was flat between 1 and 4 threads
(within ~1%), so a single value is shown.

| File | Size | Scans | SHA-1 | old @1 | old @4 | new | vs old @4 |
| --- | --: | --: | --: | --: | --: | --: | --: |
| ALS Biogen (bottom-up) | 339 MB | 24,721 | 1.0 s | 2643 | 2053 | ~1500 | 1.38× |
| HFX (bottom-up) | 978 MB | 48,017 | 3.0 s | 6848 | 5431 | ~4000 | 1.33× |
| Jurkat top-down | 1562 MB | 5,327 | 4.8 s | 7863 | 6709 | ~5600 | 1.18× |

The new code is faster on every file at every thread count — 18–38% vs the 4-thread
config, 38–77% vs single-threaded. The overlap behaves as `new ≈ max(sequential read,
SHA-1)`: on ALS the read dominates and the 1 s hash disappears under it; on Jurkat the
4.8 s hash dominates and the read disappears under it.

Note the trend: the larger the file, the more the whole-file SHA-1 dominates, and the
smaller the *relative* win (Jurkat's 4.8 s hash is an unavoidable floor once the read is
already hidden under it). This is exactly the case a **lazy `GetSourceFile`** (hash only
when `SourceFile` is actually read) would help most — deferred, since it changes when the
field is populated and needs a consumer audit.

## Reproducing

`Development/FileReading/BenchmarkFileReading.cs` — an `[Explicit]` NUnit benchmark that
sweeps `LoadAllStaticData` across thread counts for a local file (any supported format).
Point it at a file via the `MZLIB_BENCHMARK_FILE` environment variable or by editing
`DefaultBenchmarkFile`, then:

```
dotnet test mzLib/Development/Development.csproj -c Release --filter "FullyQualifiedName~Benchmark_FileReading" -l "console;verbosity=detailed"
```

It reports best/median per thread count with a speedup column relative to the
single-threaded (sequential) baseline, and asserts every configuration reads an
identical peak count.
