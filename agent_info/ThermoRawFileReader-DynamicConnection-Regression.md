# ThermoRawFileReader — `GetOneBasedScanFromDynamicConnection` Regression

**Added:** 2026-08-05
**Code:** `mzLib/Readers/Thermo/ThermoRawFileReader.cs` (method at ~line 162)
**Breaking commit:** `732c0033` — "Raw file connection not closing (#709)", Edwin Laboy, 2023-09-22
**Last good commit:** `e4c48c6a` (2023-09-22)

## Symptom

`GetOneBasedScanFromDynamicConnection` opens a brand-new Thermo file accessor on **every call**
instead of reusing the connection established by `InitiateDynamicConnection()`. The accessor is
never disposed, so each call leaks a native file handle. In a scan-by-scan loop over a run this is
one `RawFileReaderAdapter.FileFactory` open per scan.

## The bug

Current shape:

```csharp
public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
{
    if (CheckIfScansLoaded() && oneBasedScanNumber <= Scans.Length)
        return GetOneBasedScan(oneBasedScanNumber);

    var dymConnection = RawFileReaderAdapter.FileFactory(FilePath);   // <-- new connection, every call
    dymConnection.SelectInstrument(Device.MS, 1);

    if (dymConnection == null)                                        // <-- dead code; just assigned
    {
        throw new MzLibException("The dynamic connection has not been created yet!");
    }

    if (oneBasedScanNumber > dymConnection.RunHeaderEx.LastSpectrum ||
        oneBasedScanNumber < dymConnection.RunHeaderEx.FirstSpectrum)
    {
        return null;                                                  // <-- leaks on this path too
    }

    return ThermoRawFileReader.GetOneBasedScan(dymConnection, filterParams, oneBasedScanNumber);
}
```

Two consequences beyond the leak:

- The `null` guard can never fire — it sits below the assignment, so the "dynamic connection has not
  been created yet" error is unreachable. Callers who skip `InitiateDynamicConnection()` get silent
  success rather than the intended exception.
- The private `dynamicConnection` field is now only touched by `InitiateDynamicConnection()`,
  `CloseDynamicConnection()`, and `GetMsOrderByScanInDynamicConnection()`. The read path and the
  lifecycle path have drifted apart: `CloseDynamicConnection()` disposes a connection that
  `GetOneBasedScanFromDynamicConnection` never used.

## Timeline

| Commit | Date | State |
|---|---|---|
| `c3b032a2` non-static thermo dynamic connection (#498) | — | ✅ correct, in `ThermoRawFileReader/ThermoDynamicData.cs` |
| `d25500cc` beginning refactoring | 2022-09-27 | ✅ `Readers` project created; logic copied verbatim |
| `fc86427c` Reader refactoring, take 2 (#698) | 2023-04-10 | ✅ still correct in `Readers/Thermo/ThermoRawFileReader.cs` |
| `85e9e426` Prepare Release — Readers backward compatibility | 2023-04-26 | ✅ |
| **`e4c48c6a`** | 2023-09-22 | ✅ **last good commit** |
| **`732c0033`** Raw file connection not closing (#709) | 2023-09-22 | ❌ **regression introduced** |
| `6c77ca01` scan descriptions from Raw Files (#736) | 2023-11-27 | broken |
| `e33a478e` Tims tof reader (#812) | 2025-01-19 | broken |
| `8cea9fcc` Eliminate Need to call `Loaders.LoadElements` (#833) | 2025-02-18 | broken |
| `bee61eb8` Made Dynamic reading thread safe (#936) | 2025-08-21 | field → local `dymConnection`; behavior unchanged |
| `c8b2bf9c` Check if statically loaded before dynamic access (#948) | 2025-09-04 | added `CheckIfScansLoaded()` early return — sidesteps the leak only when data is already statically loaded |

**Note:** the Readers refactor is *not* the culprit. It preserved the original behavior intact.
The break landed ~5 months later, in a PR whose stated purpose was to fix connection leaks.

## Reference implementation (pre-break)

`git show fc86427c:mzLib/Readers/Thermo/ThermoRawFileReader.cs`

```csharp
public override MsDataScan GetOneBasedScanFromDynamicConnection(int oneBasedScanNumber, IFilteringParams filterParams = null)
{
    if (dynamicConnection == null)
        throw new MzLibException("The dynamic connection has not been created yet!");

    if (oneBasedScanNumber > dynamicConnection.RunHeaderEx.LastSpectrum
        || oneBasedScanNumber < dynamicConnection.RunHeaderEx.FirstSpectrum)
        return null;

    return ThermoRawFileReader.GetOneBasedScan(dynamicConnection, filterParams, oneBasedScanNumber);
}
```

This shape goes back unchanged to `ThermoDynamicData.cs` in `c3b032a2`.

## Caveat on fixing it

A straight revert to the `fc86427c` shape reintroduces the race that `bee61eb8` was fixing:
`IRawDataPlus` instances are not safe to share across threads, and the per-call `FileFactory` was
the (costly) way that PR bought thread safety.

The fix needs to preserve both properties — one connection per thread, created once, disposed on
close. Two viable routes:

1. **`ThreadLocal<IRawDataPlus>`** backing the connection, disposed in `CloseDynamicConnection()`.
2. **Thermo's own thread manager** — `RawFileReaderFactory.CreateThreadManager(FilePath)` plus
   `CreateThreadAccessor()` per thread. This pattern is already used in `LoadAllStaticData` in the
   same file, so it is consistent with existing code.

Either way, restore the `null` guard to a position where it can actually fire, and make the
out-of-range `return null` path not leak.

## Useful commands

```powershell
# trace the method's history
git log --oneline --date=short --format='%h %ad %s' `
  -L '/public override MsDataScan GetOneBasedScanFromDynamicConnection/,/^        }/:mzLib/Readers/Thermo/ThermoRawFileReader.cs'

# the diff that broke it
git show 732c0033 -- mzLib/Readers/Thermo/ThermoRawFileReader.cs

# the last correct version
git show fc86427c:mzLib/Readers/Thermo/ThermoRawFileReader.cs
```
