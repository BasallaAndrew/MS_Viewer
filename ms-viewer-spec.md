# Lightweight Thermo .raw Viewer for macOS — Spec

## Goal

A fast, minimal desktop app for browsing Thermo `.raw` LC-MS files on a Mac. It should feel like a stripped-down Xcalibur Qual Browser: open a file, see the chromatogram, click around, and have the spectra update instantly. No conversion step, no processing pipeline, no project management.

## Layout

Three stacked panes sharing one window, with a thin toolbar/status bar:

```
┌──────────────────────────────────────────────┐
│ Toolbar: Open | TIC/BPC toggle | file name   │
├──────────────────────────────────────────────┤
│ 1. Chromatogram (TIC or base peak) vs RT     │
│    - vertical line = selected scan           │
│    - shaded region = selected RT range       │
├──────────────────────────────────────────────┤
│ 2. MS1 spectrum (single scan or averaged)    │
│    - title: scan # / RT, or "avg of N scans, │
│      RT a–b"                                  │
│    - precursor m/z of current MS2 marked     │
├──────────────────────────────────────────────┤
│ 3. MS2 spectrum                              │
│    - header: "MS2 3 of 17 | scan 12345 |     │
│      RT 23.41 | precursor 745.38 (z=2) |     │
│      HCD 28" + ◀ ▶ buttons                    │
├──────────────────────────────────────────────┤
│ Status bar: cursor m/z / intensity / RT      │
└──────────────────────────────────────────────┘
```

Panes should be resizable (splitter). Default heights roughly 30 / 35 / 35.

## Interactions

### Chromatogram pane
- **Click**: select the nearest MS1 scan at that RT. MS1 pane shows that scan. The MS2 list becomes the MS2 scans triggered by that MS1 (i.e., between this MS1 and the next MS1).
- **Click-and-drag** (horizontal): select an RT range. MS1 pane shows the **average** of all MS1 scans in the range. The MS2 list becomes all MS2 scans in that RT range, ordered by scan number.
- Use a modifier to distinguish zoom from select. Suggested: plain drag = select range; scroll wheel = zoom x; shift+scroll = zoom y; right-drag or middle-drag = pan; double-click = reset zoom.
- **Left/Right arrow keys** (when chromatogram has focus): step to previous/next MS1 scan (single-scan mode).

### MS2 stepping
- **◀ / ▶ buttons** and **Up/Down arrow keys** step through the current MS2 list.
- When the MS2 changes, mark its precursor m/z and isolation window on the MS1 pane (a vertical line plus a shaded band).
- Optional nice-to-have: clicking a peak in the MS1 jumps to the MS2 in the current list whose precursor is closest to that m/z (within the isolation window). If none, do nothing and flash the status bar.

### Spectrum panes (MS1 and MS2)
- Scroll = zoom m/z; shift+scroll = zoom intensity; drag = pan; double-click = reset.
- Y-axis auto-rescales to the tallest peak in the visible m/z window.
- Hover shows m/z and intensity of nearest peak in the status bar.
- Label the top ~10 peaks in the visible window with m/z (4 decimals), avoiding overlapping labels if easy.

### General
- Open via File → Open, Cmd+O, drag-and-drop onto the window, or a command-line argument (`msview file.raw`).
- Right-click on any pane: "Copy image" and "Export spectrum as CSV" (m/z, intensity).

## Tech stack

- **Python 3.11+**
- **PySide6** (Qt) for the window and layout
- **pyqtgraph** for plotting — chosen for speed with large arrays; avoid matplotlib for interactive panes
- **pythonnet** (CoreCLR runtime) to call Thermo's **RawFileReader** .NET assemblies directly
- **.NET 8 runtime** (arm64 on Apple Silicon)
- **NumPy** for averaging/binning

### Reading .raw files

Thermo publishes RawFileReader as cross-platform .NET assemblies: https://github.com/thermofisherlsms/RawFileReader. Its license does not allow redistributing the DLLs, so the app should look for them in a configurable folder (e.g. `~/.msview/rawfilereader/`) and show a clear error with the download link if they're missing. Add a `setup` note in the README on where to put them.

Rough API sketch (verify against the RawFileReader example code in that repo):

- `RawFileReaderAdapter.FileFactory(path)` → raw file object; check `IsError`.
- `raw.SelectInstrument(Device.MS, 1)`.
- `raw.RunHeaderEx.FirstSpectrum` / `LastSpectrum` for scan range.
- Per scan: `RetentionTimeFromScanNumber`, `GetFilterForScanNumber` (MS order, polarity, analyzer), `GetScanEventForScanNumber(...).GetReaction(0)` for precursor mass, isolation width, activation type, collision energy.
- Trailer extra data (`GetTrailerExtraInformation`) for "Monoisotopic M/Z" and "Charge State" — prefer these over the filter precursor when present.
- Spectra: use the centroid stream (`GetCentroidStream`) when available (Orbitrap data); otherwise `GetSegmentedScanFromScanNumber` for profile data. Draw centroids as sticks, profile as a line.
- Chromatogram: `GetChromatogramData` with `ChromatogramTraceSettings(TraceType.TIC)` or base peak, filtered to MS1.
- RawFileReader has built-in scan-averaging helpers; use them if they work cleanly through pythonnet, otherwise implement averaging in NumPy (below).

### Fallback path
If pythonnet/RawFileReader proves painful, keep the data layer behind an interface (`RawSource` with `scan_index()`, `chromatogram()`, `spectrum(scan)`, `average(scans)`) and add an mzML backend using `pyteomics` or `pymzml`. Users can convert with ThermoRawFileParser (https://github.com/compomics/ThermoRawFileParser). The GUI should not care which backend is in use.

## Data model

On open, build a **scan index** once (a NumPy structured array or small dataframe):

| field | notes |
|---|---|
| scan | scan number |
| rt | minutes |
| ms_level | 1, 2, ... |
| precursor_mz | NaN for MS1 |
| charge | 0 if unknown |
| isolation_width | |
| activation | HCD / CID / ETD ... |
| collision_energy | |
| parent_ms1 | scan number of the preceding MS1 (for DDA grouping) |
| tic, base_peak | per scan |

Spectra are loaded lazily on demand and cached with an LRU cache (e.g. last ~500 scans).

## Averaging

- For centroided data: pool all peaks from the selected scans, bin by m/z with a ppm-based tolerance (default 5 ppm, configurable), sum intensities within a bin, divide by number of scans, and report the intensity-weighted mean m/z per bin.
- For profile data: interpolate each scan onto a common m/z grid and average.
- Cap averaging at a sensible number of scans (e.g. 500) with a status bar warning, so an accidental whole-run drag doesn't hang the UI. Ideally run averaging in a worker thread.

## Performance targets

- Opening a ~1–2 GB Orbitrap DDA file and showing the chromatogram: a few seconds at most. Build the scan index in a background thread and show the chromatogram as soon as it's ready.
- Clicking a scan: spectrum update feels instant (<100 ms).
- Dragging a range of ~100 MS1 scans: result in well under a second.
- Downsample the chromatogram for display (pyqtgraph `setDownsampling` / `setClipToView`) so zooming stays smooth.

## Suggested project structure

```
msview/
  pyproject.toml
  README.md
  msview/
    __main__.py        # CLI entry: msview [file]
    app.py             # QApplication, main window, layout
    panes.py           # ChromatogramPane, SpectrumPane, MS2Pane
    data/
      base.py          # RawSource interface
      thermo.py        # RawFileReader via pythonnet
      mzml.py          # optional fallback backend
    averaging.py
    index.py           # scan index construction
  tests/
```

## Milestones

1. **Data layer**: load a `.raw` file via pythonnet, build the scan index, fetch single spectra, fetch TIC. Test from a script with a sample file.
2. **Static three-pane window**: chromatogram + MS1 + MS2 plotted for a hard-coded scan.
3. **Chromatogram click** → MS1 single scan + MS2 list for that cycle, with stepping.
4. **Drag to average** MS1 + MS2 list for the range.
5. Precursor markers on MS1, peak labels, hover readout, keyboard navigation.
6. Polish: drag-and-drop open, export CSV/image, TIC/BPC toggle, threading, caching.
7. Optional: package as a double-clickable `.app` (PyInstaller or briefcase).

## Acceptance criteria

- Opens a Thermo Orbitrap DDA `.raw` file on an Apple Silicon Mac without converting it.
- Click on the chromatogram updates MS1 and resets the MS2 list to that duty cycle.
- Drag on the chromatogram averages MS1 across the range and lists all MS2s in the range.
- MS2 stepping works with buttons and arrow keys, and the header shows scan #, RT, precursor m/z, charge, activation, and energy.
- The current MS2's precursor and isolation window are marked on the MS1.
- Interaction stays responsive on a full-length run.

## Out of scope (for now)

Peak picking, identification/search, XIC extraction by m/z (good future feature), multi-file overlay, DIA-specific views, non-Thermo vendor formats.
