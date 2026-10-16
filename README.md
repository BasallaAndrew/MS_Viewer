# msview

A fast, minimal desktop viewer for Thermo `.raw` LC-MS files on macOS, in the spirit of a stripped-down Xcalibur Qual Browser. Open a file, click the chromatogram, and the spectra update immediately. There is no conversion step. mzML files work too.

Built with PySide6 (Qt) and pyqtgraph. It reads `.raw` files directly with Thermo's RawFileReader through pythonnet. See [ms-viewer-spec.md](ms-viewer-spec.md) for the design.

## Setup

Requires Python 3.11+ (tested on 3.12, Apple Silicon).

```bash
python3.12 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

That is enough for mzML files. Thermo `.raw` files need two more things.

### 1. .NET 8 runtime

A user-local install that doesn't need sudo:

```bash
curl -sSL https://dot.net/v1/dotnet-install.sh -o /tmp/dotnet-install.sh
bash /tmp/dotnet-install.sh --channel 8.0 --runtime dotnet --install-dir ~/.dotnet
```

msview looks for .NET in `~/.dotnet`, `/usr/local/share/dotnet` and `/opt/homebrew/share/dotnet`. Set `DOTNET_ROOT` if yours is somewhere else.

### 2. Thermo RawFileReader DLLs

Thermo's license does not allow redistributing these, so each user downloads them from <https://github.com/thermofisherlsms/RawFileReader>. That also means accepting Thermo's license. Put them in `~/.msview/rawfilereader/`, or in a folder named by `MSVIEW_RAWFILEREADER_DIR`.

The .NET 8 build ships as NuGet packages (zip files) in `Libs/NetCore/Net8`. You also need its `OpenMcdf` dependency from nuget.org:

```bash
git clone --depth 1 https://github.com/thermofisherlsms/RawFileReader.git /tmp/rfr
mkdir -p ~/.msview/rawfilereader
cd /tmp/rfr/Libs/NetCore/Net8
for p in data rawfilereader; do
  unzip -o -j thermofisher.commoncore.$p.*.nupkg 'lib/net8.0/*.dll' -d ~/.msview/rawfilereader
done
curl -sSL https://www.nuget.org/api/v2/package/OpenMcdf/3.1.4 -o /tmp/openmcdf.nupkg
unzip -o -j /tmp/openmcdf.nupkg 'lib/net8.0/OpenMcdf.dll' -d ~/.msview/rawfilereader
```

The folder should then contain `ThermoFisher.CommonCore.Data.dll`, `ThermoFisher.CommonCore.RawFileReader.dll` and `OpenMcdf.dll`. If anything is missing, msview shows an error that says what.

## Usage

```bash
.venv/bin/msview                   # a file picker opens (it remembers the last folder)
.venv/bin/msview path/to/run.raw   # or .mzML; opens it directly
```

Run these from the project folder. While the viewer is open you can also open files with ⌘O, or by dragging them onto the window.

### Controls

| Where | Action | Effect |
|---|---|---|
| Chromatogram | click | Show the nearest MS1; the MS2 list becomes that duty cycle's MS2s |
| Chromatogram | drag | Average the MS1s in that RT range; the MS2 list becomes every MS2 in the range |
| Chromatogram | right- or middle-drag | Pan |
| Any pane | scroll / shift+scroll | Zoom x / zoom y |
| Any pane | double-click | Reset zoom |
| Any pane | right-click | Copy image, export CSV |
| Spectra | drag | Pan m/z (the y-axis follows the tallest visible peak) |
| MS1 | click a peak | Jump to the MS2 of that precursor (see below) |
| MS1 / MS2 | right-click a peak | Show the XIC for that m/z |
| Keyboard | ← / → | Previous / next MS1; after clicking the MS2 pane, previous / next MS2 |
| Keyboard | ↑ / ↓ (or ◀ ▶) | Previous / next MS2 in the list, from any pane |
| Toolbar | Lock isolation window | Keep the current MS2's isolation window fixed on the MS1 while you step through time |

The pane you last clicked has a highlighted header and decides what ← / → do. The chromatogram also shows a dashed orange line at the RT of the MS2 being shown.

Finding the MS2 for an MS1 peak: click a peak in the MS1 pane. If an MS2 in the current list isolated it, msview jumps there. Otherwise it searches the run within ±1 min of the MS1 being shown, or of the selected range. A precursor match within 20 ppm, or a match to one of its 13C isotopes, ranks above a peak that merely fell inside the isolation window. The MS2 list then becomes every MS2 of that precursor, and the one closest in RT is shown, so ↑ / ↓ step through its repeat fragmentations. Clicking the chromatogram goes back to the normal duty-cycle list.

Lock isolation window: with an MS2 showing, tick **Lock isolation window** in the toolbar. That MS2's precursor line and isolation band on the MS1 stay where they are, and turn purple, while you step left/right, click the chromatogram or average a range. This makes it easy to follow one precursor's isotope envelope over time. While locked, each new cycle's MS2 list opens on the MS2 that fragmented the locked m/z, if there is one. Untick it to go back to normal.

Toolbar: switch the chromatogram between TIC, base peak and XIC. The status bar shows the values under the cursor.

XIC: choose **XIC** in the dropdown, type an m/z and press Return (or click Extract). Set the tolerance in ppm (default 10) or Da (default 0.01). The m/z box starts out filled with the current MS2's precursor. The XIC sums the MS1 intensity inside the window for every scan. It takes about 1 s for a Thermo file and a few seconds for a large mzML. Results are cached, so switching back to an XIC you've already extracted is instant.

### Method view

Set **View** in the toolbar to **Method & run info** to replace the MS1/MS2 panes with:

- **Gradient plot.** The programmed %B is read from the LC pump method and drawn dashed. The %B the pump actually logged is drawn solid. **Right axis** plots any other logged trace: pressure by default, or flow, column temperature, spray current, capillary temperature and so on. The chromatogram selection is mirrored here, and hovering shows every value at that time.
- **Parameters tab.** A searchable table with an acquisition summary computed from the scans (cycle time, MS2 per cycle, activation and energy, isolation width, charges), plus sample, instrument, tune and each LC/MS device's method settings. Type in the filter box to narrow it down (for example "resolution" or "AGC"). ⌘C copies the selected rows; right-click to copy everything.
- **Method text tab.** The full instrument method exactly as stored in the file.

How much there is depends on the file. Thermo `.raw` files carry the method, tune data and device status logs. The gradient table parser understands Agilent pump methods and Dionex/Vanquish-style `%B` scripts. mzML files keep only the instrument, software and source file, plus any pump-pressure traces the converter kept.

Notes:
- Averaging uses at most 500 scans. A wider range is sampled evenly across the range, with a warning in the status bar. Centroid scans are merged into 5 ppm bins. Profile scans are interpolated onto a shared m/z grid.
- For Thermo files, spectra come from the centroid stream when one exists, as the spec requires. Otherwise the profile data is used.

### UniDec view (intact mass)

Set **View** to **UniDec (intact mass)** to deconvolve the MS1 you've selected with [UniDec](https://github.com/michaelmarty/UniDec):

- **Top pane:** the input spectrum. This is a single scan when you click the chromatogram, or an average when you drag a range across the protein peak, which usually works better. UniDec's fit is overlaid as a red dashed line, and coloured ▼ markers with charge labels show the charge states of the masses it found.
- **Bottom pane:** the deconvolved zero-charge mass spectrum, labelled with the picked masses. The header lists the top five with their relative heights.
- **Controls:** charge range, mass range, mass bin, peak FWHM (auto-estimated unless you untick "auto"), peak window (masses closer than this merge; the 20 Da default keeps glycoforms apart) and threshold. **Only zoomed m/z range** deconvolves just what's visible in the top pane. With **Auto-run** on, it reruns whenever the selection or a setting changes, but only while the view is showing. Otherwise click **Deconvolve**.

A run takes well under a second for a typical protein. The very first run after installing takes longer (up to ~30 s) while UniDec compiles itself. Right-click either pane to copy the image or export a CSV of the masses. UniDec is installed from PyPI without its GUI (`unidec` in `requirements.txt`). The macOS build includes its native engine, so nothing else is needed.

## Project layout

```
msview/
  __main__.py      CLI entry point
  app.py           main window, selection state, background work
  panes.py         chromatogram / MS1 / MS2 panes and mouse handling
  method_pane.py   Method view: gradient plot, parameter table, method text
  method.py        method/run info model and method-text / gradient parsers
  unidec_pane.py   UniDec view: input spectrum + fit, deconvolved masses, settings
  deconv.py        UniDec wrapper (runs in a temp folder, returns arrays)
  index.py         scan-index dtype and selection logic (pure functions)
  averaging.py     centroid and profile averaging
  data/
    base.py        RawSource interface, Spectrum, LRU cache
    thermo.py      Thermo .raw via RawFileReader + pythonnet
    mzml.py        mzML reader (offset index, header-only parsing, on-demand decoding)
tests/
```

The GUI only talks to `RawSource`. To add a new format, implement `_build_index()` and `_read_spectrum()`.

## Tests

```bash
.venv/bin/pytest
```

The tests use small generated mzML files. When the sample files are present in `data/` (git ignores them), extra tests also run against the real `.raw` and `.mzML` files.
