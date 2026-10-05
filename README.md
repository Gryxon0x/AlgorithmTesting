# Step Detection Test Harness - Peak Detector Reference

This project contains a plain-C recreation of the peak-detection path from the supplied Python/SciPy implementation.

## What the detector reproduces

The detector performs:

1. Convert the CSV `ax_mg`, `ay_mg`, `az_mg` values to g.
2. Compute acceleration magnitude: `sqrt(ax^2 + ay^2 + az^2)`.
3. Design a 2nd-order Butterworth band-pass using the recording sample rate.
4. Apply a zero-phase forward/backward SOS filter equivalent to `scipy.signal.sosfiltfilt` for this filter.
5. Find local maxima.
6. Enforce the same minimum-distance ordering used by `scipy.signal.find_peaks` (larger peaks have priority).
7. Compute peak prominence and apply the configured minimum prominence.

Default configuration matches the Python `DetectionConfig`:

- low cutoff: 0.9 Hz
- high cutoff: 8.0 Hz
- minimum event interval: 0.28 s
- minimum prominence: 0.02 g
- no minimum height
- Butterworth order: 2

For `data/example_session.csv`, the Python reference produces 11 peaks at:

`14, 38, 67, 93, 120, 146, 172, 196, 214, 229, 249`

The C executable is expected to produce the same list.

## Important MCU note

This version is intentionally an **offline reference implementation**. `sosfiltfilt` filters forward and backward, so it needs future samples/the complete recording. Peak prominence can also inspect samples on both sides of a peak.

That makes this version ideal for validating the C math against the Python implementation, but it is not yet the final real-time CC1352R1/nRF52840 algorithm. A later streaming version should use causal filtering plus delayed peak confirmation, and we can compare its detections against this reference.

## Build

From the project root in the VS Code terminal:

```sh
cmake -S . -B build -G Ninja
cmake --build build
```

## Run

Windows PowerShell:

```powershell
.\build\step_test.exe data\example_session.csv
```

Linux/macOS:

```sh
./build/step_test data/example_session.csv
```

## Optional debug CSV

The second argument writes every filtered sample and marks detected peaks:

Windows:

```powershell
.\build\step_test.exe data\example_session.csv output\peak_debug.csv
```

Linux/macOS:

```sh
./build/step_test data/example_session.csv output/peak_debug.csv
```

The debug file contains:

- sample index
- local `t_ms`
- `phone_time_est_ms`
- acceleration magnitude in g
- filtered signal in g
- peak marker
- prominence in g

## Source layout

- `algorithms/peak_detector.c/.h` - peak detection algorithm
- `host/csv_reader.c/.h` - host-only CSV parser
- `host/main.c` - host-side test runner
- `include/accel_sample.h` - shared sample representation
- `data/example_session.csv` - supplied BMA400 recording
