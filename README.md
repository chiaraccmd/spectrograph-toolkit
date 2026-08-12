# 🔧 Spectrograph Design Toolkit

*A toolkit for astronomical instrumentation development, including spectrograph design, fibre optics analysis, optical performance modelling, Zemax optimization, and experimental spectral data processing.*

---

## 🚀 Overview

This repository collects MATLAB and Zemax tools developed for the design and analysis of astronomical spectrographs, with particular emphasis on fibre-fed and integral field spectroscopy.

The workflow covers several stages of instrument development:

* first-order spectrograph sizing;
* geometrical and diffraction-limited performance analysis;
* fibre-image overlap and crosstalk modelling;
* Zemax optical optimization;
* detector data extraction and wavelength calibration;
* experimental resolving-power measurements.

The tools were originally developed during MSc thesis work on an integral field spectrograph for exoplanet science.

---

## 🚀 Quick Start

### MATLAB Examples

```matlab
% Multi-band spectrograph parameter analysis
[optimal_params, analysis_data] = spectrograph_parameter_sweep({'Y','J','H'}, ...
    'resolving_power', [5000,5000,5000], ...
    'name', 'MCIFU_5000_950');

% Geometrical spectrograph analysis
[performance_metrics, geometric_params] = spectrograph_geometric_analysis(...
    'R_Y', 7880, ...
    's1', 7.3e-6, ...
    'nPix', 2000, ...
    'pix', 18e-6);

% Fibre crosstalk analysis
[crosstalk_results, analysis_data] = fibre_crosstalk_simulator('airy', ...
    'fibre_separation', 25e-6, ...
    'wavelength', 1.55e-6);

% Geometric-to-diffraction-limited transition
[transition_data, performance_metrics] = diffraction_limit_analysis(...
    'grating_density', 650e3, ...
    'beam_size', 14.8e-3, ...
    'f_number', 3.57);
```

---

## 📁 Repository Structure

```text
spectrograph-toolkit/
├── 📊 MATLAB/
│   ├── Optical_Geometry/
│   │   ├── spectrograph_parameter_sweep.m
│   │   ├── spectrograph_geometric_analysis.m
│   │   ├── diffraction_limit_analysis.m
│   │   ├── geometric_spectrograph_sweep.m
│   │   └── geometric_spectrograph_evaluation.m
│   │
│   ├── Fibre_Optics/
│   │   └── fibre_crosstalk_simulator.m
│   │
│   └── Data_Processing/
│       ├── extract_1d_spectrum.m
│       ├── wavelength_calibration.m
│       ├── apply_wavelength_calibration.m
│       ├── measure_spectral_resolution.m
│       └── integrate_spectral_peak.m
│
├── 🔍 Zemax_Templates/
│   ├── Merit_Functions/
│   │   ├── collimator_optimization.MF
│   │   └── spectrograph_optimization.MF
│   │
│   └── Macros/
│       └── glass_substitution_tool.zpl
│
└── 🧪 Examples/
    ├── airy_psf_example.png
    ├── spectrograph_transition_example.png
    └── calibrated_spectrum_example.png
```

---

## 🧰 Tool Categories

### ✅ Optical System Analysis

* **spectrograph_parameter_sweep.m** — Multi-band spectrograph parameter analysis for Y, J, and H bands, including grating-density matching between spectral channels.
* **spectrograph_geometric_analysis.m** — First-order analysis of resolving power, detector coverage, diffraction sampling, fibre-image separation, and optical geometry.
* **geometric_spectrograph_sweep.m** — Explores slit-limited spectrograph configurations over resolving power, slit width, grating density, detector geometry, and optical parameters.
* **geometric_spectrograph_evaluation.m** — Evaluates dispersion, sampling, detector coverage, and resolving power for a selected geometrical spectrograph configuration.
* **diffraction_limit_analysis.m** — Compares slit-limited and diffraction-limited resolving power and determines the transition between the two regimes.

---

### ✅ Fibre Optics & IFS

* **fibre_crosstalk_simulator.m** — Models leakage between neighbouring fibre images using Airy, Gaussian, and dispersed-spectrum PSFs.
* Includes both continuous and detector pixel-integrated calculations.
* Supports analysis of fibre separation and PSF evolution with wavelength.

---

### ✅ Experimental Data Processing

* **extract_1d_spectrum.m** — Extracts a one-dimensional spectrum from a raw detector image using spatial-region detection and local background subtraction.
* **wavelength_calibration.m** — Derives a polynomial wavelength solution from Neon reference lines and evaluates calibration residuals.
* **apply_wavelength_calibration.m** — Applies a previously determined wavelength calibration to an extracted spectrum.
* **measure_spectral_resolution.m** — Measures the FWHM of an isolated spectral line and estimates the experimental resolving power.
* **integrate_spectral_peak.m** — Performs local baseline subtraction and numerical integration of spectral peaks.

---

### ✅ Zemax Integration

* **collimator_optimization.MF** — Merit function for collimator optimization using focal-length, throughput, and angular-aberration constraints.
* **spectrograph_optimization.MF** — Merit function for full spectrograph optimization using first-order constraints and RMS spot performance.
* **glass_substitution_tool.zpl** — ZPL macro for testing alternative optical materials during optimization.

---

## 🔬 Design & Analysis Workflow

The repository is organized around a typical spectrograph development sequence:

```text
Instrument requirements
        ↓
First-order parameter sweep
        ↓
Geometrical spectrograph sizing
        ↓
Zemax optical optimization
        ↓
Diffraction and fibre crosstalk analysis
        ↓
Experimental detector acquisition
        ↓
1D spectral extraction
        ↓
Wavelength calibration
        ↓
Measured spectral performance
```

The MATLAB models are intended primarily for first-order design studies and parameter exploration.

Detailed imaging performance should subsequently be evaluated using a complete optical model and, where available, experimental measurements.

---

## 🔭 Spectrograph Models

The analytical calculations use standard relations from geometrical and Fourier optics.

### Grating Geometry

For first-order Littrow operation:

$$
m\lambda = 2d\sin\alpha
$$

where $d$ is the groove spacing and $\alpha$ is the Littrow angle.

### Geometrical Resolving Power

The slit-limited resolving power is evaluated from the projected slit width and spectrograph geometry.

The general model used in the toolkit is:

$$
R_{\mathrm{geom}} = \frac{mG\lambda F W}{s}
$$

where:

- $m$ is the diffraction order;
- $G$ is the grating line density;
- $F$ is the relevant focal ratio;
- $W$ is the illuminated grating width;
- $s$ is the slit width.

### Diffraction-Limited Resolving Power

The diffraction limit is estimated using the Airy-disk criterion:

$$
R_{\mathrm{diff}} = \frac{mGW}{1.22}
$$

The transition between slit-limited and diffraction-limited behaviour occurs when:

$$
R_{\mathrm{geom}} = R_{\mathrm{diff}}
$$

---

## 🔍 Detector Sampling & Spectral Coverage

Detector sampling is evaluated from the wavelength dispersion and projected slit or PSF size.

The local wavelength dispersion is derived from the grating equation and camera focal length.

The tools evaluate:

- wavelength coverage on the detector;
- number of detector pixels required;
- pixels per resolution element;
- resolving power as a function of wavelength;
- spectral extent across the focal plane.

These calculations can be used to compare detector formats and grating configurations during preliminary design.

---

## 🌈 Fibre Crosstalk Modelling

The fibre-analysis routine estimates contamination between neighbouring fibre spectra using several PSF models.

Supported models include:

- Airy PSF;
- Gaussian PSF;
- dispersed spectral PSF.

For each model, crosstalk is estimated from the fraction of the source PSF falling inside a neighbouring extraction region.

Detector effects can also be included through pixel integration.

The analysis is intended as a first-order estimate of fibre-image leakage and extraction overlap rather than a complete detector or extraction-pipeline simulation.

---

## 🧪 Experimental Spectral Analysis

The experimental routines connect the analytical design to measurements from a real detector.

A typical reduction sequence is:

```text
raw detector image
        ↓
extract_1d_spectrum.m
        ↓
wavelength_calibration.m
        ↓
apply_wavelength_calibration.m
        ↓
measure_spectral_resolution.m
        ↓
integrate_spectral_peak.m
```

### Spectrum extraction

The detector image is reduced by:

1. identifying the illuminated spatial region;
2. selecting a nearby dark region;
3. estimating the column-dependent background;
4. subtracting the background;
5. integrating the signal along the spatial direction.

### Wavelength Calibration

Detected Neon emission lines are matched to known reference wavelengths.

A polynomial relation

$$
\lambda = \lambda(x)
$$

is fitted between detector pixel coordinate $x$ and wavelength.

Calibration quality is evaluated through:

- RMS residual in nanometres;
- RMS residual in detector pixels;
- maximum wavelength residual;
- local spectral dispersion.

### Experimental Resolving Power

For an isolated spectral line, the measured resolving power is estimated from:

$$
R = \frac{\lambda_0}{\Delta\lambda_{\mathrm{FWHM}}}
$$

where $\lambda_0$ is the measured line centre and $\Delta\lambda_{\mathrm{FWHM}}$ is its full width at half maximum.

This allows direct comparison between analytical predictions and experimental performance.

---

## 🔍 Zemax Optimization

The Zemax merit functions complement the MATLAB first-order calculations.

### Collimator optimization

The collimator merit function includes:

* effective focal-length constraints;
* throughput constraints;
* angular aberration minimization.

Angular aberration operands are used because the main objective of the collimator is to produce a well-collimated output beam.

### Spectrograph optimization

The full spectrograph merit function includes:

* subsystem focal-length constraints;
* throughput requirements;
* RMS transverse aberration optimization at the detector plane.

This allows the first-order geometry obtained in MATLAB to be refined using full ray tracing.

---

## 📊 Example Applications

The toolkit can be used for:

* astronomical spectrograph sizing;
* integral field spectrograph design;
* fibre-fed spectroscopy;
* detector selection and sampling analysis;
* diffraction-performance budgeting;
* fibre-crosstalk estimation;
* Zemax optimization;
* laboratory spectrograph characterization;
* wavelength calibration and spectral-resolution measurements.

---

## ⚠️ Scope & Assumptions

The analytical tools use simplified first-order optical models.

They are intended for:

* conceptual design;
* parameter exploration;
* comparison of candidate configurations;
* preliminary detector and grating selection.

They are **not substitutes for complete optical modelling**.

Final instrument performance should be validated through full optical simulation, tolerance analysis, and experimental measurements.

---

## 🔄 Planned Extensions

* **VPH Grating Design** — Efficiency calculations and Bragg-condition optimization.
* **Tolerance Analysis** — Sensitivity to manufacturing and alignment errors.
* **Additional Zemax Tools** — Multi-configuration and automated optimization utilities.
* 
---

## 📝 License

This toolkit is available under the **MIT License** for academic and research use.

See the `LICENSE` file for details.

---

## 📚 Background

The tools in this repository were developed during MSc thesis work on:

**“Development of an Integral Field Spectrograph for Exoplanet Science”** at **Politecnico di Milano** and **INAF – Osservatorio Astronomico di Brera**.
