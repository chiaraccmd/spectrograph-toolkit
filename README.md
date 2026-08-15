# 🔧 Spectrograph Design Toolkit

*A MATLAB toolkit for first-order spectrograph design, diffraction and fibre-optics analysis, and experimental spectral characterization.*

---

## Overview

This repository collects MATLAB tools for the design and analysis of spectrographs, with particular emphasis on astronomical instrumentation, fibre-fed systems, and integral field spectroscopy.

The toolkit covers several stages of instrument development:

- first-order spectrograph sizing and parameter sweeps;
- multi-band spectrograph design and cross-band camera matching;
- slit-limited and diffraction-limited performance analysis;
- detector sampling and spectral coverage;
- fibre-image overlap and crosstalk modelling;
- experimental spectrum extraction and wavelength calibration;
- measured spectral-resolution and peak analysis.

The analytical models are intended primarily for preliminary design, parameter exploration, and physical performance assessment before detailed optical modelling.

---

## 🚀 Quick Start

Clone the repository and add the MATLAB folders to your path:

```matlab
addpath(genpath('MATLAB'));
```

The examples below illustrate the main analysis workflows. Additional functions and options are described in the sections that follow.

### Multi-band spectrograph design

```matlab
[matched, data] = multiband_spectrograph_sweep({'Y','J','H'});
```

The routine evaluates spectrograph parameters over a range of grating densities and performs **cross-band camera matching**.

A reference configuration is selected for one band, and the remaining bands are matched to its camera focal length. The resulting collimator focal lengths, beam sizes, grating densities, and incidence angles can then be inspected to evaluate the compatibility of the different spectral channels.

### Slit-limited spectrograph sweep

```matlab
[T, best, results] = slit_limited_spectrograph_sweep(...
    'R_list', 1300:50:1500, ...
    'grating_densities', 1511e3);
```

The routine explores candidate slit-limited configurations and evaluates their optical geometry, detector sampling, spectral coverage, and feasibility.

### Diffraction-limit analysis

```matlab
[transition_data, performance_metrics] = diffraction_limit_analysis(...
    'grating_density', 650e3, ...
    'beam_size', 14.8e-3, ...
    'f_number', 3.57);
```

The routine compares slit-limited and diffraction-limited resolving power and identifies the transition between the two regimes.

### Fibre crosstalk analysis

```matlab
[crosstalk_results, analysis_data] = fibre_crosstalk_simulator('airy', ...
    'fibre_separation', 25e-6, ...
    'wavelength', 1.55e-6);
```

The simulator provides Airy, Gaussian, and wavelength-dependent dispersed PSF models for estimating contamination between neighbouring fibre images.

---

## 📁 Repository Structure

```text
spectrograph-toolkit/
├── MATLAB/
│   ├── Optical_Geometry/
│   │   ├── multiband_spectrograph_sweep.m
│   │   ├── slit_limited_spectrograph_sweep.m
│   │   ├── diffraction_limited_multiband_analysis.m
│   │   └── diffraction_limit_analysis.m
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
├── Examples/
│   ├── multiband_spectrograph_design.png
│   ├── diffraction_limit_transition.png
│   ├── fibre_crosstalk_airy.png
│   └── experimental_spectrum_calibration.png
│
├── .gitignore
├── LICENSE
└── README.md
```

---

## 🔭 Spectrograph Design

### Multi-band parameter analysis

`multiband_spectrograph_sweep.m` explores first-order spectrograph parameters for multiple wavelength bands.

For each channel, the routine evaluates quantities such as:

- grating density;
- Littrow angle;
- collimator focal length;
- camera focal length;
- beam diameter;
- illuminated grating width.

The cross-band matching procedure uses one band as a reference and searches for configurations in the remaining bands that reproduce its camera focal length.

This provides a first-order way to investigate whether multiple spectral channels can share similar camera requirements.

The remaining optical parameters are resulting design quantities rather than additional matching constraints and should therefore be inspected to assess the overall compatibility of the channels.

### Slit-limited parameter sweep

`slit_limited_spectrograph_sweep.m` explores spectrograph configurations over user-defined ranges of:

- resolving power;
- slit width;
- grating line density;
- detector format;
- Littrow or off-Littrow geometry.

For every valid configuration, the routine calculates:

- incidence and diffraction angles;
- camera and collimator focal lengths;
- beam dimensions;
- slit-image magnification;
- detector sampling;
- spectral coverage;
- slit-limited resolving power.

Feasibility checks include:

- full spectral-band diffraction;
- detector coverage;
- camera focal-length limits;
- minimum camera focal ratio.

Valid configurations are returned as a sortable table, with feasible configurations ranked using a simple figure of merit.

### Diffraction-limited multi-band analysis

`diffraction_limited_multiband_analysis.m` evaluates resolving power across multiple spectral bands while distinguishing between geometrical and diffraction-limited performance.

For each wavelength, the physically relevant resolving power is determined by the limiting mechanism:

```math
R_{\mathrm{eff}}(\lambda)
=
\min\left[
R_{\mathrm{geom}}(\lambda),
R_{\mathrm{diff}}
\right]
```

The routine therefore shows both the theoretical slit-limited resolving power and the diffraction ceiling, while identifying which one actually determines instrument performance.

### Diffraction-limit transition analysis

`diffraction_limit_analysis.m` examines the transition between the slit-limited and diffraction-limited regimes.

It is intended for sensitivity studies and for identifying the point at which increasing geometrical resolving power ceases to improve the effective resolution of the instrument.

---

## 📐 Spectrograph Models

The analytical calculations use standard first-order relations from geometrical and Fourier optics.

### Grating geometry

For diffraction order $m$,

```math
m\lambda = d(\sin\alpha + \sin\beta)
```

where $d$ is the groove spacing, $\alpha$ is the incidence angle, and $\beta$ is the diffraction angle.

For Littrow operation,

```math
\alpha = \beta
```

and therefore

```math
m\lambda = 2d\sin\alpha
```

Using the grating line density $G = 1/d$,

```math
\sin\alpha = \frac{mG\lambda}{2}
```

### Angular dispersion

Differentiating the grating equation gives

```math
\frac{d\beta}{d\lambda}
=
\frac{mG}{\cos\beta}
```

For a camera of focal length $f_2$, the local linear wavelength dispersion is

```math
\frac{d\lambda}{dx}
=
\frac{\cos\beta}{mGf_2}
```

This relation is used for detector sampling and first-order spectral-coverage calculations.

### Slit-limited resolving power

The geometrical resolving power used in the toolkit is

```math
R_{\mathrm{geom}}
=
\frac{mG\lambda F W}{s}
```

where:

- $m$ is the diffraction order;
- $G$ is the grating line density;
- $\lambda$ is the wavelength;
- $F$ is the relevant focal ratio;
- $W$ is the illuminated grating width;
- $s$ is the entrance-slit width.

For fixed spectrograph geometry,

```math
R_{\mathrm{geom}} \propto \lambda
```

### Diffraction-limited resolving power

Using the Airy-disk criterion adopted in the toolkit,

```math
R_{\mathrm{diff}}
=
\frac{mGW}{1.22}
```

For fixed grating illumination, the diffraction-limited resolving power is independent of wavelength.

The transition between the two regimes occurs when

```math
R_{\mathrm{geom}}
=
R_{\mathrm{diff}}
```

The effective resolving power is therefore

```math
R_{\mathrm{eff}}
=
\min\left(
R_{\mathrm{geom}},
R_{\mathrm{diff}}
\right)
```

This distinction is important when interpreting a calculated geometrical resolving power: once $R_{\mathrm{geom}}$ exceeds $R_{\mathrm{diff}}$, diffraction rather than the slit determines the effective spectral resolution.

---

## 🔍 Detector Sampling & Spectral Coverage

Detector sampling is evaluated from the local wavelength dispersion and the projected slit or PSF size.

For a spectral resolution element,

```math
\Delta\lambda
=
\frac{\lambda}{R}
```

The corresponding number of detector pixels per resolution element is

```math
N_{\mathrm{pix,res}}
=
\frac{\Delta\lambda}
{\left(\frac{d\lambda}{dx}\right)p}
```

where $p$ is the detector pixel pitch.

The toolkit evaluates:

- wavelength coverage on the detector;
- detector pixels required by the spectrum;
- pixels per resolution element;
- spectral extent across the focal plane;
- resolving power as a function of wavelength.

For broad wavelength ranges, the focal-plane extent can also be evaluated from the wavelength-dependent diffraction angles rather than assuming constant dispersion across the full detector.

---

## 🌈 Fibre Crosstalk Modelling

`fibre_crosstalk_simulator.m` estimates contamination between neighbouring fibre images at the detector.

Supported models include:

- Airy PSF;
- Gaussian PSF;
- wavelength-dependent dispersed PSF.

For each model, crosstalk is estimated from the fraction of the source PSF falling within a neighbouring extraction region.

Detector pixel integration can also be included, allowing the continuous optical PSF and detector-sampled response to be compared.

The Airy analysis additionally visualizes the projected fibre-core region relative to the diffraction pattern.

The dispersed model extends the analysis across wavelength, allowing wavelength-dependent PSF displacement and overlap to be investigated.

The simulator is intended as a first-order estimate of fibre-image leakage and extraction overlap rather than a complete detector or spectral-extraction simulation.

---

## 🧪 Experimental Spectral Analysis

The `Data_Processing` folder contains reusable functions for experimental spectrograph characterization.

A typical processing sequence is:

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

The individual functions can also be used independently.

### Spectrum extraction

`extract_1d_spectrum.m` converts a two-dimensional detector image into a one-dimensional spectrum.

The routine can:

1. identify the illuminated spatial region;
2. select a nearby background region;
3. estimate the column-dependent detector background;
4. subtract the background;
5. integrate the signal along the spatial direction.

Both slit-like and compact input geometries can be analysed.

### Wavelength calibration

`wavelength_calibration.m` derives a polynomial mapping between detector position and wavelength:

```math
\lambda = \lambda(x)
```

Reference emission lines are matched to measured detector peaks and fitted with a polynomial wavelength solution.

Calibration diagnostics include:

- wavelength residuals;
- RMS calibration residual;
- maximum residual;
- local spectral dispersion.

### Calibration application

`apply_wavelength_calibration.m` applies a previously derived calibration polynomial to an extracted spectrum.

This separates calibration generation from calibration application, allowing a wavelength solution to be reused for compatible measurements.

### Experimental resolving power

`measure_spectral_resolution.m` estimates the width of an isolated spectral line and its corresponding experimental resolving power.

For a measured line centred at $\lambda_0$,

```math
R
=
\frac{\lambda_0}
{\Delta\lambda_{\mathrm{FWHM}}}
```

where $\Delta\lambda_{\mathrm{FWHM}}$ is the measured full width at half maximum.

This provides an experimental performance metric that can be compared directly with analytical resolving-power predictions.

### Peak integration

`integrate_spectral_peak.m` performs local baseline estimation and numerical integration of an isolated spectral feature.

After baseline subtraction, the integrated signal is evaluated as

```math
A
=
\int_{x_1}^{x_2}
\left[
I(x)-I_{\mathrm{base}}(x)
\right]\,dx
```

with the numerical integral evaluated using the trapezoidal rule.

---

## 📊 Examples

The `Examples/` folder contains representative outputs from the different analysis stages:

- `multiband_spectrograph_design.png` — cross-band camera-matching example;
- `diffraction_limit_transition.png` — transition between slit-limited and diffraction-limited resolving power;
- `fibre_crosstalk_airy.png` — Airy PSF and fibre-core geometry;
- `experimental_spectrum_calibration.png` — example wavelength-calibration workflow.

The example outputs use generalized or synthetic parameters and are intended to demonstrate the analysis routines rather than reproduce a specific instrument configuration.

---

## 🔄 Planned Extensions

Possible future extensions include:

- VPH grating efficiency and Bragg-condition modelling;
- optical tolerance and sensitivity analysis;
- Zemax macros and functions.

---

## 📝 License

This project is licensed under the **MIT License** for academic and research use.

See [`LICENSE`](LICENSE) for details.

---

## Background

The toolkit originated from MSc thesis work on:

**“Development of an Integral Field Spectrograph for Exoplanet Science”**

at **Politecnico di Milano** and **INAF – Osservatorio Astronomico di Brera**.

The repository has subsequently been reorganized into reusable functions for spectrograph design, fibre-optics analysis, and experimental spectral characterization.
