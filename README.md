# 🔧 Spectrograph Design Toolkit

*A toolkit for astronomical instrumentation development, featuring spectrograph design, fibre optics analysis, and optical performance modelling. Developed for integral field spectrograph design and optimization.*

---

## 🚀 Quick Start

### MATLAB Examples

```matlab
% Spectrograph parameter optimization
[optimal_params, analysis_data] = spectrograph_parameter_sweep({'Y','J','H'}, ...
    'resolving_power', [5000,5000,5000], 'name', 'MCIFU_5000_950');

% Comprehensive geometric analysis
[performance_metrics, geometric_params] = spectrograph_geometric_analysis(...
    'R_Y', 7880, 's1', 7.3e-6, 'nPix', 2000, 'pix', 18e-6);

% Fibre crosstalk analysis
[crosstalk_results, analysis_data] = fibre_crosstalk_simulator('airy', ...
    'fibre_separation', 25e-6, 'wavelength', 1.55e-6);

% Diffraction limit analysis
[transition_data, performance_metrics] = diffraction_limit_analysis(...
    'grating_density', 650e3, 'beam_size', 14.8e-3, 'f_number', 3.57);
```

---

## 📁 Repository Structure

```
astrophotonics-toolkit/
├── 📊 MATLAB/
│   ├── Optical_Geometry/           # Spectrograph layout & analysis
│   │   ├── spectrograph_parameter_sweep.m
│   │   ├── spectrograph_geometric_analysis.m  
│   │   └── diffraction_limit_analysis.m
│   └── Fibre_Optics/               # Fibre bundle & crosstalk analysis
│       └── fibre_crosstalk_simulator.m
├── 🔍 Zemax_Templates/             # Optical design templates
│   ├── Merit_Functions/
│   │   └── collimator_optimization.MF
│   │   └── spectrograph_optimization.MF
│   └── Macros/
│       └── glass_substitution_tool.zpl
└── 🧪 Examples/
    ├── airy_psf_example.png
    └── spectrograph_transition_example.png
```

---

## 🧰 Tool Categories

### ✅ Optical System Analysis

* **spectrograph_parameter_sweep.m** — Multi-band parameter optimization with cross-band consistency analysis. Supports Y, J, H bands with automatic grating density matching.
* **spectrograph_geometric_analysis.m** — Comprehensive performance analysis including resolving power vs wavelength, detector coverage verification, and fibre crosstalk assessment.
* **diffraction_limit_analysis.m** — Identifies performance transition between geometric and diffraction-limited regimes with parameter sensitivity analysis.

### ✅ Fibre Optics & IFS

* **fibre_crosstalk_simulator.m** — Multi-model PSF analysis supporting Airy disk, Gaussian beam, and dispersed spectrum models with pixel integration for detector effects.

### ✅ Zemax Integration

* **collimator_optimization.MF** — Merit function template for collimator optimization in Zemax.
* **spectrograph_optimization.MF** — Merit function template for spectroscopic system optimization in Zemax.
* **glass_substitution_tool.zpl** — ZPL macro for automated material optimization in spectrograph designs.

---

## 🔄 Planned Extensions

* **VPH Grating Design** — Efficiency calculations and Bragg condition optimization.
* **Data Processing Utilities** — IFS datacube handling and spectral extraction.
* **Additional Zemax Templates** — Multi-configuration analysis and tolerance tools.

---

## 📝 License

This toolkit is available under the **MIT License** for academic and research use.

---

## 🤝 Contributing

We welcome contributions and enhancements! Please feel free to:

* Submit issues for bugs or feature requests
* Suggest additional tools or improvements
* Share your own spectrograph design utilities

*Tools developed during MSc thesis work on "Development of an Integral Field Spectrograph for Exoplanet Science" at Politecnico di Milano and INAF - Osservatorio Astronomico di Brera.*

