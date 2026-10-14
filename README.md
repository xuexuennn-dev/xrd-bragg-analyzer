# XRD Bragg Analyzer: Diffractometer Signal Processing & Lattice Spacing Toolkit

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Scientific Computing](https://img.shields.io/badge/Stack-NumPy%20%7C%20SciPy%20%7C%20Matplotlib-orange.svg)](https://scipy.org/)

A modular scientific computing toolkit tailored for materials science researchers, crystallographers, and automated LLM-assisted laboratory pipelines. It automates raw X-ray diffractometer data ingestion (including proprietary Rigaku binary formats and columnar ASCII datasets), applies Savitzky-Golay noise filtering, computes interplanar lattice spacings ($d$-spacings) using Bragg's Law, and outputs publication-ready crystallographic peak indexing profiles.

---

## 🔬 Core Scientific Principles

### 1. Bragg's Law Formulation
The fundamental relationship governing constructive interference from crystallographic lattice planes is defined by **Bragg's Law**:

$$n\lambda = 2d \sin\theta \quad \Longrightarrow \quad d = \frac{\lambda}{2 \sin\theta}$$

Where:
* $d$: Interplanar lattice spacing ($\text{Å}$)
* $\lambda$: Incident X-ray radiation wavelength (Default: $\text{Cu } K\alpha = 1.54060\text{ Å}$)
* $\theta$: Bragg diffraction angle ($\text{rad}$, derived from experimental $2\theta$ scan values)
* $n$: Order of reflection (typically $n = 1$)

### 2. Digital Filtering & Peak Detection
Raw diffractometer scans often exhibit high-frequency detector noise and instrument-specific background scattering. This toolkit implements:
* **Savitzky-Golay Digital Filter**: Local polynomial convolution that preserves higher-moment peak features (FWHM, skewness, amplitude) without introducing artificial peak-broadening artifacts.
* **Topological Prominence Finding**: Dynamically identifies genuine diffraction reflections against local background baselines using relative amplitude thresholds and minimum angular distances.

---

## 🚀 Key Features

* **Dual-Mode Diffractometer Parser**:
  * Native byte-level extraction for Rigaku binary diffractometer files (`.raw` header offsets and float32 arrays).
  * Robust tabular parser for standard 2-column delimited text/ASCII exports (`.xy`, `.txt`, `.csv`, `.dat`).
* **Automated Peak Indexing**:
  * Precision angular location of $2\theta$ reflection centroids.
  * Instantaneous conversion of observed $2\theta$ values into crystallographic $d$-spacing ($\text{Å}$).
  * Calculation of normalized relative intensities ($I/I_{\max} \times 100\%$).
* **Publication-Quality Visualization**:
  * Generates vector-grade figures featuring raw diffraction profiles, smoothed envelopes, and automatically annotated peak indices ($2\theta$ and $d$-values).
* **CLI & Python API**: Usable both as an importable module in scientific scripts or as a standalone terminal utility.

---

## 📊 Sample Visualization

Below is an automated analysis generated from an experimental diffraction scan:

![XRD Analysis Example](xrd_example.png)

---

## 🛠️ Installation

```bash
git clone https://github.com/xuexuennn-dev/xrd-bragg-analyzer.git
cd xrd-bragg-analyzer
pip install -r requirements.txt
```

---

## 💻 Usage

### Command Line Interface (CLI)

Run an automated peak search and generate an annotated plot:

```bash
python xrd_analyzer.py sample_xrd_data.txt --wavelength 1.5406 --prominence 50 --output result.png
```

#### CLI Arguments:
* `input_file`: Path to raw diffractometer data file (`.txt`, `.xy`, `.csv`, or `.raw`).
* `--wavelength`, `-w`: X-ray wavelength in $\text{Å}$ (default: `1.5406` for $\text{Cu } K\alpha$).
* `--prominence`, `-p`: Minimum peak prominence for peak detection thresholding (default: `50.0`).
* `--distance`, `-d`: Minimum angular separation (in data points) between adjacent peaks (default: `15`).
* `--no-smooth`: Disable Savitzky-Golay filtering.
* `--output`, `-o`: Filepath to save the resulting figure (e.g. `analysis.png`).

---

### Python API

Integrate directly into automated data pipelines:

```python
from xrd_analyzer import load_xrd_data, analyze_xrd, plot_xrd

# 1. Load data (ASCII or Rigaku binary)
two_theta, intensity = load_xrd_data("sample_xrd_data.txt")

# 2. Execute peak indexing and Bragg calculations
results = analyze_xrd(
    two_theta, 
    intensity, 
    wavelength=1.5406, 
    prominence=50.0
)

# 3. Inspect detected crystallographic reflections
for peak in results['peaks']:
    print(f"2θ: {peak['2theta']:.3f}° | d-spacing: {peak['d_spacing']:.4f} Å | Relative Intensity: {peak['rel_intensity']:.1f}%")

# 4. Generate annotated visualization
plot_xrd(results, save_path="xrd_indexed.png")
```

---

## 📁 Repository Structure

```
xrd-bragg-analyzer/
├── sample_xrd_data.txt    # Representative experimental diffractometer scan
├── requirements.txt       # Project runtime dependencies (NumPy, SciPy, Matplotlib)
├── xrd_analyzer.py        # Core processing engine, binary parser & CLI
├── xrd_example.png        # Sample rendered crystallographic analysis
└── README.md              # Technical documentation & usage guide
```

---

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 👤 Author

**Ye Xu**
* GitHub: [@xuexuennn-dev](https://github.com/xuexuennn-dev)
* LinkedIn: [Ye Xu](https://www.linkedin.com/in/ye-xu-guet/)