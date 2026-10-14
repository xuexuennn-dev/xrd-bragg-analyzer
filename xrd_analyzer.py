"""
XRD Bragg Analyzer
==================
A lightweight scientific computing toolkit for X-ray diffraction (XRD) data processing,
automated peak finding, Bragg's Law interplanar lattice spacing (d-spacing) calculation,
and crystallographic phase analysis.

Author: Ye Xu (github.com/xuexuennn-dev)
License: MIT
"""

import os
import sys
import struct
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter


def calculate_d_spacing(two_theta_deg, wavelength=1.5406):
    """
    Calculate interplanar spacing (d) using Bragg's Law:
        n * lambda = 2 * d * sin(theta)  =>  d = lambda / (2 * sin(theta))

    Parameters:
    -----------
    two_theta_deg : float or np.ndarray
        Diffraction angle (2-theta) in degrees.
    wavelength : float
        X-ray wavelength in Angstroms (default: 1.5406 A for Cu K-alpha).

    Returns:
    --------
    float or np.ndarray
        Interplanar lattice spacing (d) in Angstroms.
    """
    theta_rad = np.radians(two_theta_deg / 2.0)
    return wavelength / (2.0 * np.sin(theta_rad))


def read_rigaku_raw(filepath):
    """
    Parse proprietary Rigaku binary diffractometer files (.raw).
    Extracts start angle, stop angle, step size, and counts array.
    """
    with open(filepath, 'rb') as f:
        f.seek(0x0B92)
        start_angle = struct.unpack('<f', f.read(4))[0]
        f.seek(0x0B96)
        end_angle = struct.unpack('<f', f.read(4))[0]
        f.seek(0x0B9A)
        step = struct.unpack('<f', f.read(4))[0]
        num_points = int(round((end_angle - start_angle) / step)) + 1
        two_theta = np.arange(start_angle, end_angle + step / 2.0, step)[:num_points]
        f.seek(0x0C56)
        intensity = np.frombuffer(f.read(num_points * 4), dtype='<f')
    return two_theta, intensity


def parse_ascii_xrd(filepath):
    """
    Parse standard 2-column tabular/text XRD files (2-theta vs intensity).
    Supports comma-, space-, or tab-delimited formats.
    """
    two_theta = []
    intensity = []
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('//'):
                continue
            parts = line.replace(',', ' ').split()
            if len(parts) >= 2:
                try:
                    tt = float(parts[0])
                    inte = float(parts[1])
                    two_theta.append(tt)
                    intensity.append(inte)
                except ValueError:
                    continue
    if not two_theta:
        raise ValueError(f"No valid numeric XRD data found in {filepath}")
    return np.array(two_theta), np.array(intensity)


def load_xrd_data(filepath):
    """
    Load XRD data from either binary (.raw) or standard ASCII text format.
    """
    if filepath.lower().endswith('.raw'):
        try:
            return read_rigaku_raw(filepath)
        except Exception as e:
            print(f"[!] Warning: Falling back to ASCII parsing for {filepath}: {e}")
    return parse_ascii_xrd(filepath)


def analyze_xrd(two_theta, intensity, wavelength=1.5406, prominence=50.0, distance=15, smooth=True):
    """
    Execute automated XRD processing pipeline:
      1. Noise smoothing with Savitzky-Golay filter
      2. Peak detection via scipy.signal.find_peaks
      3. d-spacing calculation using Bragg's Law
      4. Relative intensity (I/I0) calculation
    """
    if smooth and len(intensity) > 15:
        # Window length 11, polynomial order 3
        smoothed_intensity = savgol_filter(intensity, window_length=11, polyorder=3)
        smoothed_intensity = np.maximum(smoothed_intensity, 0)
    else:
        smoothed_intensity = intensity

    # Peak detection
    peaks, properties = find_peaks(
        smoothed_intensity,
        prominence=prominence,
        distance=distance
    )

    peak_2theta = two_theta[peaks]
    peak_intensity = smoothed_intensity[peaks]
    d_spacings = calculate_d_spacing(peak_2theta, wavelength=wavelength)

    max_intensity = np.max(peak_intensity) if len(peak_intensity) > 0 else 1.0
    relative_intensity = (peak_intensity / max_intensity) * 100.0

    peak_records = []
    for tt, d, counts, rel_i in zip(peak_2theta, d_spacings, peak_intensity, relative_intensity):
        peak_records.append({
            '2theta_deg': float(tt),
            'd_spacing_A': float(d),
            'intensity': float(counts),
            'relative_intensity_pct': float(rel_i)
        })

    return {
        'smoothed_intensity': smoothed_intensity,
        'peaks': peaks,
        'peak_records': peak_records
    }


def plot_diffractogram(two_theta, raw_intensity, smoothed_intensity, peak_records, output_image="xrd_plot.png"):
    """
    Generate publication-ready XRD diffractogram with peak annotations.
    """
    plt.figure(figsize=(9, 5.5), dpi=300)
    plt.plot(two_theta, raw_intensity, color='#94a3b8', alpha=0.45, label='Raw Intensity')
    plt.plot(two_theta, smoothed_intensity, color='#0f172a', linewidth=1.2, label='Smoothed Profile')

    # Annotate prominent peaks
    for rec in peak_records:
        tt = rec['2theta_deg']
        d = rec['d_spacing_A']
        inte = rec['intensity']
        rel = rec['relative_intensity_pct']
        
        if rel >= 10.0:  # Only label peaks with >=10% relative intensity to prevent clutter
            plt.scatter(tt, inte, color='#dc2626', s=24, zorder=5)
            plt.annotate(
                f"{tt:.2f}°\n({d:.3f}Å)",
                xy=(tt, inte),
                xytext=(0, 7),
                textcoords="offset points",
                ha='center',
                fontsize=7.5,
                color='#1e293b',
                weight='bold'
            )

    plt.xlabel(r'Diffraction Angle $2\theta$ (degrees)', fontsize=11, fontweight='bold')
    plt.ylabel('Intensity (counts)', fontsize=11, fontweight='bold')
    plt.title(r'X-Ray Diffractogram & Peak Indexing ($\lambda = 1.5406\ \mathrm{\AA}$)', fontsize=12, fontweight='bold')
    plt.legend(frameon=True, facecolor='#ffffff', edgecolor='#cbd5e1')
    plt.grid(True, linestyle='--', alpha=0.35)
    plt.tight_layout()
    plt.savefig(output_image)
    plt.close()
    print(f"[+] Diffractogram exported: {output_image}")


def print_summary_table(peak_records):
    """
    Print formatted peak indexing summary.
    """
    print("\n" + "=" * 65)
    print(f"{'Peak #':<8}{'2-Theta (°)':<15}{'d-spacing (Å)':<18}{'Counts':<12}{'I/I0 (%)':<10}")
    print("-" * 65)
    for idx, r in enumerate(peak_records, 1):
        print(f"{idx:<8}{r['2theta_deg']:<15.3f}{r['d_spacing_A']:<18.4f}{r['intensity']:<12.1f}{r['relative_intensity_pct']:<10.1f}")
    print("=" * 65 + "\n")


def main():
    parser = argparse.ArgumentParser(description="XRD Bragg Analyzer - Automated Peak Indexing & d-spacing Calculator")
    parser.add_argument("input_file", help="Path to XRD data file (.raw binary or .txt/.csv ASCII)")
    parser.add_argument("--wavelength", type=float, default=1.5406, help="Radiation wavelength in Angstroms (default: 1.5406 A)")
    parser.add_argument("--prominence", type=float, default=50.0, help="Peak prominence threshold (default: 50.0)")
    parser.add_argument("--distance", type=int, default=15, help="Minimum peak distance in points (default: 15)")
    parser.add_argument("--out-plot", default="xrd_analysis_plot.png", help="Output plot filename")
    parser.add_argument("--out-txt", default=None, help="Optional text export for indexed peaks")

    args = parser.parse_args()

    if not os.path.exists(args.input_file):
        print(f"[-] Error: File not found: {args.input_file}")
        sys.exit(1)

    print(f"[*] Loading XRD dataset: {args.input_file}")
    two_theta, intensity = load_xrd_data(args.input_file)
    print(f"[+] Loaded {len(two_theta)} data points from 2-theta = {two_theta[0]:.2f}° to {two_theta[-1]:.2f}°")

    results = analyze_xrd(
        two_theta,
        intensity,
        wavelength=args.wavelength,
        prominence=args.prominence,
        distance=args.distance
    )

    print_summary_table(results['peak_records'])
    plot_diffractogram(
        two_theta,
        intensity,
        results['smoothed_intensity'],
        results['peak_records'],
        output_image=args.out_plot
    )

    if args.out_txt:
        with open(args.out_txt, 'w', encoding='utf-8') as f:
            f.write("# Peak#\t2-Theta(deg)\td-spacing(Angstrom)\tIntensity\tRelative_Intensity(%)\n")
            for idx, r in enumerate(results['peak_records'], 1):
                f.write(f"{idx}\t{r['2theta_deg']:.4f}\t{r['d_spacing_A']:.4f}\t{r['intensity']:.2f}\t{r['relative_intensity_pct']:.2f}\n")
        print(f"[+] Peak table exported: {args.out_txt}")


if __name__ == '__main__':
    main()
