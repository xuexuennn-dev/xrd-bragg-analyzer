"""
XRD Bragg Analyzer
==================
A lightweight scientific computing toolkit for X-ray diffraction (XRD) data processing,
candidate peak detection, Bragg's Law interplanar lattice spacing (d-spacing) calculation,
and crystallographic profile visualization.

Maintainer: @xuexuennn-dev (github.com/xuexuennn-dev)
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


def analyze_xrd(two_theta, intensity, wavelength=1.5406, prominence=None, distance=20, smooth=True):
    """
    Execute XRD processing pipeline:
      1. Noise smoothing with Savitzky-Golay filter
      2. Peak detection via scipy.signal.find_peaks (adaptive or user-specified prominence)
      3. d-spacing calculation using Bragg's Law
      4. Net relative intensity (I/I0) calculation
    """
    intensity = np.array(intensity, dtype=float)
    two_theta = np.array(two_theta, dtype=float)

    if smooth and len(intensity) > 15:
        # Window length 11, polynomial order 3
        smoothed_intensity = savgol_filter(intensity, window_length=11, polyorder=3)
        smoothed_intensity = np.maximum(smoothed_intensity, 0)
    else:
        smoothed_intensity = intensity

    # Adaptive prominence if not explicitly specified:
    # 5% of dynamic range (max - min) of smoothed intensity, minimum 20.0
    dynamic_range = float(np.max(smoothed_intensity) - np.min(smoothed_intensity))
    if prominence is None:
        effective_prominence = max(20.0, 0.05 * dynamic_range)
    else:
        effective_prominence = float(prominence)

    # Peak detection
    peaks, _ = find_peaks(
        smoothed_intensity,
        prominence=effective_prominence,
        distance=distance
    )

    peak_2theta = two_theta[peaks]
    peak_intensity = smoothed_intensity[peaks]
    d_spacings = calculate_d_spacing(peak_2theta, wavelength=wavelength)

    # Net relative intensity: subtract baseline minimum to reflect true peak signal
    baseline = float(np.min(smoothed_intensity))
    net_intensity = peak_intensity - baseline
    max_net = float(np.max(net_intensity)) if len(net_intensity) > 0 and np.max(net_intensity) > 0 else 1.0
    relative_intensity = (net_intensity / max_net) * 100.0

    peak_records = []
    for tt, d, counts, rel_i in zip(peak_2theta, d_spacings, peak_intensity, relative_intensity):
        peak_records.append({
            '2theta_deg': float(tt),
            'd_spacing_A': float(d),
            'intensity': float(counts),
            'relative_intensity_pct': float(rel_i)
        })

    return {
        'two_theta': two_theta,
        'raw_intensity': intensity,
        'smoothed_intensity': smoothed_intensity,
        'peaks': peaks,
        'peak_records': peak_records,
        'prominence_used': effective_prominence
    }


def plot_diffractogram(two_theta, raw_intensity, smoothed_intensity, peak_records, output_image="xrd_plot.png", max_labeled_peaks=10):
    """
    Generate publication-ready XRD diffractogram with clear, non-overlapping peak annotations.
    """
    plt.figure(figsize=(9, 5.5), dpi=300)
    plt.plot(two_theta, raw_intensity, color='#94a3b8', alpha=0.45, linewidth=0.8, label='Raw Intensity')
    plt.plot(two_theta, smoothed_intensity, color='#0f172a', linewidth=1.2, label='Smoothed Profile')

    # Filter candidate peaks for annotation to prevent clutter:
    # Select peaks with relative intensity >= 15%, up to max_labeled_peaks
    annotated_peaks = [rec for rec in peak_records if rec['relative_intensity_pct'] >= 15.0]
    if len(annotated_peaks) > max_labeled_peaks:
        annotated_peaks = sorted(annotated_peaks, key=lambda x: x['relative_intensity_pct'], reverse=True)[:max_labeled_peaks]
    # Re-sort by 2theta for consistent left-to-right presentation
    annotated_peaks = sorted(annotated_peaks, key=lambda x: x['2theta_deg'])

    y_max = float(np.max(raw_intensity))
    y_min = float(np.min(raw_intensity))
    y_range = y_max - y_min if y_max > y_min else 1.0

    for i, rec in enumerate(annotated_peaks):
        tt = rec['2theta_deg']
        d = rec['d_spacing_A']
        inte = rec['intensity']
        
        # Stagger label heights to avoid overlap
        offset_y = 12 if (i % 2 == 0) else 26
        
        plt.scatter(tt, inte, color='#dc2626', s=28, zorder=5)
        plt.annotate(
            f"{tt:.2f}°\n({d:.3f}Å)",
            xy=(tt, inte),
            xytext=(0, offset_y),
            textcoords="offset points",
            ha='center',
            fontsize=8,
            color='#1e293b',
            weight='bold',
            arrowprops=dict(arrowstyle="-", color='#dc2626', alpha=0.45, lw=0.7)
        )

    plt.xlabel(r'Diffraction Angle $2\theta$ (degrees)', fontsize=11, fontweight='bold')
    plt.ylabel('Intensity (counts)', fontsize=11, fontweight='bold')
    plt.title(r'X-Ray Diffractogram & Candidate Peaks ($\lambda = 1.5406\ \mathrm{\AA}$)', fontsize=12, fontweight='bold')
    plt.legend(frameon=True, facecolor='#ffffff', edgecolor='#cbd5e1')
    plt.grid(True, linestyle='--', alpha=0.35)
    plt.ylim(y_min - 0.05 * y_range, y_max + 0.16 * y_range)
    plt.tight_layout()
    plt.savefig(output_image)
    plt.close()
    print(f"[+] Diffractogram exported: {output_image}")


def plot_xrd(results_or_two_theta, *args, save_path="xrd_analysis_plot.png", **kwargs):
    """
    Convenient plotting function supporting both results dict and raw array inputs:
      Usage 1: plot_xrd(results, save_path="plot.png")
      Usage 2: plot_xrd(two_theta, raw_intensity, smoothed_intensity, peak_records, output_image="plot.png")
    """
    if isinstance(results_or_two_theta, dict):
        res = results_or_two_theta
        plot_diffractogram(
            two_theta=res['two_theta'],
            raw_intensity=res['raw_intensity'],
            smoothed_intensity=res['smoothed_intensity'],
            peak_records=res['peak_records'],
            output_image=kwargs.get('output_image', save_path),
            max_labeled_peaks=kwargs.get('max_labeled_peaks', 10)
        )
    else:
        plot_diffractogram(results_or_two_theta, *args, output_image=save_path, **kwargs)


def print_summary_table(peak_records):
    """
    Print formatted candidate peak summary.
    """
    print("\n" + "=" * 65)
    print(f"{'Peak #':<8}{'2-Theta (°)':<15}{'d-spacing (Å)':<18}{'Counts':<12}{'I/I0 (%)':<10}")
    print("-" * 65)
    for idx, r in enumerate(peak_records, 1):
        print(f"{idx:<8}{r['2theta_deg']:<15.3f}{r['d_spacing_A']:<18.4f}{r['intensity']:<12.1f}{r['relative_intensity_pct']:<10.1f}")
    print("=" * 65 + "\n")


def main():
    parser = argparse.ArgumentParser(description="XRD Bragg Analyzer - Candidate Peak Detection & d-spacing Calculator")
    parser.add_argument("input_file", help="Path to XRD data file (.raw binary or .txt/.csv ASCII)")
    parser.add_argument("--wavelength", "-w", type=float, default=1.5406, help="Radiation wavelength in Angstroms (default: 1.5406 A)")
    parser.add_argument("--prominence", "-p", type=float, default=None, help="Peak prominence threshold (default: adaptive 5%% of dynamic range)")
    parser.add_argument("--distance", "-d", type=int, default=20, help="Minimum peak distance in points (default: 20)")
    parser.add_argument("--output", "-o", "--out-plot", dest="out_plot", default="xrd_analysis_plot.png", help="Output plot filename")
    parser.add_argument("--out-txt", default=None, help="Optional text export for candidate peaks")

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

    print(f"[+] Detected {len(results['peak_records'])} candidate peaks (adaptive prominence = {results['prominence_used']:.1f})")
    print_summary_table(results['peak_records'])
    
    plot_xrd(results, save_path=args.out_plot)

    if args.out_txt:
        with open(args.out_txt, 'w', encoding='utf-8') as f:
            f.write("# Peak#\t2-Theta(deg)\td-spacing(Angstrom)\tIntensity\tRelative_Intensity(%)\n")
            for idx, r in enumerate(results['peak_records'], 1):
                f.write(f"{idx}\t{r['2theta_deg']:.4f}\t{r['d_spacing_A']:.4f}\t{r['intensity']:.2f}\t{r['relative_intensity_pct']:.2f}\n")
        print(f"[+] Peak table exported: {args.out_txt}")


if __name__ == '__main__':
    main()
