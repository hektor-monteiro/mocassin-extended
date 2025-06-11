#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 13:27:01 2025

@author: hmonteiro
"""

import numpy as np
import matplotlib.pyplot as plt
import cmath # For complex mathematical operations like log

# --- Configuration ---
# IMPORTANT: Replace 'your_nk_file.txt' with the actual path to your N/K data file.
FILE_PATH = '/home/hmonteiro/Downloads/mocassin-rw_changes_2023/dustData/forst-cryst-perp1-kyoto-uvx.nk'
# IMPORTANT: Adjust the particle radius (in microns) as needed for your analysis.
PARTICLE_RADIUS_MICRONS = 1. # Example: 0.1 microns

# --- Constants ---
PI = cmath.pi # Use cmath.pi for consistency with complex numbers
# 1 micron = 1e-4 cm

# --- Mathematical Functions for Mie and CDE Theory ---

def calculate_mie_cabs(k, V, alpha):
    """
    Calculates Mie absorption cross-section (Cabs).
    Cabs = k * V * imag(3*alpha)
    """
    return k * V * alpha.imag * 3

def calculate_mie_csca(k, V, alpha):
    """
    Calculates Mie scattering cross-section (Csca).
    Csca = k^4 * V^2 * 3 * absolute(alpha)^2 / (2*pi)
    Note: Original email had 2*pi in denominator.
    """
    return (k**4) * (V**2) * 3 * (abs(alpha)**2) / (2 * PI)

def calculate_cde_cabs(k, V, alpha):
    """
    Calculates CDE absorption cross-section (Cabs).
    Cabs = k * V * alpha.imag
    """
    return k * V * alpha.imag

def calculate_cde_csca(k, V, alpha, sigma):
    """
    Calculates CDE scattering cross-section (Csca).
    Csca = k * V * imag(alpha) / sigma
    Handles potential division by zero, NaN, or Inf in sigma using an epsilon check.
    """
    # Define a small epsilon to check if sigma is effectively zero.
    # This helps catch floating-point numbers that are very close to zero but not exactly zero.
    EPSILON = 1e-18

    if np.isnan(sigma) or np.isinf(sigma) or abs(sigma) < EPSILON:
        # Return NaN if sigma is problematic (NaN, Inf, or effectively zero)
        # This will prevent the ZeroDivisionError and the NaN values will not be plotted.
        return np.nan
    return k * V * alpha.imag / sigma

# --- Main Data Processing and Calculation ---

def process_nk_file_and_calculate(file_path, particle_radius_microns):
    """
    Reads the N/K file, performs Mie and CDE calculations, and returns results.
    """
    wavelengths_micron = []
    n_values = []
    k_values = []

    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()

        # Skip header lines (first two lines)
        data_lines = lines[2:]

        for line in data_lines:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    wavelengths_micron.append(float(parts[0]))
                    n_values.append(float(parts[1]))
                    k_values.append(float(parts[2]))
                except ValueError as e:
                    print(f"Skipping malformed data line: {line.strip()} - Error: {e}")
            else:
                print(f"Skipping incomplete data line: {line.strip()}")

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return None, None, None, None, None, None, None, None, None
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None, None, None, None, None, None, None, None, None

    if not wavelengths_micron:
        print("No valid data found in the file after skipping headers.")
        return None, None, None, None, None, None, None, None, None

    # Convert particle radius to cm for volume calculation
    particle_radius_cm = particle_radius_microns * 1e-4 # 1 micron = 1e-4 cm
    # Assuming spherical particle volume
    V = (4/3) * PI * (particle_radius_cm**3)
    # Calculate geometrical cross-section (Ag)
    Ag = PI * (particle_radius_cm**2)


    # Lists to store calculated cross-sections and efficiencies
    mie_cabs_list = []
    mie_csca_list = []
    cde_cabs_list = []
    cde_csca_list = []
    mie_qabs_list = [] # New list for Mie Absorption Efficiency
    cde_qabs_list = [] # New list for CDE Absorption Efficiency


    for i in range(len(wavelengths_micron)):
        lam_micron = wavelengths_micron[i]
        n_val = n_values[i]
        k_val = k_values[i]

        # Convert wavelength from microns to cm for wavenumber calculation
        lam_cm = lam_micron * 1e-4 # 1 micron = 1e-4 cm
        k_wavenumber = (2 * PI) / lam_cm # k is wavenumber (2pi/lambda(cm))

        # Complex refractive index 'm'
        m = complex(n_val, k_val)

        # --- Mie Theory Calculations ---
        # alpha for Mie theory
        mie_alpha = (m**2 - 1) / (m**2 + 2)
        mie_cabs = calculate_mie_cabs(k_wavenumber, V, mie_alpha)
        mie_csca = calculate_mie_csca(k_wavenumber, V, mie_alpha)
        mie_cabs_list.append(mie_cabs)
        mie_csca_list.append(mie_csca)
        # Calculate Mie Absorption Efficiency (Qabs)
        mie_qabs_list.append(mie_cabs / Ag if Ag != 0 else np.nan)


        # --- CDE Theory Calculations ---
        # alpha for CDE theory
        # Handle log(m^2) carefully, especially if m^2 is zero or negative real.
        # cmath.log handles principal value (n=0) by default.
        try:
            cde_alpha = (2 * m**2 * cmath.log(m**2)) / (m**2 - 1) - 2
        except ZeroDivisionError:
            print(f"Warning: Division by zero when calculating CDE alpha for m^2-1 at wavelength {lam_micron} micron. Setting alpha to NaN.")
            cde_alpha = complex(np.nan, np.nan)
        except Exception as e:
            print(f"Warning: Error calculating CDE alpha for wavelength {lam_micron} micron: {e}. Setting alpha to NaN.")
            cde_alpha = complex(np.nan, np.nan)

        # sigma for CDE theory
        # Handle potential zero in denominator abs(m^2 - 1)^2
        try:
            sigma_denominator = (abs(m**2 - 1))**2
            if sigma_denominator == 0:
                print(f"Warning: Division by zero when calculating CDE sigma for abs(m^2-1)^2 at wavelength {lam_micron} micron. Setting sigma to NaN.")
                sigma = np.nan
            else:
                sigma = (6 * PI / (V * k_wavenumber**3)) * (m**2).imag / sigma_denominator
        except ZeroDivisionError: # This might not catch all cases, explicit check added above
            print(f"Warning: Division by zero when calculating CDE sigma (V*k^3) at wavelength {lam_micron} micron. Setting sigma to NaN.")
            sigma = np.nan
        except Exception as e:
            print(f"Warning: Error calculating CDE sigma for wavelength {lam_micron} micron: {e}. Setting sigma to NaN.")
            sigma = np.nan

        cde_cabs = calculate_cde_cabs(k_wavenumber, V, cde_alpha)
        cde_csca = calculate_cde_csca(k_wavenumber, V, cde_alpha, sigma)
        cde_cabs_list.append(cde_cabs)
        cde_csca_list.append(cde_csca)
        # Calculate CDE Absorption Efficiency (Qabs)
        cde_qabs_list.append(cde_cabs / Ag if Ag != 0 else np.nan)


    return wavelengths_micron, mie_cabs_list, mie_csca_list, cde_cabs_list, cde_csca_list, mie_qabs_list, cde_qabs_list, V, k_wavenumber # Return k_wavenumber for reference

# --- Plotting Function ---

def plot_comparison(wavelengths, mie_cabs, mie_csca, cde_cabs, cde_csca, mie_qabs, cde_qabs):
    """
    Plots the comparison of Mie and CDE cross-sections and efficiencies.
    """
    plt.figure(figsize=(18, 6)) # Adjust figure size for better readability (wider for 3 plots)

    # Plot Cabs
    plt.subplot(1, 3, 1) # 1 row, 3 columns, first plot
    plt.plot(wavelengths, mie_cabs, label='Mie $C_{abs}$', color='blue', linestyle='-')
    plt.plot(wavelengths, cde_cabs, label='CDE $C_{abs}$', color='red', linestyle='--')
    plt.title('Absorption Cross-Section ($C_{abs}$) Comparison')
    plt.xlabel('Wavelength (microns)')
    plt.ylabel('Cross-Section')
    plt.legend()
    plt.grid(True)
    plt.yscale('log') # Often cross-sections vary widely, log scale can be useful
    plt.xscale('log')
    plt.minorticks_on()

    # Plot Csca
    plt.subplot(1, 3, 2) # 1 row, 3 columns, second plot
    plt.plot(wavelengths, mie_csca, label='Mie $C_{sca}$', color='blue', linestyle='-')
    plt.plot(wavelengths, cde_csca, label='CDE $C_{sca}$', color='red', linestyle='--')
    plt.title('Scattering Cross-Section ($C_{sca}$) Comparison')
    plt.xlabel('Wavelength (microns)')
    plt.ylabel('Cross-Section')
    plt.legend()
    plt.grid(True)
    plt.yscale('log') # Often cross-sections vary widely, log scale can be useful
    plt.xscale('log')
    plt.minorticks_on()

    # Plot Qabs (Absorption Efficiency)
    plt.subplot(1, 3, 3) # 1 row, 3 columns, third plot
    plt.plot(wavelengths, mie_qabs, label='Mie $Q_{abs}$', color='blue', linestyle='-')
    plt.plot(wavelengths, cde_qabs, label='CDE $Q_{abs}$', color='red', linestyle='--')
    plt.title('Absorption Efficiency ($Q_{abs}$) Comparison')
    plt.xlabel('Wavelength (microns)')
    plt.ylabel('Efficiency (dimensionless)')
    plt.legend()
    plt.grid(True)
    # Qabs is typically dimensionless and can be <= 1, so log scale might not always be necessary
    # but can be useful if values vary widely. Let's keep it for consistency for now.
    #plt.yscale('log')
    #plt.xscale('log')
    plt.xlim(5,40)
    plt.minorticks_on()

    plt.tight_layout() # Adjust layout to prevent overlapping
    plt.show()

# --- Execution ---
if __name__ == "__main__":
    wavelengths, mie_cabs, mie_csca, cde_cabs, cde_csca, mie_qabs, cde_qabs, V, k_wavenumber_ref = \
        process_nk_file_and_calculate(FILE_PATH, PARTICLE_RADIUS_MICRONS)

    if wavelengths is not None:
        print(f"\n--- Calculation Summary ---")
        print(f"Particle Radius: {PARTICLE_RADIUS_MICRONS} microns")
        print(f"Calculated Particle Volume (V): {V:.2e} cm^3")
        print(f"Calculated Geometrical Cross-Section (Ag): {PI * (PARTICLE_RADIUS_MICRONS * 1e-4)**2:.2e} cm^2") # Print Ag for reference

        # Display first few calculated values for verification
        print("\nFirst 5 Mie Cabs:", [f"{x:.2e}" for x in mie_cabs[:5]])
        print("First 5 CDE Cabs:", [f"{x:.2e}" for x in cde_cabs[:5]])
        print("First 5 Mie Csca:", [f"{x:.2e}" for x in mie_csca[:5]])
        print("First 5 CDE Csca:", [f"{x:.2e}" for x in cde_csca[:5]])
        print("First 5 Mie Qabs:", [f"{x:.2e}" for x in mie_qabs[:5]]) # New print for Qabs
        print("First 5 CDE Qabs:", [f"{x:.2e}" for x in cde_qabs[:5]]) # New print for Qabs

        plot_comparison(wavelengths, mie_cabs, mie_csca, cde_cabs, cde_csca, mie_qabs, cde_qabs)
