#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 21:59:17 2025

@author: hmonteiro
"""

import numpy as np

def convert_draine_OptData_file(input_filename="callindex.out_sil.D03", output_filename="astrosil.nk"):
    """
    Converts a Draine's astronomical silicate optical constant file
    to a new file with wavelength, n, and k columns.

    Args:
        input_filename (str): The name of the input data file (e.g., "callindex.out_sil.D03").
        output_filename (str): The name of the output file (e.g., "astrosil_nk.txt").
    """
    try:
        # Read the header to find the number of wavelengths
        with open(input_filename, 'r') as f_in:
            lines = f_in.readlines()
            # Find the line containing "number of wavelengths"
            num_wavelengths_line = [line for line in lines if "number of wavelengths" in line][0]
            num_wavelengths = int(num_wavelengths_line.split('=')[0].strip())

        # Determine the starting line of the data by skipping header lines
        # We know there are 4 header lines before the column names, and then the data starts.
        # So, the data starts on the 6th line (index 5) if we count from 0.
        # Alternatively, find the line that starts with "wave(um)" and the next line is data.
        data_start_line = 0
        for i, line in enumerate(lines):
            if "wave(um)" in line:
                data_start_line = i + 1
                break
        
        # Load the data using numpy.genfromtxt, skipping the header lines
        # We need to skip data_start_line lines.
        data = np.genfromtxt(input_filename, skip_header=data_start_line)

        # Extract the relevant columns
        # Column 0: wavelength (um)
        # Column 3: Re(n)-1
        # Column 4: Im(n) which is k
        wavelength = data[:, 0]
        re_n_minus_1 = data[:, 3]
        im_n = data[:, 4]

        # Calculate n by adding 1 to Re(n)-1
        n = re_n_minus_1 + 1.0
        k = im_n # Im(n) is already k

        # Combine into a new array for saving
        output_data = np.column_stack((wavelength, n, k))

        # Save to the new file
        header = "wavelength        n          k"
        np.savetxt(output_filename, output_data, fmt='%.6e', header=header, comments='# ')

        print(f"Successfully converted '{input_filename}' to '{output_filename}'.")

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    
    # Example usage:
    convert_draine_OptData_file(input_filename="callindex.out_silD03", output_filename="astroSilD2003.nk")

