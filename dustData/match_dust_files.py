#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 26 13:25:29 2025

@author: hmonteiro
"""

import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

############################################################
def interpolate_dust_properties(base_particle, target_particle):
    """
    Interpolates the optical constants (n and k) of a target dust particle
    onto the wavelength grid of a base particle.

    If the target particle's wavelength range is smaller than the base particle's,
    the base particle's n and k values are used for the missing ranges.

    Args:
        base_particle (np.ndarray): A NumPy array with shape (N, 3)
                                    representing [lambda, n, k] for the base particle.
        target_particle (np.ndarray): A NumPy array with shape (M, 3)
                                      representing [lambda, n, k] for the target particle to be interpolated.

    Returns:
        np.ndarray: A new NumPy array with shape (N, 3) containing the interpolated
                    n and k values on the base particle's lambda grid.
    """
    # Unpack the arrays for clarity
    lambda_base = base_particle[:, 0]
    n_base = base_particle[:, 1]
    k_base = base_particle[:, 2]

    lambda_target = target_particle[:, 0]
    n_target = target_particle[:, 1]
    k_target = target_particle[:, 2]

    # Perform linear interpolation of n and k values from the target particle
    # onto the base particle's wavelength grid.
    # For wavelengths outside the target's range, np.interp will use the
    # specified 'left' and 'right' arguments. We set these to be the
    # corresponding values from the base particle itself.
    n_interpolated = np.interp(
        lambda_base,
        lambda_target,
        n_target,
        left=n_base[0],
        right=n_base[-1]
    )

    k_interpolated = np.interp(
        lambda_base,
        lambda_target,
        k_target,
        left=k_base[0],
        right=k_base[-1]
    )
    
    # check intervals and adjust accordingly
    # if target min > base min
    if (lambda_target[0] > lambda_base[0]):
        ind = np.where(lambda_base < lambda_target[0])
        n_interpolated[ind] = (n_base+n_interpolated)[ind]/2
        k_interpolated[ind] = (k_base+k_interpolated)[ind]/2
        
    # if target max < base max
    if (lambda_target[-1] < lambda_base[-1]):
        ind = np.where(lambda_base > lambda_target[-1])
        n_interpolated[ind] = (n_base+n_interpolated)[ind]/2
        k_interpolated[ind] = (k_base+k_interpolated)[ind]/2

    # Combine the base lambda values with the new interpolated n and k values.
    # The result is stacked horizontally to create the final [lambda, n, k] array.
    n_interpolated = (n_interpolated)
    interpolated_properties = np.vstack((lambda_base, n_interpolated, k_interpolated)).T

    return interpolated_properties

    
###############################################################################
# Base file of Astrophysical Silicates

astrosilfile = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/astroSilD2003.nk'
astrosildata = np.loadtxt(astrosilfile,max_rows=837,skiprows=2)

# increase sampling of astro silicate data
new_lambda = np.logspace(np.log10(astrosildata[:,0].min()),np.log10(astrosildata[:,0].max()), 2000)
new_Sil_n = np.interp(new_lambda,astrosildata[:,0],astrosildata[:,1])
new_Sil_k = np.interp(new_lambda,astrosildata[:,0],astrosildata[:,2])

astrosildata = np.vstack((new_lambda,new_Sil_n,new_Sil_k)).T

header = "nk\nSsil_dl_base  1400. 3.6  0.588 20.077"
output_file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/astroSilD2003_base.nk'
np.savetxt(output_file, astrosildata, header=header, comments='', fmt='%.6E')

###############################################################################
# Water ice from Warren and Brandt 2008: Ice; n,k 0.00443–2e6 µm
# S. G. Warren and R. E. Brandt. Optical constants of ice from the ultraviolet to the microwave: A revised compilation. J. Geophys. Res. 113, D14220 (2008)

h2ofile08 = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/h20_Warren08.nk'
h2O08data = np.loadtxt(h2ofile08,max_rows=486,skiprows=2)

# Call the function to get the interpolated data
new_dust_prop = interpolate_dust_properties(astrosildata, h2O08data)

plt.figure()
plt.xscale('log')
# plt.yscale('log')
plt.title('H2O ice')
plt.xlabel('wavelenth (microns)')

plt.plot(astrosildata[:,0], astrosildata[:,1],'C0',label='astrosil n')
plt.plot(astrosildata[:,0], astrosildata[:,2],'C0',label='astrosil k')

plt.plot(h2O08data[:,0], h2O08data[:,1],'C1--',label='h2O ice n')
plt.plot(h2O08data[:,0], h2O08data[:,2],'C1--',label='h2O ice k')

plt.plot(new_dust_prop[:,0], new_dust_prop[:,1],'C2--',label='h2O ice new n')
plt.plot(new_dust_prop[:,0], new_dust_prop[:,2],'C2--',label='h2O ice new k')

plt.legend()

header = "nk\nSh2O_ice 273 0.917 0.588 6.005"
output_file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/H2O_ice.nk'
np.savetxt(output_file, new_dust_prop, header=header, comments='', fmt='%.6E')

###############################################################################
# Dity Ice from https://ui.adsabs.harvard.edu/abs/1993A%26A...279..577P/abstract

dirtyfile = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/dirty-ice.nk'
dirtydata = np.loadtxt(dirtyfile,max_rows=90,skiprows=2)

# Call the function to get the interpolated data
new_dust_prop = interpolate_dust_properties(astrosildata, dirtydata)

plt.figure()
plt.xscale('log')
# plt.yscale('log')
plt.title('dirty H2O ice')
plt.xlabel('wavelenth (microns)')

plt.plot(astrosildata[:,0], astrosildata[:,1],'C0',label='astrosil n')
plt.plot(astrosildata[:,0], astrosildata[:,2],'C0',label='astrosil k')

plt.plot(dirtydata[:,0], dirtydata[:,1],'C1--',label='dirty h2O ice n')
plt.plot(dirtydata[:,0], dirtydata[:,2],'C1--',label='dirty h2O ice k')

plt.plot(new_dust_prop[:,0], new_dust_prop[:,1],'C2--',label='dirty h2O ice new n')
plt.plot(new_dust_prop[:,0], new_dust_prop[:,2],'C2--',label='dirty h2O ice new k')

plt.legend()


header = "nk\nSh2O_dirtyice 175 1.6 0.588 6.005"
output_file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/H2O_dirtyice.nk'
np.savetxt(output_file, new_dust_prop, header=header, comments='', fmt='%.6E')

###############################################################################
# Optical constants for forsterite (Mg2SiO4) - crystalline
# Sogawa et al. 2006, Astron. Astrophys. 451, 357

forstp1file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/forst-cryst-perp1-kyoto-uvx.nk'
forstp2file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/forst-cryst-perp2-kyoto-uvx.nk'
forstp3file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/forst-cryst-perp3-kyoto-uvx.nk'

forstp1data = np.loadtxt(forstp1file,max_rows=6769,skiprows=2)
forstp2data = np.loadtxt(forstp2file,max_rows=6783,skiprows=2)
forstp3data = np.loadtxt(forstp3file,max_rows=6814,skiprows=2)


# Call the function to get the interpolated data
fost1 = interpolate_dust_properties(astrosildata, forstp1data)
fost2 = interpolate_dust_properties(astrosildata, forstp2data)
fost3 = interpolate_dust_properties(astrosildata, forstp3data)

# combine polarizations
forstdata = (fost1+fost2+fost3)/3

# plot results
plt.figure()
plt.xscale('log')
# plt.yscale('log')
plt.title('Fosterite')
plt.xlabel('wavelenth (microns)')

plt.plot(astrosildata[:,0], astrosildata[:,1],'C0',label='astrosil n')
plt.plot(astrosildata[:,0], astrosildata[:,2],'C0',label='astrosil k')

plt.plot(forstp1data[:,0], forstp1data[:,1],'C1--',label='Fost. crystal n')
plt.plot(forstp1data[:,0], forstp1data[:,2],'C1--',label='Fost. crystal k')

plt.plot(forstdata[:,0], forstdata[:,1],'C2--',label='new n')
plt.plot(forstdata[:,0], forstdata[:,2],'C2--',label='new k')

plt.legend()


header = 'nk\nSforst-cryst-kyoto-uvx 1000 3.27 0.588 20.10'
output_file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/forsterite.nk'
np.savetxt(output_file, forstdata, header=header, comments='', fmt='%.6E')

###############################################################################
# Optical constants for enstatite (MgSiO3) - crystalline
# Jager et al. 1998, Astron. Astrophys. 339, 904

enstap1file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/enst-cryst-par-jena-uvx.nk'
enstap2file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/enst-cryst-perp1-jena-uvx.nk'
enstap3file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/enst-cryst-perp2-jena-uvx.nk'

enstap1data = np.loadtxt(enstap1file,max_rows=1651,skiprows=2)
enstap2data = np.loadtxt(enstap2file,max_rows=1652,skiprows=2)
enstap3data = np.loadtxt(enstap3file,max_rows=1652,skiprows=2)


# Call the function to get the interpolated data
ensta1 = interpolate_dust_properties(astrosildata, enstap1data)
ensta2 = interpolate_dust_properties(astrosildata, enstap2data)
ensta3 = interpolate_dust_properties(astrosildata, enstap3data)

# combine polarizations
enstadata = (ensta1+ensta2+ensta3)/3

# plot results
plt.figure()
plt.xscale('log')
# plt.yscale('log')
plt.title('Enstatite')
plt.xlabel('wavelenth (microns)')

plt.plot(astrosildata[:,0], astrosildata[:,1],'C0',label='astrosil n')
plt.plot(astrosildata[:,0], astrosildata[:,2],'C0',label='astrosil k')

plt.plot(enstap1data[:,0], enstap1data[:,1],'C1--',label='ensta. crystal n')
plt.plot(enstap1data[:,0], enstap1data[:,2],'C1--',label='ensta. crystal k')

plt.plot(enstadata[:,0], enstadata[:,1],'C2--',label='new n')
plt.plot(enstadata[:,0], enstadata[:,2],'C2--',label='new k')

plt.legend()


header = 'nk\nSenst-cryst-jena-uvx 1000 3.2 0.588 20.08'
output_file = '/home/hmonteiro/Google Drive/work/PN/software/mocassin_HM_2025/dustData/new_dust_files/enstatite.nk'
np.savetxt(output_file, enstadata, header=header, comments='', fmt='%.6E')


