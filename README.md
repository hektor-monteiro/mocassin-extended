# Mocassin-Extended: Monte Carlo Simulations of Ionised Nebulae

**MOCASSIN** is a fully 3D or 2D photoionisation and dust radiative transfer code which employs a Monte Carlo approach to the transfer of radiation through media of arbitrary geometry and density distribution. It was originally developed by Barbara Ercolano (https://mocassin.nebulousresearch.org/publications) for the modelling of photoionised regions like HII regions and planetary nebulae and has since expanded and been applied to a variety of astrophysical problems, including modelling clumpy dusty supernova envelopes, star forming galaxies, protoplanetary disks and inner shell fluorescence emission in the photospheres of stars and disk atmospheres.

The code can deal with arbitrary Cartesian grids of variable resolution, it has successfully been used to model complex density fields from SPH calculations and can deal with ionising radiation extending from Lyman edge to the X-ray. The dust and gas microphysics is fully coupled both in the radiation transfer and in the thermal balance.

The code is detailed in https://mocassin.nebulousresearch.org/

This repository is a fork of the official Mocassin repo located at: https://github.com/rwesson/mocassin

## Overview

The project provides four main executables:

- `mocassin`: The main MOCASSIN driver that is used to start a new simulation.
- `mocassinWarm`: Resumes an interrupted simulation using existing grid files.
- `mocassinOutput`: Runs the output routines using the current grid files in the output subdirectory.
- `mocassinPlot`: Uses the current grid files in the output/ subdirectory to create 3d-emission maps.

## Building and Installation

Building the project requires MPI and Fortran compilers. On Debian/Ubuntu, install them via:
```bash
sudo apt-get install -y openmpi-bin libopenmpi-dev gfortran
```

The project is built using `make`:
```bash
make clean
make
```
Compiler flags can be overridden, for example: `make FCFLAGS="-Wall -Wextra"`.

To install locally:
```bash
make install
```


## Running the Benchmarks

### Directory Setup
To run the models, the code must be properly installed wit the `data` and `dustData` directories copied or linked to `/usr/share/mocassin/` or available locally depending on your environment. 

### Pure Photoionisation Benchmarks
You can try to run one or more of the available benchmark problems located in the `benchmarks/` directory.

To run the standard test problem (Meudon standard HII region), navigate to `benchmarks/test_Problem`:
```bash
cd benchmarks/test_Problem
mkdir -p output/
mpirun -np 1 ../../mocassin
```

Alternatively, to setup a benchmark manually:
1. Copy or link the benchmark `input.in` to your working `input` directory.
2. Ensure you have the corresponding abundance file.
3. Run using MPI: `mpirun -np <num_procs> ./mocassin`

### Pure Dust Benchmarks
There are also 1D and 2D benchmark models included under `benchmarks/dust/1D` and `benchmarks/dust/2D`. Copy the required input file to `input/input.in` and execute MOCASSIN.

## Input and Output Files

- **Input files**: Placed in the `input/` directory (e.g. `input.in`, abundance files, density distribution files).
- **Data files**: Located in the `data/` directory (atomic data files) and `dustData/` (dust optical data library).
- **Output files**: Saved in the `output/` directory (e.g. `ionratio.out`, `lineFlux.out`, `temperature.out`, `SED.out`, `grid0.out`, `grid1.out`, `grid2.out`, `grid3.out`, `dustGrid.out`).

See `man/mocassin.1` for a complete list of keywords that can be configured in the `input.in` file.

## Licence and credits

Mocassin-Extended is distributed under the **GNU General Public License, version 3** ([LICENSE](LICENSE)).

- **MOCASSIN** is Copyright (C) 2003-2005 Barbara Ercolano and is licensed under GPL version 2 or any
  later version. This fork starts from R. Wesson's version ([rwesson/mocassin](https://github.com/rwesson/mocassin)).
  Files changed here keep their original copyright notice and carry a modification note.
Work that uses MOCASSIN-PDR should cite:

- MOCASSIN: Ercolano et al. 2003, MNRAS 340, 1136; Ercolano, Barlow & Storey 2005, MNRAS 362, 1038;
  Ercolano et al. 2008, ApJS 175, 534

