# Density-functional-theory workflow

This directory contains an exploratory ASE/GPAW workflow for studying twist-dependent MoS2-Au interfacial structures. It is a research snapshot rather than a turnkey software package: the reference structure and primary relaxation input are included, while some structure-building and plotting scripts still point to uncommitted local inputs.

## Included reference calculation

`moire_supercell_MoS2-3DAu_30.cif` is a P1 supercell with:

- composition: `Mo56 S112 Au128` (296 atoms);
- cell: 22.1032 x 22.3322 x 28.0 Angstrom; and
- orthogonal cell angles.

`job.py` relaxes this structure with the following settings:

| Setting | Value in the script |
|---|---|
| Calculator | GPAW |
| Basis / cutoff | Plane waves, 600 eV |
| Exchange-correlation | vdW-DF-cx through `libvdwxc` |
| k-point mesh | 1 x 1 x 1 |
| Occupations | Fermi-Dirac, 0.03 eV |
| Eigensolver | Davidson(3) |
| Geometry optimizer | ASE LBFGS |
| Force stopping criterion | 0.05 eV/Angstrom |

The calculation writes `MoS2-3DAu_30.txt`, `MoS2-3DAu_30.log`, and `MoS2-3DAu_30.traj`.

## Files and roles

| File | Purpose | Checkout status |
|---|---|---|
| `moire_supercell_MoS2-3DAu_30.cif` | Reference MoS2-Au structure | Ready to read |
| `job.py` | GPAW geometry relaxation | Input included; GPAW environment required |
| `submit.sh` | SLURM launcher for `job.py` | Cluster-specific paths, modules, and resources must be edited |
| `rotate.py` | Rotates the Au atoms about the surface normal | Hard-coded working directory, rotation center, and angle must be checked |
| `merge.py` | Translates and combines independently stored layers | Required `energy calculation/Nov. 5` source files are not committed |
| `plot E_interfacial & E_tot.py` | Fits and plots twist-dependent interfacial-energy data | Required `plot.csv` is not committed |

## Software

The relaxation requires:

- Python;
- ASE;
- GPAW with PAW setup files; and
- a working `libvdwxc` installation for vdW-DF-cx.

The structure utilities require NumPy and ASE. The plotting script additionally imports pandas, Matplotlib, SciPy, and `xlrd`.

Because GPAW builds and launch configurations depend on the computing environment, install and validate GPAW separately before submitting a production calculation. Record the exact Python, ASE, GPAW, PAW-dataset, MPI, BLAS, FFTW, ScaLAPACK, and `libvdwxc` versions used for reproducible results.

## Suggested workflow

### 1. Inspect and prepare a structure

Start with the included CIF or adapt `merge.py` to combine your own MoS2 and Au structures. Before using `rotate.py`, replace its absolute `path`, verify `shift_vec` is the intended rotation center, and set the desired angle.

~~~bash
cd DFT
python3 rotate.py
~~~

The checked-in script uses a 38-degree rotation and writes `moire_supercell_MoS2-3DAu_.cif`. Treat these values as one recorded configuration, not general defaults.

### 2. Validate the electronic-structure setup

Before a production relaxation, check at minimum:

- cell vectors, vacuum, stoichiometry, periodic boundaries, and atomic overlaps;
- plane-wave cutoff, k-point sampling, smearing, and force tolerance;
- spin and relativistic treatment appropriate to the target observable; and
- convergence with respect to the Au and MoS2 slab sizes.

The included `1 x 1 x 1` k-point mesh and 600 eV cutoff document the original input; they are not a substitute for convergence tests.

### 3. Run a relaxation

For a configured local GPAW environment:

~~~bash
cd DFT
python3 job.py
~~~

For SLURM, edit `submit.sh` to match the cluster account, partition, module stack, GPAW setup path, Python environment, and resource policy, then submit:

~~~bash
sbatch submit.sh
~~~

### 4. Analyze twist-dependent energies

The plotting script expects a headerless CSV at the path encoded in the script, with twist angle in the first column and energy in the second. Update the path or stage the data before running:

~~~bash
python3 "plot E_interfacial & E_tot.py"
~~~

It fits the interfacial-energy points to a cosine series and evaluates a two-interface energy of the form `f(x) + f(m - x)`. The present script fixes selected total twist angles near 10 degrees; inspect these values before applying it to a different system.

## Important reproducibility notes

- `rotate.py`, `merge.py`, and the plotting script contain hard-coded paths and geometry-specific translations. Convert these to explicit arguments before treating the workflow as reusable.
- `merge.py` and the plotting script are not independently runnable from this checkout because their source structure and CSV files are absent.
- `submit.sh` records one historical cluster configuration and overrides `HOME`; remove that override or replace it with site-appropriate paths before use.
- Preserve raw input files, relaxed trajectories, calculator text output, convergence tests, and exact software versions alongside any reported energies.
- Compare total energies only between calculations with consistent cells, atom counts, boundary conditions, calculator settings, and numerical convergence.

For the continuum-scale confinement model and interactive parameter sweep, return to the [project README](../README.md) or [launch the explorer](https://cuiyist.github.io/MoS2-Au-MoS2/).
