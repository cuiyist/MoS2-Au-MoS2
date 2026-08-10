# DFT structure rotation

This historical branch preserves the original utility used to rotate an Au layer relative to an MoS2 lattice. It contains one extensionless Python script, `rotate`.

The more complete DFT snapshot - including a reference CIF, GPAW relaxation input, SLURM launcher, structure utilities, and energy-fitting script - is documented in [`DFT/` on the `main` branch](https://github.com/cuiyist/MoS2-Au-MoS2/tree/main/DFT).

## What the script does

`rotate`:

1. reads `moire_supercell_MoS2-3DAu_30.cif` with ASE;
2. identifies atoms whose chemical symbol is `Au`;
3. translates the Au coordinates to a selected rotation center;
4. rotates them about the surface-normal z axis;
5. translates them back; and
6. writes `moire_supercell_MoS2-3DAu_.cif`.

The checked-in configuration applies a 38-degree rotation. The rotation matrix is generated from an axis-angle representation. Only Au atomic coordinates are changed; the cell and MoS2 coordinates are left unchanged.

## Requirements

- Python
- NumPy
- ASE

For example, in an isolated Python environment:

~~~bash
python3 -m pip install numpy ase
~~~

## Before running

This branch does not include the input CIF and the script is not portable without editing. Update the following values in `rotate`:

- `path` - replace the original absolute working directory;
- input CIF filename - point it to the intended MoS2-Au structure;
- `shift_vec` - verify the in-plane center of rotation;
- `theta` - set the intended angle and sign convention; and
- output filename - avoid overwriting an existing structure.

The input structure is available in the [maintained DFT directory on `main`](https://github.com/cuiyist/MoS2-Au-MoS2/blob/main/DFT/moire_supercell_MoS2-3DAu_30.cif).

After adapting the paths and parameters:

~~~bash
python3 rotate
~~~

## Validation

Before using a generated structure in an electronic-structure calculation:

- inspect the rotated cell visually;
- check minimum interatomic distances and periodic-boundary crossings;
- confirm the rotation center and angle convention;
- wrap atoms into the periodic cell if required by the downstream workflow; and
- retain the exact input structure and rotation parameters with the calculation record.

This script prepares geometry only; it does not run DFT or calculate an energy. See the [main project README](https://github.com/cuiyist/MoS2-Au-MoS2) for the confinement-pressure analysis and the fuller DFT workflow.
