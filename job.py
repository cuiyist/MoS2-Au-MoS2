from gpaw import GPAW, PW, FermiDirac
from gpaw.xc.libvdwxc import vdw_df_cx
from ase.optimize.lbfgs import LBFGS
from ase.io import read
from gpaw import Davidson

# define a structure
mos2 = read("moire_supercell_MoS2-3DAu_30.cif")

# define a calculator
calc = GPAW(mode=PW(600),
            xc=vdw_df_cx(),
            convergence={'eigenstates': 1e-8},
            eigensolver=Davidson(3),
            #eigensolver='cg',
            dtype=complex,
            kpts=(1, 1, 1),
            parallel={'augment_grids': True, 'sl_auto': True},
            occupations=FermiDirac(0.03),
            txt='MoS2-3DAu_30.txt')

mos2.calc = calc

# geometry optimization
opt = LBFGS(mos2, logfile="MoS2-3DAu_30.log", trajectory="MoS2-3DAu_30.traj")
opt.run(0.05)

