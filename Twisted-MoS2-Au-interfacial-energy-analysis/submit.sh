#!/bin/bash
#
#SBATCH --job-name=MoS2-3D_30
#SBATCH --partition=normal
#SBATCH -N 6
#SBATCH --time=48:00:00
#SBATCH --ntasks=60
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5G

export HOME=/home/users/cuiy
export GPAW_SETUP_PATH=$HOME/gpaw-setups-0.9.20000

module load math
module load openmpi/4.1.2
module load openblas/0.3.10
module load fftw/3.3.10
module load scalapack/2.2.0

srun $HOME/miniconda3/envs/gpaw/bin/python job.py 2> error.txt
