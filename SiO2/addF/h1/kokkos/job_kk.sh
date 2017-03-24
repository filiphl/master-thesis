#!/bin/bash
#SBATCH --partition=kif
#SBATCH --ntasks=4
#SBATCH --time=01:00:00
#SBATCH --gres=gpu:4
#SBATCH --job-name=normal
echo $CUDA_VISIBLE_DEVICES
mpirun lmp_kokkos_cuda_openmpi -echo both -k on g 4 -sf kk -in system.run
