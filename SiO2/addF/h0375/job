#!/bin/bash
#SBATCH --partition=smaug-c
#SBATCH --ntasks=64
#SBATCH --time=5-01:00:00
#SBATCH --job-name=MD
mpirun ~/scratch/lammps/src/lmp_mpi -in system.run2
#python notify.py $SLURM_JOB_ID
