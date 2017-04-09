#!/bin/sh
# Job name:
#SBATCH --job-name=shear
#Project:
#SBATCH --account=uio
##SBATCH --account=trocks
# Wall clock limit:
#SBATCH --time='1-00:00:00'
# Max memory usage per task:
#SBATCH --mem-per-cpu=3600M
# Number of tasks (MPI ranks):
#SBATCH --ntasks=16
#SBATCH --mail-type=ALL
## Set up job environment
source /cluster/bin/jobsetup
#module load openmpi.gnu/2.0.1
#module load openmpi.gnu
module load intel
module load intelmpi.intel

## Run command
mpirun ~/smaug/lammps/src/lmp_intel_cpu_intelmpi -in system.run4
