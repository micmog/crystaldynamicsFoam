#!/bin/bash --login
###SBATCH -p multicore  # (or --partition=) 
###SBATCH -n 2 # -c (or --cpus-per-task=) Number of cores to use for OpenMP (2--40)

module purge
module load openfoam/10-foss-2021a

source $FOAM_BASH
rm -r 0* 1* 2* 3* 5* 6* 7* 8* 9* #processor* 
cp -r initial 0
blockMesh &> blockMesh.log
# decomposePar &> decomposePar.log
# mpirun -np 2 EulerSolidMechFoam -parallel &> EulerSolidMechFoam.log
Euler_Solid_MechFoam &> Euler_Solid_MechFoam.log