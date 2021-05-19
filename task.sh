#!/bin/sh
#PBS -N million_256_2
#PBS -e errores.out
#PBS -o output.out
#PBS -M german.chaparro@udea.edu.co
#PBS -m b
#PBS -m e
#PBS -q long
#PBS -l nodes=2:ppn=32
cd $PBS_O_WORKDIR
python3 parallel_sketch.py no_pert_nb million_256_2
