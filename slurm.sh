#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -A phys5630
#SBATCH --output=slurm.log
#SBATCH --partition=instructional
##SBATCH --partition=standard  # default

# Simple "executable script for testing slurm batch system"
# Normally you would use a script like this to do the usual things
# you need to do when working with programs:
# o) cd to a working directory
# o) setup links, variables, copy data files, etc
# o) run a program with any necessary command line parameters
# o) copy or process results in another program
# o) ...

# For PHYS5630, we want to setup the same environment that is present
# for interactive logins

source ~/.bashrc
conda activate phys56xx

# Below we can run any program(s), eg:
# cd ~/examples/fitsys
# ./fitsys
# cd ~/homework/really_long_problem
# ./really_long_problem input.dat
# ...

# By default SLURM changes to the directory from which the job was
# submitted, so the SLURM_SUBMIT_DIR environment variable is usually not needed.
#cd $SLURM_SUBMIT_DIR

#NMILLION_PTS=10000

#./Integral5d -b $NMILLION_PTS  
#./simpix Claude_Monet_The_Cliffs_at_Etretat.png Claude_Monet,_Le_Grand_Canal.png out_slurm.png
./simpix Claude_Monet,_Le_Grand_Canal.png Claude_Monet_The_Cliffs_at_Etretat.png out_slurm2.png

echo "goodbye"
