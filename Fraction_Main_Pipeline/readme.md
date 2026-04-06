## Main pipeline for running iBME reweighting on structure fractions
### More efficient compute time

Main file is iBME_full.sh. Run this file in slurm to start process
***

frame_fraction.py - main working script. Run sbatch for entire process. Define grid here

iBME_grid_scan.py - simulate structure fractions through Pepsi and concatenate them by grid point

iBME_reweight.py- take combined grid point files and perform iBME reweighting and heat map plotting