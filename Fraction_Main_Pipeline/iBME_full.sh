#!/bin/sh

#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --job-name=iBME_full
#SBATCH --time=12:00:00
#SBATCH --mem=8G
##SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=/users/t/j/tjaglal/output/ibme_out/logs_%j.out
#SBATCH --error=/users/t/j/tjaglal/output/ibme_out/logs_%j.err
#SBATCH --mail-type=BEGIN,END

cd ~/Projects/iBME
source iBME_env/bin/activate
python frame_fraction.py
deactivate


