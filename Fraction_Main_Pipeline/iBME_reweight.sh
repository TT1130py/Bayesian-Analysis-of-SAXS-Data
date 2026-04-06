#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --job-name=iBME_scan
#SBATCH --time=24:00:00
#SBATCH --mem=8G
##SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --output=/users/t/j/tjaglal/output/ibme_out/logs_%j.out
#SBATCH --error=/users/t/j/tjaglal/output/ibme_out/logs_%j.err
#SBATCH --mail-type=BEGIN,END

THETA=$1
SAVE_PATH=$2

cd ~/Projects/iBME
source iBME_env/bin/activate

python iBME_reweight.py $THETA $SAVE_PATH

deactivate
