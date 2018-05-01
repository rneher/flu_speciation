#!/bin/sh
#
# Reserve 1 CPUs for this job
#
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#
# Request it to run this for HH:MM:SS with ?G per core
#
#SBATCH --time=00:59:00
#
#  Run job from current working directory
#$ -cwd
#

source ~/.bashrc
python $@
