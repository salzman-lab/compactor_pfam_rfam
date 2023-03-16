#!/bin/sh
#SBATCH -p horence,normal,quake
#SBATCH --time=6:00:00
#SBATCH --mem=32000      # In MB
#SBATCH --job-name=prfam     # job name

source /oak/stanford/groups/horence/george/dog/bin/activate

run_name=" "
path_to_compactors=" "

python3 bin/rfam_pfam.py "$run_name" "$path_to_compactors"
