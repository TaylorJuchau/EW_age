#!/bin/bash
#SBATCH --account=galaxies
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tjuchau@uwyo.edu
            
#SBATCH --job-name=generating_small_slug_library
#SBATCH --output=generating_small_slug_library.out
#SBATCH --error=generating_small_slug_library.err
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=30
#SBATCH --nodes=1

source ~/.bashrc
conda activate Modeling
cd /cluster/medbow/project/galaxies/tjuchau/projects/EW_vs_Age/modeling_files
conda run -n Modeling /project/galaxies/tjuchau/software/Slug/slug2/bin/slug /project/galaxies/tjuchau/projects/EW_vs_Age/modeling_files/make_library.slugin