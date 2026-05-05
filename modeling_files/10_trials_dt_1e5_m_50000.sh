#!/bin/bash
#SBATCH --account=galaxies
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tjuchau@uwyo.edu
            
#SBATCH --job-name=10_trials_dt_1e5_m_50000
#SBATCH --output=10_trials_dt_1e5_m_50000.out
#SBATCH --error=10_trials_dt_1e5_m_50000.err
#SBATCH --time=8:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=5
#SBATCH --nodes=1

source ~/.bashrc
conda activate Modeling
cd /cluster/medbow/project/galaxies/tjuchau/projects/EW_vs_Age/modeling_files
conda run -n Modeling /project/galaxies/tjuchau/software/Slug/slug2/bin/slug /project/galaxies/tjuchau/projects/EW_vs_Age/modeling_files/10_trials_dt_1e5_m_50000.slugin