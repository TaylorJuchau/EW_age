#!/bin/bash
#SBATCH --account=galaxies
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tjuchau@uwyo.edu
            
#SBATCH --job-name=Full_functionality_test_dt5_m_5000
#SBATCH --output=Full_functionality_test_dt5_m_5000.out
#SBATCH --error=Full_functionality_test_dt5_m_5000.err
#SBATCH --time=8:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

source ~/.bashrc
conda activate Modeling
cd /cluster/medbow/project/galaxies/tjuchau/projects/EW_vs_Age/modeling_files
conda run -n Modeling /project/galaxies/tjuchau/software/Slug/slug2/bin/slug /project/galaxies/tjuchau/projects/EW_vs_Age/modeling_files/Full_functionality_test_dt5_m_5000.slugin