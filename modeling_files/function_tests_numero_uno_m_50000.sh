#!/bin/bash
#SBATCH --account=galaxies
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tjuchau@uwyo.edu
            
#SBATCH --job-name=function_tests_numero_uno_m_50000
#SBATCH --output=function_tests_numero_uno_m_50000.out
#SBATCH --error=function_tests_numero_uno_m_50000.err
#SBATCH --time=8:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=5
#SBATCH --nodes=1

source ~/.bashrc
conda activate Modeling
cd /cluster/medbow/project/galaxies/tjuchau/projects/EW_vs_Age/modeling_files
conda run -n Modeling /project/galaxies/tjuchau/software/Slug/slug2/bin/slug /project/galaxies/tjuchau/projects/EW_vs_Age/modeling_files/function_tests_numero_uno_m_50000.slugin