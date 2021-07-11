#!/bin/bash
#SBATCH --job-name=kBET
#SBATCH -c 4
#SBATCH --time=08:10:00
##SBATCH --mem=48G
#SBATCH --mem-per-cpu=25G
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=all
# send email on start, end and fault
#SBATCH --mail-user=jdeng@umcutrecht.nl
cd /hpc/dla_lti/jdeng/scRNASeq/PBMC
alias loadR4013.10="export PATH=/hpc/local/CentOS7/dla_lti/R-4.0.1/bin:$PATH && export R_LIBS_USER=/hpc/dla_lti/jdeng/.Rlib/R-4.0.1-Bioc-3.10 && module load rstudio/1.0.136"
loadR4013.10
Rscript kBET.R
