#!/bin/bash -l
#SBATCH --job-name=verve
#SBATCH --output=outputr%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=qib-compute
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=100G

bash scripts/erv_locus_check.sh CMEH_Vag HERVH_20p11.23b
bash scripts/erv_locus_check.sh CM28_Sa3 MER34B_1q23.3a
