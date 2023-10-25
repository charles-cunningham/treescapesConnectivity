#!/bin/bash
#SBATCH --job-name=ST_SDMs                    			# Job name
#SBATCH --time=09:55:00                           		# Time limit hrs:min:sec 167:55:00
#SBATCH --ntasks=1                        			# How many tasks
#SBATCH --cpus-per-task=4                 			# ...with two cores
#SBATCH --mem=99950                         			# Minimum memory required per allocated CPU
#SBATCH --output Treescapes/ST_SDMs/run_logs/run_o%j.log 	# Standard output and error log
#SBATCH --error Treescapes/ST_SDMs/run_logs/run_e%j.log      	# Standard output and error log
#SBATCH --account=biol-distmdl-2018
#SBATCH --array=11-12   #1-65

cd /users/cac567/scratch

module load lang/R/4.2.1-foss-2022a

Rscript Treescapes/ST_SDMs/Scripts/09_INLABRU_ST_SDM_cluster.R ${SLURM_ARRAY_TASK_ID}
