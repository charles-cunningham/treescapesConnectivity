#!/bin/bash
#SBATCH --job-name=setup_omniscape                    		# Job name
#SBATCH --time=00:20:00                           		# Time limit hrs:min:sec 167:55:00
#SBATCH --ntasks=1	                        		# How many tasks on each node
#SBATCH --mem=199950                        			# Minimum memory required per allocated CPU
#SBATCH --output Treescapes/Omniscape/run_logs/run_o%j.log     	# Standard output and error log
#SBATCH --error Treescapes/Omniscape/run_logs/run_e%j.log      	# Standard output and error log
#SBATCH --account=biol-distmdl-2018

cd /users/cac567/scratch

module load lang/R/4.2.1-foss-2022a

Rscript Treescapes/Omniscape/06_Omniscape_setup.R ${SLURM_ARRAY_TASK_ID}
