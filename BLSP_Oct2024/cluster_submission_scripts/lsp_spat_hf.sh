#!/bin/bash

#BSUB -J burn_80_80                 # job name AND job array
#BSUB -n 20                         # number of cores: generally only need 1 to 5
#BSUB -R span[hosts=1]              # restricts new cores to one node
##BSUB -x
#BSUB -R "rusage[mem=25GB]"
##BSUB -R select[model=Plat8358]
#BSUB -W 15:00                      # walltime limit: hours:minutes
#BSUB -o ./outs/out_burn_80_80.out  # filepath for the output
#BSUB -e ./err/err_burn_80_80.err   # filepath for error
source ~/.bashrc
conda activate /usr/local/usrapps/bjreich/mpshisle/lsp_spat_sim_env
Rscript lsp_spat_hf_burn_mcmc.R 80 80 2 1000 1 20 50
conda deactivate
