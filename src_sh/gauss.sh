#!/bin/bash
#SBATCH --job-name DPPD
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#SBATCH --mem-per-cpu=24G
#SBATCH --time=8:00:00
#SBATCH --account=synthesis
#SBATCH --error=std.err_%j
#SBATCH --output=std.out_%j
#SBATCH --partition shared

# Load Gaussian module to set environment
module load gaussian python


cd $SLURM_SUBMIT_DIR
echo $PWD
INPUT_BASENAME=optim_nodiffuse_DPPD
GAUSSIAN_EXEC=g16

# Run gaussian NREL script (performs much of the Gaussian setup)
g16_nrel

#Setup Linda parameters
if [ $SLURM_JOB_NUM_NODES -gt 1 ]; then 
export GAUSS_LFLAGS='-vv -opt "Tsnet.Node.lindarsharg: ssh"' 
export GAUSS_EXEDIR=$g16root/g16/linda-exe:$GAUSS_EXEDIR 
fi 

# Run Gaussian job 
$GAUSSIAN_EXEC < $INPUT_BASENAME.com >& $INPUT_BASENAME.log 

echo "Completed successfully .. :)"
