
# load modules
env2lmod
module load gcc/6.3.0 openmpi/3.0.1 fenics/2019.1.0

# run job with email message, but without parameters
bsub -N -B -W 24:00 -n 48 -R fullnode 
