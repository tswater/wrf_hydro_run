#SBATCH --ntasks-per-node=32
#SBATCH --nodes=1
#SBATCH --job-name="hydro-1"
#SBATCH --output="log.txt"
cd ../
mpiexec -n 32 wrf_hydro.exe
