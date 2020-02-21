import subprocess
import sys

argv = sys.argv[1:]
n = argv[0]

subprocess.run('cp tables/* .',shell=True)
subprocess.run('cp namelists/namelist.hrldas .',shell=True)
subprocess.run('cp namelists/hydro.namelist .',shell=True)
subprocess.run('mpiexec -n '+n+' wrf_hydro.exe',shell=True)
subprocess.run('rm namelist.hrldas',shell=True)
subprocess.run('rm hydro.namelist',shell=True)
subprocess.run('rm *.TBL',shell=True)
subprocess.run('mv diag* OUT/',shell=True)

