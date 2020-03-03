# Completely clean out the run folder and return it to an 'empty' state
import subprocess

run_dir ='/home/tsw35/soteria/wrf_hydro/wrf_hydro_run'
w_dir = run_dir+'/workspace/'
namelist_dir = run_dir+'/namelists/'
log_dir = run_dir+'/logs/'
dom_dir = run_dir+'/DOMAIN/'
forc_dir = run_dir+'/FORCING/'
subprocess.run('rm -r '+log_dir,shell=True)
subprocess.run('rm '+namelist_dir+'*',shell=True)
subprocess.run('rm -r '+w_dir,shell=True)
subprocess.run('rm -r '+dom_dir,shell=True)
subprocess.run('rm -r '+forc_dir,shell=True)
subprocess.run('rm -r '+run_dir+'/OUT/',shell=True)
subprocess.run('rm -r '+run_dir+'/hydro2dtbl*',shell=True)
subprocess.run('rm '+run_dir+'/*.nc',shell=True)
