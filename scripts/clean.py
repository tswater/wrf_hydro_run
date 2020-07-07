# Completely clean out the run folder and return it to an 'empty' state
import subprocess
import os

run_dir = os.getcwd()
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
subprocess.run('rm -r '+run_dir+'/*.TBL',shell=True)
subprocess.run('rm '+run_dir+'/hydro.namelist',shell=True)
subprocess.run('rm '+run_dir+'/namelist.hrldas',shell=True)
subprocess.run('rm '+run_dir+'/diag*',shell=True)
subprocess.run('rm '+run_dir+'/*.nc',shell=True)
