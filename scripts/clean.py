# Completely clean out the run folder and return it to an 'empty' state
import subprocess

run_dir ='../'
w_dir = '../workspace'
namelist_dir = '../namelists/'
log_dir = '../logs'
dom_dir = '../DOMAIN'
forc_dir = '../FORCING'
try :
	subprocess.run('rm -r '+log_dir,shell=True)
except:
	print('No such Files in '+log_dir)
try:
	subprocess.run('rm '+namelist_dir+'*',shell=True)
except:
	print('No such Files in '+namelist_dir)
try:
	subprocess.run('rm -r '+workspace,shell=True)
except:
	print('ERROR: cannot remove workspace or does not exist')
try: 
	subprocess.run('rm -r -i '+dom_dir,shell=True)
except: 
	print('ERROR: cannot remove DOMAIN or does not exist')
try:
	subprocess.run('rm -r -i '+forc_dir,shell=True)
except:
	print('ERROR: cannot remove FORCING or does not exist')

