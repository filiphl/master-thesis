import os
from sys import argv
from time import sleep


sleep(5)
SLURM_JOB_ID = int(argv[1])
slurmFile = 'slurm-%d.out'%SLURM_JOB_ID

c = 0
while slurmFile not in os.listdir('.'):
	sleep(1)
	c+=1
	if c == 1000:
		break

os.system('scp slurm-%d.out kontor:Other/random/slurm.out'%SLURM_JOB_ID)
command = "python Other/random/sendMail.py 'Job %d done!' 'Other/random/slurm.out'; rm Other/random/slurm.out; exit"%SLURM_JOB_ID
sleep(10)
os.system('ssh kontor "%s"'%command)


