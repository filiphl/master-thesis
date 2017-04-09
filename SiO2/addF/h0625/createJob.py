import os, sys, re

inJob           = open('jobs/job.sh', 'r')
inInputScript   = open('inputScripts/system.run', 'r')

try:
    step = int(sys.argv[1])
except:
    print 'First argument represents step, and should be an int.'

try:
    velocity = float(sys.argv[2])
except:
    print 'Second argument represents velocity, and should be an int.'

try:
    N = int(sys.argv[3])
except:
    print 'Third argument represents # time steps, and should be an int.'


job             = 'jobs/job%dv%d.sh'%(step,velocity*100)
inputScript     = 'inputScripts/system.run.%dv%d'%(step,velocity*100)
outJob          = open(job, 'w')
outInputScript  = open(inputScript, 'w')

inputContent    = inInputScript.read()
jobContent      = inJob.read()
inJob.close()
inInputScript.close()

jobContent = re.sub(r'(-in\s)(.*)', r'\1%s'%inputScript, jobContent)
jobContent = re.sub(r'job-name=filip', r'job-name=filip%d-%d'%(step/1000,velocity), jobContent)

inputContent = re.sub(r'(variable\s*vel\s*equal\s*)(-\d*)', r'\1 -%.02f'%velocity,   inputContent)
inputContent = re.sub(r'(variable\s*loadStep\s*equal\s*[^\d])(\d*)', r'\1 %d'%step, inputContent)
inputContent = re.sub(r'(variable\s*N\s*equal\s*[^\d])(\d*)', r'\1 %d'%N, inputContent)

outJob.write(jobContent)
outInputScript.write(inputContent)

outJob.close()
outInputScript.close()

sbatch = 'sbatch %s'%job
os.system(sbatch)
