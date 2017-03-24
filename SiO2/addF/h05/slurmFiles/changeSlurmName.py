import re, os

def renameByLoadTime(filepath):
    content = open(filepath, 'rw').read()
    vel = re.findall(r'vel=-([^\s]*)', content)
    loadStep = re.findall(r'Current step  : ([^s\n]*)', content)
    if (len(vel)>0 and len(loadStep)>0):
        return vel[0], loadStep[0]
    elif (len(vel)==0 and len(loadStep)>0):
        return 0, loadStep[0]
    else:
        return 'error'

if __name__ == '__main__':
    folderPath = '.'
    directory = os.listdir(folderPath)
    for filepath in directory:
        if 'slurm' in filepath:
            matches = renameByLoadTime(filepath)
            print filepath,  matches
            if matches=='error':
                os.system('rm %s'%filepath)
            newName = 't%dv%d.out'%(int(matches[1]), int(matches[0]))
            os.system('mv %s %s'%(filepath, newName))
