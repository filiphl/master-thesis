from os import listdir, system

directory=sorted(listdir('.'))
ts = []
for f in directory:
    try:
        t = int(f[10:-4])
        ts.append(t)
        #system('mv %s %s'%(f,newName))
    except:
        continue
for t in ts:
    newName = 'hysteresis%d.png'%(t-10)
    cmd = 'mv hysteresis%d.png %s'%(t,newName)
    system(cmd)
