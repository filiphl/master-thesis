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
for t in ts[::-1]:
    newName = 'hysteresis%d.png'%(t/3+10)
    if t < 10:
        cmd = 'mv hysteresis000%d.png %s'%(t,newName)
    elif t<100:
        cmd = 'mv hysteresis00%d.png %s'%(t,newName)
    else:
        cmd = 'mv hysteresis0%d.png %s'%(t,newName)
    system(cmd)
