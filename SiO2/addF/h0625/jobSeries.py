import os


#v = [0.5, 1, 2, 3, 4]
v = [1,2,3,4,5]

#N = [80000, 50000, 40000, 30000, 20000] 
#N = [60000, 50000, 40000, 30000, 20000]
N = [80000, 60000, 50000, 40000, 30000]

#steps = [50000, 60000, 70000, 80000, 90000]
#steps = [100000, 110000, 120000, 130000, 140000]
steps = [150000, 160000, 170000, 180000, 190000, 200000]

for step in steps:
    for i in range(len(v)):
        os.system('python createJob.py %d %.02f %d'%(step, v[i], N[i]))