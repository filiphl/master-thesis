class Single:

    def __init__(self, path):
        self.path = path
        infile = open(path, 'r')
        self.file = infile.read()
        infile.close()

        self.data , self.key = self.retrieveData(path)

        self.F = np.asarray(self.data['f_smdID[1]'])
        self.N = np.mean(self.data['f_addforce[3]'])
        self.v = float(re.findall(r'vel=-([^\s]*)', self.file)[0])
        self.t = np.asarray(self.data['Step'])

        self.dt = 0.002 #ps

        self.Fmax = self.findFmax(self.F)

    def retrieveData(self, path):
        infile = open(path, 'r')
        slurmData = {}
        key = []
        # -------------------------- Fill key -------------------------------- #
        for line in infile:
            col = line.split()
            if str(col[0]) == 'Step':
                nCols = len(col)
                for i in xrange(nCols):
                    slurmData[str(col[i])] = []
                    key.append(str(col[i]))
                break
        # ------------------------ Fill slurmData ---------------------------- #
        for line in infile:
            col = line.split()
            try:
                for i in xrange(nCols):
                    slurmData[key[i]].append(float(col[i]))
            except:
                for j in xrange(i):
                    slurmData[key[j]] = slurmData[key[j]][:-1]
                continue
        infile.close()

        return slurmData, key



    def findFmax(self, F):
        for i in xrange(11,len(F)-11):
            c=0
            for j in xrange(10):
                if F[i-j]< F[i] > F[i+j]:
                    c+=1
            if c>8:
                break   # Exit at first top
        return F[i]



    def printObject(self):
        print '{0:<10}:  {1:<30}'.format('File', self.path)
        print '{0:<10}:  {1:<30}'.format('Entries', len(self.F))
        print '{0:<10}:  {1:<30}'.format('Fmax', self.Fmax)
        print '{0:<10}:  {1:<30}'.format('N', self.N)
        print '{0:<10}:  {1:<30}'.format('loadStep', self.t)
        print '{0:<10}:  {1:<30}'.format('v', self.v)



class SameTime:
    def __init__(self, data):
        self.data = data
        self.F      = [i[0] for i in data]
        self.T      = [i[1] for i in data]
        self.T      = [T - T[0] for T in self.T]
        self.N      = np.mean([i[2] for i in data])
        self.Fmax   = [i[3] for i in data]
        self.v      = [i[4] for i in data]
        self.dt     = 0.002 #ps

        self.c  = ['#D7DFC0','#AECFA2','#82BC92','#5FA38E','#49848B','#3F5F7F','#383D65','#2C1E3E']
        self.cc = ['#AF97A9','#E5B3BB','#FCB7AE','#F7B490','#E3AA75','#C49859','#937E47','#5B5728']

        self.u = self.findUofV(1)
        self.scale = lambda F, v: F/(self.u[0]*v + self.u[1])   # u(v)


    def findUofV(self, plot=False):
        pf = np.polyfit(self.v,self.Fmax/self.N,1)
        l = np.linspace(0, max(self.v), 100)
        p = np.polyval(pf, l)
        if plot:
            plt.figure()
            print pf
            plt.plot(l,p,
            '--',
            linewidth=3,
            color='#666666',
            label='$\\mu(v) = %.03fv + %.03f$'%(pf[0], pf[1]))

            plt.hold('on')

            plt.plot(self.v,self.Fmax/self.N,
            'o',
            markersize=10,
            color=self.c[6])

            l= plt.legend(loc='upper left')
            l.legendHandles[0].set_linewidth(3)

            plt.title('$N=%.02f$'%self.N)



        return pf

    def plotFNT(self):
        T = self.T
        dt = self.dt
        v = self.v
        F = self.F
        N = self.N
        scale = self.scale
        d = [T[i]*dt*v[i] for i in xrange(len(T))]
        plt.figure()
        for i in xrange(len(F)):
            plt.plot(scale(d[i], v[i]), scale(F[i], v[i])/N,
            '-',
            linewidth=3,
            color=self.c[i+2],
            label='$v=%d$'%int(v[i]))
            plt.hold('on')
        plt.legend()
        plt.ylabel('$\\frac{F}{N\\mu(v)}$', rotation=0, fontsize=28)#'Force [eV/\AA]', fontsize=16)
        plt.xlabel('$\\frac{t v}{\\mu(v)}$', fontsize=28)
        plt.tick_params(axis='both', which='major', labelsize=16)
        ax = plt.gca()
        ax.yaxis.set_label_coords(-0.12, 0.455)
        plt.title('$N=%.02f$'%self.N)


################################################################################
if __name__ == '__main__':
    import os, re
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    data = {}
    time = []
    for path in sorted(os.listdir('.')):
        if path[-3:]=='out':
            s = Single(path)
            try:
                data[s.t[0]].append([s.F, s.t, s.N, s.Fmax, s.v])
            except:
                data[s.t[0]] = [[s.F, s.t, s.N, s.Fmax, s.v]]
                time.append(s.t[0])

    time = sorted(time)

mu = []
N = []

for t in time:
    st = SameTime(data[t])
    #st.plotFNT()
    #plt.show()
    mu.append(st.u[-1])
    N.append(st.N)

F = np.asarray(mu)*np.asarray(N)

print '{0:<7}  {1:<7}  {2:<7}'.format('t', 'N(t)', 'mu(t)')
for t in xrange(len(time)):
    print '{0:<7}  {1:<7.02f}  {2:<7.02f}'.format(t, N[t], mu[t])

pf = np.polyfit(N, F, 1)
l = np.linspace(min(N), max(N), 100)
p = np.polyval(pf, l)

plt.figure()

plt.plot(l,p,
'--',
linewidth=3,
color='#666666',
label="$F\'(N) = %.03fN + %.03f$"%(pf[0], pf[1]))

plt.hold('on')

plt.plot(N, F,
'o',
markersize=10,
color=st.cc[6])

l= plt.legend(loc='upper left', fontsize=16)
l.legendHandles[0].set_linewidth(3)


#plt.ylabel('$F$~[ev/\\AA]', rotation=90, fontsize=24)#'Force [eV/\AA]', fontsize=16)
plt.ylabel('$F$~[ev/\\AA]', rotation=90, fontsize=24)
plt.xlabel('$N$~[ev/\\AA]', fontsize=24)
plt.tick_params(axis='both', which='major', labelsize=16)
ax = plt.gca()
#ax.yaxis.set_label_coords(-0.12, 0.455)
plt.tight_layout()


plt.figure()
plt.plot(N, mu,
'o',
markersize=10,
color=st.cc[6])


#plt.ylabel('$F$~[ev/\\AA]', rotation=90, fontsize=24)#'Force [eV/\AA]', fontsize=16)
plt.ylabel('$\\mu', rotation=0, fontsize=28)
plt.xlabel('$N$~[ev/\\AA]', fontsize=24)
plt.tick_params(axis='both', which='major', labelsize=16)
ax = plt.gca()
ax.yaxis.set_label_coords(-0.12, 0.48)
plt.tight_layout()




plt.show()
