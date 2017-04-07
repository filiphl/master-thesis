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
    def __init__(self, data, plot=False):
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

        self.u = self.findUofV(plot)
        self.scale = lambda F, v: F/(self.u[0]*v + self.u[1])   # u(v)


    def findUofV(self, plot=False):
        pf = np.polyfit(self.v,self.Fmax/self.N,1)
        l = np.linspace(0, max(self.v), 100)
        p = np.polyval(pf, l)
        if plot:
            plt.figure()
            #print pf
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



class All:

    def __init__(self, dirPath='.', plotEach=False):

        self.majorFontSize = 28
        self.minorFontSize = 25
        self.tickFontSize  = 23

        self.c  = ['#D7DFC0','#AECFA2','#82BC92','#5FA38E','#49848B','#3F5F7F','#383D65','#2C1E3E']
        self.cc = ['#AF97A9','#E5B3BB','#FCB7AE','#F7B490','#E3AA75','#C49859','#937E47','#5B5728']

        self.data = {}
        self.time = []
        self.mu   = []
        self.N    = []
        for path in sorted(os.listdir(dirPath)):
            if path[-3:]=='out':
                s = Single(path)
                try:
                    self.data[s.t[0]].append([s.F, s.t, s.N, s.Fmax, s.v])
                except:
                    self.data[s.t[0]] = [[s.F, s.t, s.N, s.Fmax, s.v]]
                    self.time.append(s.t[0])

        self.time = sorted(self.time)

        self.findCoefficient(plotEach=plotEach)



    def findCoefficient(self, plotFinal=False, plotEach=False, printTable=False):
        for t in self.time:
            st = SameTime(self.data[t], plot=plotEach)
            #st.plotFNT()
            #plt.show()
            self.mu.append(st.u[-1])
            self.N.append(st.N)

        self.F = np.asarray(self.mu)*np.asarray(self.N)

        if printTable:
            print '{0:<2}    {1:<7}    {2:<7}'.format('t', 'N(t)', 'mu(t)')
            for t in xrange(len(time)):
                print '{0:<2}    {1:<7.02f}    {2:<7.03f}'.format(t, N[t], mu[t])

        self.pf = np.polyfit(self.N, self.F, 1)
        self.l = np.linspace(min(self.N), max(self.N), 100)
        self.p = np.polyval(self.pf, self.l)

        if plotFinal:
            self.plotAll()



    def plotAll(self):
        plt.figure()

        plt.plot(self.l, self.p,
        '--',
        linewidth=3,
        color='#666666',
        label="$F\'(N) = %.03fN + %.03f$"%(self.pf[0], self.pf[1]))

        plt.hold('on')

        plt.plot(self.N, self.F,
        'o',
        markersize=10,
        color=self.cc[6])


        lgnd = plt.legend(loc='upper left', fontsize=self.tickFontSize)
        lgnd.legendHandles[0].set_linewidth(3)
        #plt.ylabel('$F$~[ev/\\AA]', rotation=90, fontsize=24)#'Force [eV/\AA]', fontsize=16)
        plt.ylabel('$F$~[ev/\\AA]', rotation=90, fontsize=self.majorFontSize)
        plt.xlabel('$N$~[ev/\\AA]', fontsize=self.majorFontSize)
        plt.xticks(np.arange(400 , 2201 , 400))
        plt.yticks(np.arange(300 , 751 , 100))
        plt.tick_params(axis='both', which='major', labelsize=self.tickFontSize)
        #ax = plt.gca()
        #ax.yaxis.set_label_coords(-0.12, 0.455)
        plt.tight_layout()

        """
        plt.figure()
        plt.plot(self.N, self.mu,
        'o',
        markersize=10,
        color=self.cc[6])
        #plt.ylabel('$F$~[ev/\\AA]', rotation=90, fontsize=24)#'Force [eV/\AA]', fontsize=16)
        plt.ylabel('$\\mu', rotation=0, fontsize=28)
        plt.xlabel('$N$~[ev/\\AA]', fontsize=24)
        plt.tick_params(axis='both', which='major', labelsize=16)
        ax = plt.gca()
        ax.yaxis.set_label_coords(-0.12, 0.48)
        plt.tight_layout()
        """
    def residualPlot(self):
        self.residual = []

        for i in xrange(len(self.F)):
            self.residual.append(self.F[i]-np.polyval(self.pf, self.N[i]))
            #print self.N[i], self.F[i], np.polyval(self.pf, self.N[i]), self.residual[i])

        plt.figure()
        plt.stem(self.N, self.residual)
        plt.ylabel('Residual~[ev/\\AA]', rotation=90, fontsize=24)
        plt.xlabel('$N$~[ev/\\AA]', fontsize=24)
        plt.axis([min(self.N)*0.9, max(self.N)+min(self.N)*0.1, min(self.residual)*1.1, max(self.residual)*1.1 ])
        plt.tick_params(axis='both', which='major', labelsize=16)


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


    h1 = All(plotEach=0)
    h1.plotAll()
    #h1.residualPlot()
    plt.show()
