from matplotlib import pyplot as plt
import numpy as np
from matplotlib import rc




rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)
lc = ["#5FA38E", "#49848B", "#3F5F7F"]




def decay(x,c,a):
    return (2/3.0)*np.exp(c*(x[0]-x)) + a*( 1-np.exp(c*(x[0]-x)) )


def stick(x,a,b):
    return a*(x-x[0])+b

N  = 2000
N2 = N/2
x = np.linspace(0,2,N+1)
X=[]
Y = []
limit = [i-N/6 for i in range(N2,N, N/6)]

s = []
for i in xrange(len(limit)-1):
    s.append([limit[i]+50,limit[i+1]])




fig,ax = plt.subplots(1,1, figsize=(10,5))
y = stick(x[:N2-N/6], 1, 0)
ax.plot(x[:N2-N/6], y,
'-',
color=lc[0],
linewidth=4)
X.append(x[:N2-N/6])
Y.append(y)
plt.hold('on')

y = decay(x[N2-N/6:s[0][0]], 120, 0.384)
ax.plot(x[N2-N/6:s[0][0]], y,
'-',
color=lc[0],
linewidth=4)
X.append(x[N2-N/6:s[0][0]])
Y.append(y)

for i in xrange(len(s)):
    y = stick(x[s[i][0]:s[i][1]], 1, 0.384)
    ax.plot(x[s[i][0]:s[i][1]], y,
    '-',
    color=lc[0],
    linewidth=4)
    X.append(x[s[i][0]:s[i][1]])
    Y.append(y)
    try:
        y = decay(x[s[i][1]:s[i+1][0]], 120, 0.384)
        ax.plot(x[s[i][1]:s[i+1][0]], y,
        '-',
        color=lc[0],
        linewidth=4)
        X.append(x[s[i][1]:s[i+1][0]])
        Y.append(y)
    except:
        continue



#y2 = decay(x[N2:], 50, 0.7)








#ax.legend(loc=1)
ax.grid()
ax.set_ylabel(r"Friction force", fontsize=20)
ax.set_xlabel(r"Time"  , fontsize=20)
ax.xaxis.set_label_coords(0.5, -0.05)
ax.yaxis.set_label_coords(-0.025, 0.5)

ax.set_ylim([0,max(y)*1.1])











#------------------------------------------------------------------------------#

ax.set_xticks([])
ax.set_yticks([])

#ax.set_yticks([0.7, 1.0])
#ax.set_yticklabels([r'$F_k$', r'$F_s$'])



xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()

# removing the default axis on all sides:
for side in ['bottom','right','top','left']:
    ax.spines[side].set_visible(False)


# get width and height of axes object to compute
# matching arrowhead length and width
dps = fig.dpi_scale_trans.inverted()
bbox = ax.get_window_extent().transformed(dps)
width, height = bbox.width, bbox.height

# manual arrowhead width and length
hw = 1./30.*(ymax-ymin)
hl = 1./50.*(xmax-xmin)
lw = 2.5 # axis line width
ohg = 0.384 # arrow overhang

# compute matching arrowhead length and width
yhw = hw/(ymax-ymin)*(xmax-xmin)* height/width
yhl = hl/(xmax-xmin)*(ymax-ymin)* width/height

# draw x and y axis
ax.arrow(xmin, 0, xmax-xmin, 0., fc='#333333', ec='#333333', lw = lw,
         head_width=hw, head_length=hl, overhang = ohg,
         length_includes_head= True, clip_on = False)

ax.arrow(0, ymin, 0., ymax-ymin, fc='#333333', ec='#333333', lw = lw,
         head_width=yhw, head_length=yhl, overhang = ohg,
         length_includes_head= True, clip_on = False)




plt.show()
#plt.savefig('../../thesis/figures/friction/steadySlide.pdf')

Y2 = []
for y in Y:
    for i in y:
        Y2.append(i)

Y = Y2
X = x[:len(Y)]

# SMOOTH
def smoother(degree, function):
    smoothed = np.zeros(len(function))
    smoothed[0] = function[0]

    for i in xrange(1,degree+1):
        smoothed[i] = smooth(i,function, i)
    for i in xrange(degree, len(function)-degree):
        smoothed[i] = smooth(degree,function,i)
    for i in xrange(2,degree+2):
        smoothed[-i] = smooth(i-1,function, -i)
    smoothed[-1] = function[-1]
    return smoothed

def smooth(degree, function, atIndex):
    value            = 0.0
    dividor          = 0.0
    localCoeffisient = 1.0
    for i in xrange(-degree, degree+1):
        dividor += localCoeffisient
        value += localCoeffisient*function[atIndex+i]
        if i < 0:
            localCoeffisient += 1
        if i == 0:
            localCoeffisient -= 1
        if i > 0:
            localCoeffisient -= 1

    localCoeffisient +=1
    return value/dividor

Y = smoother(10,Y)

X = X[:-10]
Y = Y[:-10]

fig,ax = plt.subplots(1,1, figsize=(10,5))
ax.plot(X,Y,
'-',
color=lc[1],
linewidth=4)


#ax.legend(loc=1)
ax.grid('on')
ax.set_ylabel(r"Friction force", fontsize=20)
ax.set_xlabel(r"Time"  , fontsize=20)
ax.xaxis.set_label_coords(0.5, -0.05)
ax.yaxis.set_label_coords(-0.025, 0.5)

ax.set_ylim([0,max(y)*1.1])


ax.set_xticks([])
ax.set_yticks([])

#ax.set_yticks([0.7, 1.0])
#ax.set_yticklabels([r'$F_k$', r'$F_s$'])



xmin, xmax = ax.get_xlim()
ymin, ymax = ax.get_ylim()

# removing the default axis on all sides:
for side in ['bottom','right','top','left']:
    ax.spines[side].set_visible(False)


# get width and height of axes object to compute
# matching arrowhead length and width
dps = fig.dpi_scale_trans.inverted()
bbox = ax.get_window_extent().transformed(dps)
width, height = bbox.width, bbox.height

# manual arrowhead width and length
hw = 1./30.*(ymax-ymin)
hl = 1./50.*(xmax-xmin)
lw = 2.5 # axis line width
ohg = 0.384 # arrow overhang

# compute matching arrowhead length and width
yhw = hw/(ymax-ymin)*(xmax-xmin)* height/width
yhl = hl/(xmax-xmin)*(ymax-ymin)* width/height

# draw x and y axis
ax.arrow(xmin, 0, xmax-xmin, 0., fc='#333333', ec='#333333', lw = lw,
         head_width=hw, head_length=hl, overhang = ohg,
         length_includes_head= True, clip_on = False)

ax.arrow(0, ymin, 0., ymax-ymin, fc='#333333', ec='#333333', lw = lw,
         head_width=yhw, head_length=yhl, overhang = ohg,
         length_includes_head= True, clip_on = False)



plt.savefig('../../thesis/figures/friction/stick-slip.pdf')
plt.show()
