from matplotlib import pyplot as plt
import numpy as np
from matplotlib import rc





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



def decay(x,c,a):
    return np.exp(c*(x[0]-x)) + a*( 1-np.exp(c*(x[0]-x)) )


N  = 2000
N2 = N/2
x = np.linspace(0,2,N+1)
y = np.linspace(0,2,N+1)
y[N2:] = decay(x[N2:], 50, 0.7)

y=smoother(7,y)

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)
lc = ["#5FA38E", "#3F5F7F"]
LC = [(95/255.,163/255.,142/255.), (63/255.,95/255.,127/255.)]

fig,ax = plt.subplots(1,1, figsize=(10,5))


l1 = [max(y) , max(y)]
l2 = [y[-1], y[-1]]


ax.plot([x[N2/4], x[N2]], l1,
'--',
color='#666666',
linewidth=3)

ax.plot([x[N2/4], x[N2+N2/8]], l2,
'--',
color='#666666',
linewidth=3)

less = 15
more = 3
d = float(less+more)

ax.plot(x[3:N2-less], y[3:N2-less],
'-',
color=lc[0],
linewidth=4,
label='Loading')


plt.hold('on')

for i in xrange(less+more):
    ax.plot(x[N2-less+i:N2-less+i+2], y[N2-less+i:N2-less+i+2],
    '-',
    color=( (i/d) * (LC[1][0]-LC[0][0]) + LC[0][0], (i/d) * (LC[1][1]-LC[0][1]) + LC[0][1], (i/d) * (LC[1][2]-LC[0][2]) + LC[0][2]  ),
    linewidth=4,
    )


ax.plot(x[N2+more:], y[N2+more:],
'-',
color=lc[1],
linewidth=4,
label='Sliding')







plt.text(0.1, 0.96, r'$F_s$', fontsize=27)

plt.text(0.1, 0.67, r'$F_k$', fontsize=27)








ax.legend(loc=1, fontsize=25)
ax.grid()
ax.set_ylabel(r"Friction force", fontsize=30)
ax.set_xlabel(r"Time"  , fontsize=30)
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
ohg = 0.3 # arrow overhang

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



plt.savefig('steadySlide.pdf')
plt.show()
