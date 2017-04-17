import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import rc
from sys import argv





def figsize(scale):
    fig_width_pt = 469.755                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size


pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": [],
    "font.monospace": [],
    "axes.labelsize": 12,               # LaTeX default is 10pt font.
    "text.fontsize": 12,
    "legend.fontsize": 10,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
        r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
        r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)








p = [2**i for i in range(7)]
#t = [61.55, 32.53, 18.16, 10.35, 6.09, 3.33, 2.02]
t = [3741, 1883, 1061, 634, 378.5, 208.25, 120.5]

speedup = np.zeros(7)

for i in xrange(len(t)):
    t[i] = int(t[i])*60 + t[i] - int(t[i])
    speedup[i] = float(t[0])/t[i]



fig, ax = plt.subplots(figsize=figsize(0.7))

plt.plot(p,p, '--', color="#DB7477", linewidth=3)
plt.hold('on')
plt.plot(p,speedup,
 '-h',
 color="#478684",
 linewidth=4,
 markersize=9,
 markerfacecolor='#3b617c',
 fillstyle='full')

ax.set_xscale('log', basex=2)
ax.set_yscale('log', basey=2)
ax.set_xticks([2**i for i in range(0,7)])
ax.set_yticks([2**i for i in range(0,7)])
ax.set_xticklabels([r'$%d$'%2**i for i in range(0,7)])
ax.set_yticklabels([r'$%d$'%2**i for i in range(0,7)])

plt.axis([min(p), max(p), min(p), max(p)])
plt.xlabel("Number of processors")
plt.ylabel("Speedup")
plt.grid('on')
ax.xaxis.set_label_coords(0.5, -0.1)

plt.tight_layout()


if len(argv)>1:
    filename = str(argv[1])
    plt.savefig(filename+'.pgf', format='pgf')
    plt.savefig(filename+'.pdf', format='pdf')

plt.show()
