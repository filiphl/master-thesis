
class rgb2hex:
    def __init__(self, rgb):
        print '"#%02x%02x%02x"'.upper() %(rgb[0]*256, rgb[1]*256, rgb[2]*256)


if __name__ == '__main__':
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sys import argv


    if len(argv) > 3:
        nSteps = int(argv[1])
        s      = float(argv[2])
        r      = float(argv[3])
    else:
        nSteps = 8
        s = 0.5
        r = -0.75
    print "balle"
    f = sns.cubehelix_palette(sns.choose_cubehelix_palette())
    print "BALLE"
    c=0
    for i in f:
        print c,
        rgb2hex(i)
        c+=1
    sns.palplot(f)
    plt.show()
