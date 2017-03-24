class rgb2hex:
    def __init__(self, rgb):
        print '"#%02x%02x%02x"'.upper() %(rgb[0]*256, rgb[1]*256, rgb[2]*256)


if __name__ == '__main__':
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt

    f = sns.cubehelix_palette(8, start=0.5, rot=-.75)
    c=0
    for i in f:
        print c,
        rgb2hex(i)
        c+=1
    sns.palplot(f)
    plt.show()
