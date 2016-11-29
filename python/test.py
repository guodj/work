#!/home/guod/anaconda3/bin/python3

import numpy as np
import matplotlib.pyplot as plt

def func1():
    x, y = np.arange(4), [0, 1, 1, 1]
    plt.plot(x, y)
    plt.xlim(0, 3)
    plt.ylim(0, 1.5)
    plt.xticks([])
    plt.yticks([])
    plt.xlabel('Time')
    plt.ylabel('Density difference')
    plt.show()


def func2():
    def f1():
        plt.plot(np.arange(10))
    def f2():
        plt.plot(np.arange(10)[::-1])
    f1()
    f2()
    plt.show()


def func3():
    import gitm
    import gitm_3D_global_plots as gpt

    a = gitm.GitmBin('/home/gdj/tmp/3DALL_t100516_113000.bin')
    gpt.gitm_single_nsglobal_3D_image('Temperature', a)
    plt.show()
    return a


def func4():
    """
    Test (lat, lt) to (mlat, mlt) conversions.

    height is needed for better conversion of mlat(qd)
    champ and grace files use qd coordinates.
    """
    import champ_grace as cg
    from apexpy import Apex as Apex
    import matplotlib.pyplot as plt
    a =  cg.ChampDensity('2005-1-1', '2005-1-2')
    mlat, mlt = Apex().convert(lat=a.lat, lon=a.long, source='geo', dest='mlt',
                               datetime=a.index, height=a.height)
    #mlat, mlon = Apex().convert(lat=a.lat, lon=a.long, source='geo', dest='qd',
    #                            datetime=a.index, height=a.height)
    mlat2 = np.array(a.Mlat)
    mlt2 = np.array(a.MLT)
    plt.plot(mlat-mlat2)
    plt.plot(abs(mlt-mlt2) % 24, '.')
    plt.show()
    return

if __name__ == '__main__':
    a = func3()
