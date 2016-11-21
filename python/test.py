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

if __name__ == '__main__':
    func2()
