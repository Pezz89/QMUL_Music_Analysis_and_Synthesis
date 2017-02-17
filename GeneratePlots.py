#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import scipy.io as io
import numpy as np
import pdb

def main():
    s1 = io.loadmat("./vars.mat")
    a = s1['a'][0]
    analysis = s1['analysis'][0]
    transience_s = s1['transience_s'][0]
    transience_e = s1['transience_e'][0]
    WLen = s1['WLen'][0]
    win_count = s1['win_count'][0]
    n1 = s1['n1'][0]


    plt.plot((np.arange(win_count)*n1)+WLen/2,analysis)
    for i in transience_s:
        plt.axvline(i, color='r', linestyle='--')
    for i in transience_e:
        plt.axvline(i, color='b', linestyle='--')
    plt.show()

if __name__ == "__main__":
    main()
