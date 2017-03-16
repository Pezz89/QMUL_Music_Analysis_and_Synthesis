#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import scipy.io as io
import numpy as np
import sppysound.audiofile as af
import pdb

def main():
    if True:
        dpi = 100
        #fig = plt.figure(figsize=(1707/dpi, (846)/dpi), dpi=dpi)
        fig = plt.figure()
        s1 = io.loadmat("./analysis.mat")
        plt.legend(
            loc='upper center',
            bbox_to_anchor=(0.5, 1.05),
            ncol=3,
            fancybox=True,
            shadow=True
        )
        s = np.array(s1['f0'], dtype=float)
        s[s == 0] = np.nan
        #plt.plot(np.mod(12*np.log2(s/440.0)+69.0, 12))
        with af.AudioFile('./media/sguitar.aiff', 'r') as a:
            b = a.read_frames()
        x = np.linspace(0, len(s1['rms']), len(b))
        plt.plot(x, b)
        s = np.array(s1['rms'], dtype=float)
        plt.plot(s)
        s = np.array(s1['sf'], dtype=float)
        plt.plot(s)
        s = np.array(s1['f0'][0], dtype=float)
        plt.plot(s)
        pdb.set_trace()
        #plt.savefig('Segmentation.png')

        s1 = io.loadmat("./loopInds.mat")
        s = np.array(s1['loopInds'][0], dtype=int)
        plt.axvline(x[s[0]], color='r', linestyle='--')
        plt.axvline(x[s[1]], color='b', linestyle='--')

        plt.show()

    if False:
        dpi = 100
        fig = plt.figure(figsize=(1707/dpi, (846)/dpi), dpi=dpi)
        s1 = io.loadmat("./SpectralFlux.mat")
        a = s1['analysis'][0]
        plt.xlabel("Window Index")
        plt.xlim(0, len(a))
        plt.plot(a)

        plt.savefig('SpectralFlux.png')
    if False:
        dpi = 100
        fig = plt.figure(figsize=(1707/dpi, (846)/dpi), dpi=dpi)
        s1 = io.loadmat("./NormalisedAnalysis.mat")
        a = s1['analysis'][0]
        plt.xlabel("Window Index")
        plt.xlim(0, len(a))
        plt.plot(a)

        plt.savefig('NormalisedAnalysis.png')

if __name__ == "__main__":
    main()
