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
        fig = plt.figure(figsize=(1707/dpi, (846)/dpi), dpi=dpi)
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
        with af.AudioFile('./media/vibrato.wav', 'r') as a:
            b = a.read_frames()
        #plt.plot(np.linspace(0, len(s1['f0']), len(b)), b)
        pdb.set_trace()
        plt.plot(s)
        #plt.plot(s1['rms'])
        #plt.plot(s1['sf'])
        plt.show()
        #plt.savefig('Segmentation.png')


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
