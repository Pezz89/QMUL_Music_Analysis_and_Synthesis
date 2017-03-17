#!/usr/bin/env python
from __future__ import division
import matplotlib.pyplot as plt
import scipy.io as io
import numpy as np
import sppysound.audiofile as af
import pdb

def main():
    # Plot analysis 1
    with af.AudioFile('./media/sguitar.aiff', 'r') as a:
        b = a.read_frames()[:, 0]

    if True:
        dpi = 100
        fig = plt.figure(figsize=(1707/dpi, (846)/dpi), dpi=dpi)
        s1 = io.loadmat("./analysis3.mat")
        x = s1['raw'][0]
        plt.plot(x, label='Original Samples')
        x = s1['clipped'][0]
        plt.plot(x, label='Segmented Samples')
        plt.xlabel('Time (samples)')
        plt.legend(
            loc='upper center',
            bbox_to_anchor=(0.5, 1.05),
            ncol=3,
            fancybox=True,
            shadow=True
        )
        plt.xlim([0, len(x)])
        plt.savefig('ZeroX.png')

    if False:
        dpi = 100
        fig = plt.figure(figsize=(1707/dpi, (846)/dpi), dpi=dpi)
        s1 = io.loadmat("./analysis1.mat")

        ax = plt.gca()
        ax2 = ax.twinx()
        #plt.plot(np.mod(12*np.log2(s/440.0)+69.0, 12))
        x = np.linspace(0, len(s1['rms']), len(b))
        lns1 = ax.plot(x, b, color='g')
        s = np.array(s1['rms'], dtype=float)
        lns2 = ax.plot(s, label='RMS', color='r')
        s = np.array(s1['sf'], dtype=float)
        lns3 = ax.plot(s, label='Spectral Flux', color='y')
        s = np.array(s1['f0'][0], dtype=float)
        s[s == 0] = np.nan
        lns4 = ax2.plot(s, label='f0', color='b')

        s1 = io.loadmat("./loopInds.mat")
        s = np.array(s1['loopInds'][0], dtype=int)
        lns5 = ax.axvline(x[s[0]], color='r', linestyle='--', label='Seg Start')
        lns6 = ax.axvline(x[s[1]], color='b', linestyle='--', label='Seg End')

        ax.set_xlim([0, 100])

        ax.set_xlabel('Time (samples)')
        ax.set_ylabel('Normalised values')
        ax2.set_ylabel('Frequency (Hz)')
        ax2.set_ylim([60, 300])

        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2,
            loc='upper center',
            bbox_to_anchor=(0.5, 1.10),
            ncol=3,
            fancybox=True,
            shadow=True
        )
        plt.savefig('Segmentation.png')

    # Plot analysis 2
    if False:
        dpi = 100
        fig = plt.figure(figsize=(1707/dpi, (846)/dpi), dpi=dpi)
        s1 = io.loadmat("./analysis2.mat")

        ax = plt.gca()
        #plt.plot(np.mod(12*np.log2(s/440.0)+69.0, 12))
        x = np.linspace(0, len(s1['rms']), len(b))
        ax.plot(x, b, color='g', label='Audio')
        s = np.array(s1['rms'], dtype=float)
        ax.plot(s, label='RMS Std-dev', color='r')
        s = np.array(s1['sf'], dtype=float)
        ax.plot(s, label='Spectral Flux Std-dev', color='y')
        s = np.array(s1['f0'][0], dtype=float)
        s[s == 0] = np.nan
        ax.plot(s, label='f0 Std-dev', color='b')

        ax.set_ylim([-0.001, 0.075])
        ax.set_xlim([0, 100])

        s1 = io.loadmat("./loopInds.mat")
        s = np.array(s1['loopInds'][0], dtype=int)

        ax.axvline(x[s[0]], color='r', linestyle='--', label='Seg Start')
        ax.axvline(x[s[1]], color='b', linestyle='--', label='Seg End')

        plt.xlabel('Time (samples)')
        plt.ylabel('Normalised values')
        plt.legend(
            loc='upper center',
            bbox_to_anchor=(0.5, 1.05),
            ncol=3,
            fancybox=True,
            shadow=True
        )
        plt.savefig('Segmentation2.png')

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
