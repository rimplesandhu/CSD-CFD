#!/usr/bin/python

import sys
import subprocess
import itertools
from numpy import *
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
# response plot
TR = loadtxt('../data/XT_S_TIP.dat')

Tlist = TR[:,0]
tip_le = TR[:,1]
tip_tr = TR[:,2]
tip_mid = TR[:,3]
heave = (tip_tr + tip_le)/2 - 0.0075
pitch = (tip_tr - tip_le)/0.3681984
fs = 5000.0
N = len(tip_le)

yf = fft.fft(tip_le)
psd_le = (2.0/N)*abs(yf[0:N/2])
yf = fft.fft(tip_tr)
psd_tr = (2.0/N)*abs(yf[0:N/2])
yf = fft.fft(tip_mid)
psd_mid = (2.0/N)*abs(yf[0:N/2])
f = linspace(0.0, fs/2, N/2)


fcc = 1
fig = plt.figure(1)
fig.set_size_inches(9,12)
plt.subplot(4,1,1)
plt.plot(Tlist,tip_tr*fcc,'-b',label="Trailing edge")
plt.plot(Tlist,tip_tr*fcc,'.b')
plt.plot(Tlist,tip_mid*fcc,'-k',label="Mid chord")
plt.plot(Tlist,tip_mid*fcc,'.k')
plt.plot(Tlist,tip_le*fcc,'-r',label="Leading edge")
plt.plot(Tlist,tip_le*fcc,'.r')
plt.ylabel('Vertical coordinates [m]')
plt.xlabel('Time [s]')
plt.title('AGARD Weakened model 3: Free vibration study (Restarting mode)')
plt.grid()
#plt.legend(loc='upper right')

plt.subplot(4,1,2)
plt.plot(f,psd_tr,'-b',label="Trailing edge")
plt.plot(f,psd_mid,'-k',label="Mid chord")
plt.plot(f,psd_le,'-r',label="Leading edge")
plt.ylabel('PSD')
plt.xlabel('Frequency (Hz)')
plt.title('')
plt.xlim([0, 200])
plt.grid()
plt.legend(loc='upper right')

plt.subplot(4,1,3)
plt.plot(Tlist,heave,'-b',label="Heave")
plt.plot(Tlist,heave,'.b')
plt.grid()
plt.xlabel('Time [s]')
plt.legend(loc='upper right')

plt.subplot(4,1,4)
plt.plot(Tlist,pitch,'-b',label="Pitch")
plt.plot(Tlist,pitch,'.b')
plt.xlabel('Time [s]')
plt.title('')
plt.grid()
plt.legend(loc='upper right')


plt.savefig('summary_response.pdf')

