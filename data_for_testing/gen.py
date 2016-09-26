#!/usr/bin/python

import sys
import subprocess
import itertools
from numpy import *
import math
import matplotlib.pyplot as plt
import csv
from mpl_toolkits.mplot3d import Axes3D
   
meshAEraw = loadtxt('P_CFD_static')
nAe = meshAEraw.shape[1]
AmatBN = meshAEraw[:,0:3]
savetxt('undefAe.dat',AmatBN, fmt=' %-24.12E %-24.12E %-24.12E')
savetxt('pressure1D_static.dat',meshAEraw[:,3], fmt=' %-24.12E')


meshAEraw = loadtxt('P_CFD_mode4')
nAe = meshAEraw.shape[1]
AmatBNd = meshAEraw[:,0:3]
savetxt('defAe_mode4.dat',AmatBNd, fmt=' %-24.12E %-24.12E %-24.12E')
savetxt('pressure1D_mode4.dat',meshAEraw[:,3], fmt=' %-24.12E')


meshAEraw = loadtxt('P_CFD_mode1to4')
nAe = meshAEraw.shape[1]
AmatBNd = meshAEraw[:,0:3]
savetxt('defAe_mode1to4.dat',AmatBNd, fmt=' %-24.12E %-24.12E %-24.12E')
savetxt('pressure1D_mode1to4.dat',meshAEraw[:,3], fmt=' %-24.12E')



myfile = open('undefAe.csv', 'wb')
wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
wr.writerow(AmatBNd)

#xxx = loadtxt('../CSD_Module/cfd_data/defAe.dat')
#AmatBNd[:,0] = xxx[:,0]
#AmatBNd[:,1] = xxx[:,2]
#AmatBNd[:,2] = xxx[:,1]
#savetxt('defAe.dat',AmatBNd)

# just for looking at the transformation
#Hmat = loadtxt('../CSD_module/cfd_data/rbfHmatrix.dat')
#P_csd = dot(transpose(Hmat),meshAEraw[:,3])
#savetxt('pressure1D_CSD_2.dat',P_csd, fmt=' %-20.10E')

#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(AmatBN[:,0],AmatBN[:,1],AmatBN[:,2])
#ax.set_xlabel('X Label')
#ax.set_ylabel('Y Label')
#ax.set_zlabel('Z Label')
#AmatBN = loadtxt('../CSD_module/cfd_data/undefAe.dat')
#AmatBNd = loadtxt('../CSD_module/cfd_data/defAe.dat')

X = AmatBNd[:,0]
Y = AmatBNd[:,1]
Z = AmatBNd[:,2] 
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(X, Y, Z, c='r', marker='o')
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
#fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')



fig = plt.figure()
plt.plot(AmatBN[1,:], AmatBN[2,:])
plt.savefig('wing_box.pdf')
