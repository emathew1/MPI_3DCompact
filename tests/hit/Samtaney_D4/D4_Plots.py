# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#%%

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

#%%

sdir = '/Users/emathew1/Google Drive/MPI_3DCompact/tests/hit/Samtaney_D4'

A = np.loadtxt(sdir+'/D4_Rhoprime.dat')

timeRhoPrime = A[:,0]
RhoPrime = A[:,1]

#%%

sdir = '/Users/emathew1/Google Drive/MPI_3DCompact/tests/hit/Samtaney_D4'

A = np.loadtxt(sdir+'/D4_Mt.dat')

timeMt = A[:,0]
Mt = A[:,1]

#%%

sdir = '/Users/emathew1/Google Drive/MPI_3DCompact/tests/hit/Samtaney_D4'

A = np.loadtxt(sdir+'/D4_K.dat')

timeK = A[:,0]
K = A[:,1]

#%%

ddir = '/Users/emathew1/Google Drive/MPI_3DCompact/tests/hit/Samtaney_D4'

A = np.loadtxt(ddir+'/turbdata.out')


#%%
plt.figure()
plt.plot(timeRhoPrime, RhoPrime,'o-', A[:,0]/A[0,1], A[:,7]/A[0,3]/A[0,3],'-')
plt.axis([0, 7.5, 0, 0.5])
plt.xlabel(r't/\tau')
plt.ylabel('rho\'/Mt(0)')

#%%
plt.figure()
plt.plot(timeMt, Mt, 'o-', A[:,0]/A[0,1], A[:,3]/A[0,3],'-')
plt.axis([0, 7.5, 0, 1.0])
plt.xlabel(r't/\tau')
plt.ylabel('Mt/Mt(0)')

#%%
plt.figure()
plt.plot(timeK, K, 'o-', A[:,0]/A[0,1], A[:,10]/A[0,10])
plt.axis([0, 7.5, 0, 1.0])
plt.xlabel(r't/\tau')
plt.ylabel('K(t)/K_0')

#%%

#plt.plot(A[:,0]/A[0,1], A[:,8]/A[0,9])
