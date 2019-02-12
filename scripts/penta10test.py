#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 21:24:14 2019

@author: edwin
"""

# %%

import numpy as np
import matplotlib.pyplot as plt


#%%

N = 25

A = np.zeros((N,N))
B = np.zeros((N,N))

#Interior points

#Penta 10
beta  =   1.0/20.0;
alpha =   1.0/2.0;
a     =  17.0/12.0;
b     = 101.0/150.0;
c     =   1.0/100.0;

#Tri6
#beta  =   0;
#alpha =   1.0/3.0;
#a     =  14.0/9.0;
#b     =   1.0/9.0;
#c     =   0.0/100.0;


#Point 0, N-1

#Penta 10
#beta1  = 28.0;
#alpha1 = 16.0;
#a1     = -1181.0/280.0;
#b1     =  -892.0/35.0;
#c1     =    77.0/5.0;
#d1     =    56.0/3.0;
#e1     =   -35.0/6.0;
#f1     =    28.0/15.0;
#g1     =    -7.0/15.0;
#h1     =     8.0/105.0;
#i1     =    -1.0/168.0;
#j1     = 0.0

#Tri 10
beta1  = 0;
alpha1 = 9.0;
a1     = -2485.0/649.0;
b1     = -1809.0/280.0;
c1     =    18.0/1.0;
d1     =   -14.0/1.0;
e1     =    21.0/2.0;
f1     =   -63.0/10.0;
g1     =    14.0/5.0;
h1     =    -6.0/7.0;
i1     =     9.0/56.0;
j1     =    -1.0/72.0;

#Penta 8
#beta1  = 15.0;
#alpha1 = 8.0;
#a1     = -79.0/20.0;
#b1     = -77.0/5.0;
#c1     =    55.0/4.0;
#d1     =    20.0/3.0;
#e1     =    -5.0/4.0;
#f1     =     1.0/5.0;
#g1     =    -1.0/60.0;
#h1     =    -0.0/7.0;
#i1     =     0.0/56.0;
#j1     =    -0.0/72.0;

#Tri 6
#beta1  = 0;
#alpha1 = 5.0;
#a1     = -197.0/60.0;
#b1     = -5.0/12.0;
#c1     =    5.0/1.0;
#d1     =   -5.0/3.0;
#e1     =    5.0/12.0;
#f1     =   -1.0/20.0;
#g1     =    0.0#14.0/5.0;
#h1     =    0.0#-6.0/7.0;
#i1     =    0.0# 9.0/56.0;
#j1     =    0.0#-1.0/72.0;

#Point 1, N-2

#Penta 10
beta2    = 5.0/3.0;
alpha2_1 = 1.0/21.0;
alpha2_2 = 3.0;
a2       = -544.0/2581.0;
b2       =  -39.0/20.0;
c2       =  -17.0/20.0;
d2       =   95.0/36.0;
e2       =    5.0/12.0;
f2       =   -1.0/20.0;
g2       =    1.0/180.0;
h2       =   -1.0/2940.0;

#Tri6
#beta2    = 0.0;
#alpha2_1 = 1.0/8.0;
#alpha2_2 = 3.0/4.0;
#a2       =    -43/96.0;
#b2       =   -5.0/6.0;
#c2       =    9.0/8.0;
#d2       =    1.0/6.0;
#e2       =   -1.0/96.0;
#f2       =   0.0#-1.0/20.0;
#g2       =   0.0# 1.0/180.0;
#h2       =   0.0#-1.0/2940.0;

#Penta6
#beta2    = 0.0;
#alpha2_1 = 1.0/8.0;
#alpha2_2 = 3.0/4.0;
#a2       =    -43/96.0;
#b2       =   -5.0/6.0;
#c2       =    9.0/8.0;
#d2       =    1.0/6.0;
#e2       =   -1.0/96.0;
#f2       =   0.0#-1.0/20.0;
#g2       =   0.0# 1.0/180.0;
#h2       =   0.0#-1.0/2940.0;

#Point 2, N-3

#Penta 10
#beta3_1  = 1.0/90.0;
#beta3_2  = 1.0;
#alpha3_1 = 4.0/15.0;
#alpha3_2 = 8.0/9.0;
#a3       =  -34.0/675.0;
#b3       = -127.0/225.0;
#c3       =   -7.0/12.0;
#d3       =   20.0/27.0;
#e3       =    4.0/9.0;
#f3       =    1.0/75.0;
#g3       =   -1.0/2700.0;
#h3 = 0
#i3 = 0

#Tri 10
#beta3_1  = 0/90.0;
#beta3_2  = 0.0;
#alpha3_1 = 1.0/7.0;
#alpha3_2 = 1.0;
#a3       =  -1.0/168.0;
#b3       = -433.0/980.0;
#c3       =   -19.0/20.0;
#d3       =   21.0/20.0;
#e3       =    5.0/12.0;
#f3       =    -1.0/20.0;
#g3       =   1.0/60.0;
#h3       =  -1.0/420.0
#i3       =  1.0/5880.0


#Tri 6 interior
#beta3_1  = 0.0;
#beta3_2  = 0.0;
#alpha3_1 = 1.0/3.0;
#alpha3_2 = 1.0/3.0;
#a3       = -1.0/9.0/4.0;
#b3       = -14/9.0/2;
#c3       =   0.0;
#d3       =  14/9/2.0;
#e3       =    1/9/4.0;
#f3       =    0;
#g3       =   0;
#h3 = 0
#i3 = 0

#Penta 8 interior
beta3_1  = 1.0/36.0;
beta3_2  = 1.0/36.0;
alpha3_1 = 4.0/9.0;
alpha3_2 = 4.0/9.0;
a3       = -25.0/54.0/4.0;
b3       = -40.0/27.0/2.0;
c3       =   0.0;
d3       =  40.0/27.0/2.0;
e3       =  25.0/54.0/4.0;
f3       =    0;
g3       =   0;
h3 = 0
i3 = 0

for i in range(0,N):
    A[i,i] = 1

for i in range(0,N-1):
    A[i,i+1] = alpha
    A[i+1,i] = alpha
    B[i,i+1] = a/2
    B[i+1,i] = -a/2
    
for i in range(0,N-2):
    A[i,i+2] = beta
    A[i+2,i] = beta
    B[i,i+2] = b/4
    B[i+2,i] = -b/4
    
for i in range(0,N-3):
    B[i,i+3] =  c/6
    B[i+3,i] = -c/6    
    
Bperiodic = np.copy(B)
Aperiodic = np.copy(A)
Bperiodic[0,-1] = -a/2
Bperiodic[0,-2] = -b/4
Bperiodic[0,-3] = -c/6    
Bperiodic[1,-1] = -b/4
Bperiodic[1,-2] = -c/6
Bperiodic[2,-1] = -c/6

Bperiodic[-1,0] = a/2
Bperiodic[-1,1] = b/4
Bperiodic[-1,2] = c/6
Bperiodic[-2,0] = b/4
Bperiodic[-2,1] = c/6
Bperiodic[-3,0] = c/6   

Aperiodic[0,-1] = alpha
Aperiodic[0,-2] = beta
Aperiodic[1,-1] = beta
Aperiodic[-1,0] = alpha
Aperiodic[-1,1] = beta
Aperiodic[-2,0] = beta 

B[0,0] = a1
B[0,1] = b1
B[0,2] = c1
B[0,3] = d1
B[0,4] = e1
B[0,5] = f1
B[0,6] = g1
B[0,7] = h1
B[0,8] = i1
B[0,9] = j1

B[-1,-1] = -a1
B[-1,-2] = -b1
B[-1,-3] = -c1
B[-1,-4] = -d1
B[-1,-5] = -e1
B[-1,-6] = -f1
B[-1,-7] = -g1
B[-1,-8] = -h1
B[-1,-9] = -i1
B[-1,-10]= -j1    

B[1,0] = a2
B[1,1] = b2
B[1,2] = c2
B[1,3] = d2
B[1,4] = e2
B[1,5] = f2
B[1,6] = g2
B[1,7] = h2

B[-2,-1] = -a2
B[-2,-2] = -b2
B[-2,-3] = -c2
B[-2,-4] = -d2
B[-2,-5] = -e2
B[-2,-6] = -f2
B[-2,-7] = -g2
B[-2,-8] = -h2

B[2,0] = a3
B[2,1] = b3
B[2,2] = c3
B[2,3] = d3
B[2,4] = e3
B[2,5] = f3
B[2,6] = g3
B[2,7] = h3
B[2,8] = i3


B[-3,-1] = -a3
B[-3,-2] = -b3
B[-3,-3] = -c3
B[-3,-4] = -d3
B[-3,-5] = -e3
B[-3,-6] = -f3
B[-3,-7] = -g3
B[-3,-8] = -h3
B[-3,-9] = -i3


A[0,1] = alpha1
A[0,2] = beta1
A[-1,-2] = alpha1
A[-1,-3] = beta1

A[1,0] = alpha2_1
A[1,2] = alpha2_2
A[1,3] = beta2

A[-2, -1] = alpha2_1
A[-2, -3] = alpha2_2
A[-2, -4] = beta2

A[2,0] = beta3_1
A[2,1] = alpha3_1
A[2,3] = alpha3_2
A[2,4] = beta3_2

A[-3, -1] = beta3_1
A[-3, -2] = alpha3_1
A[-3, -4] = alpha3_2
A[-3, -5] = beta3_2

x = np.zeros((N,1))
xperiodic = np.zeros((N,1))
for i in range(0,N):
    x[i] = np.sin((float(i)/float(N-1))*2.0*np.pi)
    xperiodic[i] = np.sin((float(i)/float(N))*2.0*np.pi)

# %%

dx = 2.0*np.pi/(N-1)
dxp = 2.0*np.pi/N
RHS = (1/dx)*np.matmul(B,x)
RHSperiodic = (1/dxp)*np.matmul(Bperiodic,xperiodic)

#%%
y = np.linalg.solve(A,RHS)
yperiodic = np.linalg.solve(Aperiodic, RHSperiodic)
