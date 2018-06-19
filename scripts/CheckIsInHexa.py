# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

A = np.transpose(np.asarray([0, 0, 0]))
B = np.transpose(np.asarray([0, 1, 0]))
C = np.transpose(np.asarray([1, 0, 0]))
D = np.transpose(np.asarray([1, 1, 0]))
E = np.transpose(np.asarray([0, 0, 1]))
F = np.transpose(np.asarray([0, 1, 1]))
G = np.transpose(np.asarray([1, 0, 1]))
H = np.transpose(np.asarray([1, 1, 1]))

P = (A + B + C + D + E + F + G + H)/8
VolP = 0.0

#Face 1 - ABCD
#Tet 1 - ABCP
Mat1 = np.matrix([[A[0], A[1], A[2], 1.0],
                  [B[0], B[1], B[2], 1.0],
                  [C[0], C[1], C[2], 1.0],
                  [P[0], P[1], P[2], 1.0]])

VolP += abs(np.linalg.det(Mat1))/6.0

#Tet 2 - BDCP
Mat2 = np.matrix([[B[0], B[1], B[2], 1.0],
                  [D[0], D[1], D[2], 1.0],
                  [C[0], C[1], C[2], 1.0],
                  [P[0], P[1], P[2], 1.0]])

VolP += abs(np.linalg.det(Mat2))/6.0


#Face 2 - EFGH
#Tet 3 - EFGP
Mat3 = np.matrix([[G[0],G[1],G[2], 1.0],
                  [F[0],F[1],F[2], 1.0],
                  [E[0],E[1],E[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP += abs(np.linalg.det(Mat3))/6.0


#Tet 4 - FGHP
Mat4 = np.matrix([[F[0],F[1],F[2], 1.0],
                  [G[0],G[1],G[2], 1.0],
                  [H[0],H[1],H[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP += abs(np.linalg.det(Mat4))/6.0


#Face 3 - ABEF
#Tet 5 - ABEP
Mat5 = np.matrix([[E[0],E[1],E[2], 1.0],
                  [B[0],B[1],B[2], 1.0],
                  [A[0],A[1],A[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP += abs(np.linalg.det(Mat5))/6.0


#Tet 6 - BEFP
Mat6 = np.matrix([[B[0],B[1],B[2], 1.0],
                  [E[0],E[1],E[2], 1.0],
                  [F[0],F[1],F[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP += abs(np.linalg.det(Mat6))/6.0


#Face 4 - CDGH
#Tet 7 - CDGP
Mat7 = np.matrix([[C[0],C[1],C[2], 1.0],
                  [D[0],D[1],D[2], 1.0],
                  [G[0],G[1],G[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP += abs(np.linalg.det(Mat7))/6.0


#Tet 8 - DGHP
Mat8 = np.matrix([[H[0],H[1],H[2], 1.0],
                  [G[0],G[1],G[2], 1.0],
                  [D[0],D[1],D[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP += abs(np.linalg.det(Mat8))/6.0


#Face 5 - ACGE
#Tet 9 - ACGP
Mat9 = np.matrix([[A[0],A[1],A[2], 1.0],
                  [C[0],C[1],C[2], 1.0],
                  [G[0],G[1],G[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP += abs(np.linalg.det(Mat9))/6.0


#Tet 10 - CGEP
Mat10 = np.matrix([[C[0],C[1],C[2], 1.0],
                  [G[0],G[1],G[2], 1.0],
                  [E[0],E[1],E[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP += abs(np.linalg.det(Mat10))/6.0


#Face 6 - BDFH
#Tet 11 - BDFP
Mat11 = np.matrix([[F[0],F[1],F[2], 1.0],
                   [D[0],D[1],D[2], 1.0],
                   [B[0],B[1],B[2], 1.0],
                   [P[0],P[1],P[2], 1.0]])

VolP += abs(np.linalg.det(Mat11))/6.0


#Tet 12 - DFHP
Mat12 = np.matrix([[D[0],D[1],D[2], 1.0],
                   [F[0],F[1],F[2], 1.0],
                   [H[0],H[1],H[2], 1.0],
                   [P[0],P[1],P[2], 1.0]])

VolP += abs(np.linalg.det(Mat12))/6.0



#Test point

P = [0.0, 0.0, 0.1]

VolP2 = 0.0

#Face 1 - ABCD
#Tet 1 - ABCP
Mat1 = np.matrix([[A[0], A[1], A[2], 1.0],
                  [B[0], B[1], B[2], 1.0],
                  [C[0], C[1], C[2], 1.0],
                  [P[0], P[1], P[2], 1.0]])

VolP2 += abs(np.linalg.det(Mat1))/6.0

#Tet 2 - BDCP
Mat2 = np.matrix([[B[0], B[1], B[2], 1.0],
                  [D[0], D[1], D[2], 1.0],
                  [C[0], C[1], C[2], 1.0],
                  [P[0], P[1], P[2], 1.0]])

VolP2 += abs(np.linalg.det(Mat2))/6.0


#Face 2 - EFGH
#Tet 3 - EFGP
Mat3 = np.matrix([[G[0],G[1],G[2], 1.0],
                  [F[0],F[1],F[2], 1.0],
                  [E[0],E[1],E[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP2 += abs(np.linalg.det(Mat3))/6.0


#Tet 4 - FGHP
Mat4 = np.matrix([[F[0],F[1],F[2], 1.0],
                  [G[0],G[1],G[2], 1.0],
                  [H[0],H[1],H[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP2 += abs(np.linalg.det(Mat4))/6.0


#Face 3 - ABEF
#Tet 5 - ABEP
Mat5 = np.matrix([[E[0],E[1],E[2], 1.0],
                  [B[0],B[1],B[2], 1.0],
                  [A[0],A[1],A[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP2 += abs(np.linalg.det(Mat5))/6.0


#Tet 6 - BEFP
Mat6 = np.matrix([[B[0],B[1],B[2], 1.0],
                  [E[0],E[1],E[2], 1.0],
                  [F[0],F[1],F[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP2 += abs(np.linalg.det(Mat6))/6.0


#Face 4 - CDGH
#Tet 7 - CDGP
Mat7 = np.matrix([[C[0],C[1],C[2], 1.0],
                  [D[0],D[1],D[2], 1.0],
                  [G[0],G[1],G[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP2 += abs(np.linalg.det(Mat7))/6.0


#Tet 8 - DGHP
Mat8 = np.matrix([[H[0],H[1],H[2], 1.0],
                  [G[0],G[1],G[2], 1.0],
                  [D[0],D[1],D[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP2 += abs(np.linalg.det(Mat8))/6.0


#Face 5 - ACGE
#Tet 9 - ACGP
Mat9 = np.matrix([[A[0],A[1],A[2], 1.0],
                  [C[0],C[1],C[2], 1.0],
                  [G[0],G[1],G[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP2 += abs(np.linalg.det(Mat9))/6.0


#Tet 10 - CGEP
Mat10 = np.matrix([[C[0],C[1],C[2], 1.0],
                  [G[0],G[1],G[2], 1.0],
                  [E[0],E[1],E[2], 1.0],
                  [P[0],P[1],P[2], 1.0]])

VolP2 += abs(np.linalg.det(Mat10))/6.0


#Face 6 - BDFH
#Tet 11 - BDFP
Mat11 = np.matrix([[F[0],F[1],F[2], 1.0],
                   [D[0],D[1],D[2], 1.0],
                   [B[0],B[1],B[2], 1.0],
                   [P[0],P[1],P[2], 1.0]])

VolP2 += abs(np.linalg.det(Mat11))/6.0


#Tet 12 - DFHP
Mat12 = np.matrix([[D[0],D[1],D[2], 1.0],
                   [F[0],F[1],F[2], 1.0],
                   [H[0],H[1],H[2], 1.0],
                   [P[0],P[1],P[2], 1.0]])

VolP2 += abs(np.linalg.det(Mat12))/6.0
