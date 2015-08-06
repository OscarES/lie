from __future__ import division # needed for 1/2 = 0.5

from sympy.parsing.sympy_parser import parse_expr
from sympy import *
import numpy as np

L_0 = Symbol('L_0')
K = Symbol('K')
L_d = Symbol('L_d')

drift = np.array([
    [1,L_0,0,0],
    [0,1,0,0],
    [0,0,1,L_0],
    [0,0,0,1]
    ])
#print drift

focus = np.array([
    [cos(K*L_d/2),sin(K*L_d/2)/K,0,0],
    [-K*sin(K*L_d/2),cos(K*L_d/2),0,0],
    [0,0,cosh(K*L_d/2),sinh(K*L_d/2)/K],
    [0,0,K*sinh(K*L_d/2),cosh(K*L_d/2)]
    ])
#print focus

defocus = np.array([
    [cosh(K*L_d),sinh(K*L_d)/K,0,0],
    [K*sinh(K*L_d),cosh(K*L_d),0,0],
    [0,0,cos(K*L_d),sin(K*L_d)/K],
    [0,0,K*sin(K*L_d),cos(K*L_d)]
    ])
#print defocus

fodof = focus*drift*defocus*drift*focus
#print fodof

m11 = fodof[0][0]
m12 = fodof[0][1]
m21 = fodof[1][0]
m22 = fodof[1][1]

m33 = fodof[2][2]
m34 = fodof[2][3]
m43 = fodof[3][2]
m44 = fodof[3][3]

#mfodof = np.array([
#    [m11,m12,0,0],
#    [m21,m22,0,0],
#    [0,0,m33,m34],
#    [0,0,m43,m44]
#    ])
#print mfodof

T_x = np.array([
    [m11**2,2*m11*m12,m12**2],
    [m11*m21,m11*m22+m12*m21,m12*m22],
    [m21**2,2*m21*m22,m22**2]
    ])
T_y = np.array([
    [m33**2,2*m33*m34,m34**2],
    [m33*m43,m33*m44+m34*m43,m34*m44],
    [m43**2,2*m43*m44,m44**2]
    ])
#print T_x

V = np.array([
    [0,0,1/2],
    [0,1,0],
    [1/2,0,0]
    ])

lh = T_x.transpose()*V*T_x
print lh
print V