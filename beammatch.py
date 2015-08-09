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

T_x = Matrix([
    [m11**2,2*m11*m12,m12**2],
    [m11*m21,m11*m22+m12*m21,m12*m22],
    [m21**2,2*m21*m22,m22**2]
    ])
T_y = Matrix([
    [m33**2,2*m33*m34,m34**2],
    [m33*m43,m33*m44+m34*m43,m34*m44],
    [m43**2,2*m43*m44,m44**2]
    ])
#print T_x

T_xmod = T_x-Matrix([
    [1,0,0],
    [0,1,0],
    [0,0,1]
    ])
#print T_xmod

T_ymod = T_y-Matrix([
    [1,0,0],
    [0,1,0],
    [0,0,1]
    ])
#print T_ymod

# Putting in correct value
kval = 1
L_defocus = 0.4
L_drift = 0.4
x = -6.16013714551e-05 # this is just for one particle, the best would be the mean value of all particles
xprime = -6.3139281937e-05
y = 1.62143430768e-05
yprime = 4.98176950809e-05
T_xmod = T_xmod.subs([(K,kval),(L_d,L_defocus),(L_0,L_drift)])
T_ymod = T_ymod.subs([(K,kval),(L_d,L_defocus),(L_0,L_drift)])

# Solve eqn system
evx = T_xmod.eigenvects()
#print evx

beta_x = evx[0][2][0][0]
print 'beta_x',beta_x
alpha_x = evx[0][2][0][1]
print 'alpha_x',alpha_x
gamma_x = evx[0][2][0][2]
print 'gamma_x',gamma_x
epsilon_x = gamma_x*x**2+2*alpha_x*x*xprime+beta_x*xprime**2
print 'epsilon_x',epsilon_x

evy = T_ymod.eigenvects()
#print evy

beta_y = evy[0][2][0][0]
print 'beta_y',beta_y
alpha_y = evy[0][2][0][1]
print 'alpha_y',alpha_y
gamma_y = evy[0][2][0][2]
print 'gamma_y',gamma_y
epsilon_y = gamma_y*y**2+2*alpha_y*y*yprime+beta_y*yprime**2
print 'epsilon_y',epsilon_y

# Save
datafiletwiss = raw_input('Enter twiss file name:')
dtb = np.dtype([('alpha_x', 'd'), ('beta_x', 'd'), ('epsilon_x', 'd'), ('alpha_y', 'd'), ('beta_y', 'd'), ('epsilon_y', 'd')])
b = np.zeros(1,dtb)
b['alpha_x'] = alpha_x
b['beta_x'] = beta_x
b['epsilon_x'] = epsilon_x
b['alpha_y'] = alpha_y
b['beta_y'] = beta_y
b['epsilon_y'] = epsilon_y
np.savetxt(datafiletwiss, b, '%10s')



V = np.array([
    [0,0,1/2],
    [0,1,0],
    [1/2,0,0]
    ])

#lh = T_x.transpose()*V*T_x
#print lh
#print V