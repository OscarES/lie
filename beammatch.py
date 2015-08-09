from __future__ import division # needed for 1/2 = 0.5

from sympy.parsing.sympy_parser import parse_expr
from sympy import *
import numpy as np

L_0 = Symbol('L_0')
K = Symbol('K')
L_d = Symbol('L_d')

drift_x = np.array([
    [1,L_0],
    [0,1]
    ])
drift_y = np.array([
    [1,L_0],
    [0,1],
    ])
#print drift

focus_x = np.array([
    [cos(K*L_d/2),sin(K*L_d/2)/K],
    [-K*sin(K*L_d/2),cos(K*L_d/2)]
    ])
focus_y = np.array([
    [cosh(K*L_d/2),sinh(K*L_d/2)/K],
    [K*sinh(K*L_d/2),cosh(K*L_d/2)]
    ])
#print focus

defocus_x = np.array([
    [cosh(K*L_d),sinh(K*L_d)/K],
    [K*sinh(K*L_d),cosh(K*L_d)]
    ])
defocus_y = np.array([
    [cos(K*L_d),sin(K*L_d)/K],
    [K*sin(K*L_d),cos(K*L_d)]
    ])
#print defocus

fodof_x = focus_x*drift_x*defocus_x*drift_x*focus_x
fodof_y = focus_y*drift_y*defocus_y*drift_y*focus_y
#print fodof

m11_x = fodof_x[0][0]
m12_x = fodof_x[0][1]
m21_x = fodof_x[1][0]
m22_x = fodof_x[1][1]

m11_y = fodof_y[0][0]
m12_y = fodof_y[0][1]
m21_y = fodof_y[1][0]
m22_y = fodof_y[1][1]

T_x = Matrix([
    [m11_x**2,2*m11_x*m12_x,m12_x**2],
    [m11_x*m21_x,m11_x*m22_x+m12_x*m21_x,m12_x*m22_x],
    [m21_x**2,2*m21_x*m22_x,m22_x**2]
    ])
T_y = Matrix([
    [m11_y**2,2*m11_y*m12_y,m12_y**2],
    [m11_y*m21_y,m11_y*m22_y+m12_y*m21_y,m12_y*m22_y],
    [m21_y**2,2*m21_y*m22_y,m22_y**2]
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

## Putting in correct value
kval = 1
L_defocus = 0.4
L_drift = 0.4
x = -6.16013714551e-05 # this is just for one particle, the best would be the mean value of all particles
xprime = -6.3139281937e-05
y = 1.62143430768e-05
yprime = 4.98176950809e-05
T_xmod = T_xmod.subs([(K,kval),(L_d,L_defocus),(L_0,L_drift)])
T_ymod = T_ymod.subs([(K,kval),(L_d,L_defocus),(L_0,L_drift)])

## Solve eqn system
evx = T_xmod.eigenvects()
#print evx
evy = T_ymod.eigenvects()
#print evy

## Extracting the twiss functions
beta_x = evx[0][2][0][0]
print 'beta_x',beta_x
alpha_x = evx[0][2][0][1]
print 'alpha_x',alpha_x
gamma_x = evx[0][2][0][2]
print 'gamma_x',gamma_x
epsilon_x = gamma_x*x**2+2*alpha_x*x*xprime+beta_x*xprime**2
print 'epsilon_x',epsilon_x

beta_y = evy[0][2][0][0]
print 'beta_y',beta_y
alpha_y = evy[0][2][0][1]
print 'alpha_y',alpha_y
gamma_y = evy[0][2][0][2]
print 'gamma_y',gamma_y
epsilon_y = gamma_y*y**2+2*alpha_y*y*yprime+beta_y*yprime**2
print 'epsilon_y',epsilon_y

## Save
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