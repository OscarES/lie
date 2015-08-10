from __future__ import division # needed for 1/2 = 0.5

from sympy.parsing.sympy_parser import parse_expr
from sympy import *
import numpy as np
from numpy import linalg as la

#L_0 = Symbol('L_0')
#K = Symbol('K')
#L_d = Symbol('L_d')

L_0 = 0.4
K = 0.8
L_d = 0.4

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
    [-K*sin(K*L_d),cos(K*L_d)]
    ])
#print defocus

fodof_x = np.dot(focus_x,np.dot(drift_x,np.dot(defocus_x,np.dot(drift_x,focus_x))))
fodof_y = np.dot(focus_y,np.dot(drift_y,np.dot(defocus_y,np.dot(drift_y,focus_y))))
#fodof_y = focus_y*drift_y*defocus_y*drift_y*focus_y
#print fodof

## Putting in correct value
kval = 0.8
L_defocus = 0.4
L_drift = 0.4

m11_x = fodof_x[0][0]#.subs([(K,kval),(L_d,L_defocus),(L_0,L_drift)])
m12_x = fodof_x[0][1]#.subs([(K,kval),(L_d,L_defocus),(L_0,L_drift)])
m21_x = fodof_x[1][0]#.subs([(K,kval),(L_d,L_defocus),(L_0,L_drift)])
m22_x = fodof_x[1][1]#.subs([(K,kval),(L_d,L_defocus),(L_0,L_drift)])

#print 'm11_x',m11_x
#print 'm12_x',m12_x
#print 'm21_x',m21_x
#print 'm22_x',m22_x
# WORKS SO FAR

m11_y = fodof_y[0][0]#.subs([(K,kval),(L_d,L_defocus),(L_0,L_drift)])
m12_y = fodof_y[0][1]#.subs([(K,kval),(L_d,L_defocus),(L_0,L_drift)])
m21_y = fodof_y[1][0]#.subs([(K,kval),(L_d,L_defocus),(L_0,L_drift)])
m22_y = fodof_y[1][1]#.subs([(K,kval),(L_d,L_defocus),(L_0,L_drift)])

#print 'm11_x',m11_x
#print 'm12_x',m12_x
#print 'm21_x',m21_x
#print 'm22_x',m22_x


x = -6.16013714551e-05 # this is just for one particle, the best would be the mean value of all particles
xprime = -6.3139281937e-05
y = 1.62143430768e-05
yprime = 4.98176950809e-05


print 'product',(2-m11_x**2-2*m12_x*m21_x-m22_x**2) # bigger than zero -> works

print 's94 i wille'
beta_x = 2*m12_x/sqrt(2 - m11_x**2 - 2*m12_x*m21_x - m22_x**2)
print 'beta_x',beta_x
alpha_x = (m11_x - m22_x)*beta_x/(2*m12_x)
print 'alpha_x',alpha_x
gamma_x = (1+alpha_x**2)/beta_x
print 'gamma_x',gamma_x
epsilon_x = gamma_x*x**2+2*alpha_x*x*xprime+beta_x*xprime**2
print 'epsilon_x',epsilon_x

beta_y = 2*m12_y/sqrt(2 - m11_y**2 - 2*m12_y*m21_y - m22_y**2)
print 'beta_y',beta_y
alpha_y = (m11_y - m22_y)*beta_y/(2*m12_y)
print 'alpha_y',alpha_y
gamma_y = (1+alpha_y**2)/beta_y
print 'gamma_y',gamma_y
epsilon_y = gamma_y*y**2+2*alpha_y*y*yprime+beta_y*yprime**2
print 'epsilon_y',epsilon_y

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
print 'end wille'
quit()

T_x = Matrix([
    [m11_x**2,-2*m11_x*m12_x,m12_x**2],
    [-m11_x*m21_x,m11_x*m22_x+m12_x*m21_x,-m12_x*m22_x],
    [m21_x**2,-2*m21_x*m22_x,m22_x**2]
    ])
T_y = Matrix([
    [m11_y**2,-2*m11_y*m12_y,m12_y**2],
    [-m11_y*m21_y,m11_y*m22_y+m12_y*m21_y,-m12_y*m22_y],
    [m21_y**2,-2*m21_y*m22_y,m22_y**2]
    ])
#print T_x
#T_x = Matrix([
#    [m11_x**2,2*m11_x*m12_x,m12_x**2],
#    [m11_x*m21_x,m11_x*m22_x+m12_x*m21_x,m12_x*m22_x],
#    [m21_x**2,2*m21_x*m22_x,m22_x**2]
#    ])
#T_y = Matrix([
#    [m11_y**2,2*m11_y*m12_y,m12_y**2],
#    [m11_y*m21_y,m11_y*m22_y+m12_y*m21_y,m12_y*m22_y],
#    [m21_y**2,2*m21_y*m22_y,m22_y**2]
#    ])
#print T_x

## WRONG
#print 'T_x\n',T_x
T_xmod = T_x-Matrix([
    [1,0,0],
    [0,1,0],
    [0,0,1]
    ])
#print 'T_xmod\n',T_xmod
#print T_xmod

T_ymod = T_y-Matrix([
    [1,0,0],
    [0,1,0],
    [0,0,1]
    ])
#print T_ymod



###### Particles, but this is a chicken and the egg problem! x,xp,y,yp -> twiss -> epsilon -> x,xp,y,yp
#x = -6.16013714551e-05 # this is just for one particle, the best would be the mean value of all particles
#xprime = -6.3139281937e-05
#y = 1.62143430768e-05
#yprime = 4.98176950809e-05
#x = -2.71243417537e-06 # mean
#xprime = -1.93168145531e-06
#y = -1.09113505722e-06
#yprime = 1.83854775315e-08
#x = 9.86402683437e-05 # std
#xprime = 9.8092029984e-05
#y = 0.00010098891309
#yprime = 0.000100293178212
#T_x = T_x.subs([(K,kval),(L_d,L_defocus),(L_0,L_drift)])
#T_y = T_y.subs([(K,kval),(L_d,L_defocus),(L_0,L_drift)])

## Solve eqn system WRONG
evx = T_xmod.eigenvects()
#evax,evex = np.linalg.eig(T_xmod)
print 'evx',evx
evy = T_ymod.eigenvects()
#evay,evey = la.eig(T_ymod)
print 'evy',evy
quit()
## Extracting the twiss functions WRONG
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

## Twiss
#print 'm11_x',m11_x
#print 'm12_x',m12_x
#print 'm21_x',m21_x
#print 'm22_x',m22_x
#print 'product',(2-m11_x**2-2*m12_x*m21_x-m22_x**2)
#beta_x = (2*m12_x)/sqrt(2-m11_x**2-2*m12_x*m21_x-m22_x**2)
#print 'beta_x',beta_x
#alpha_x = (m11_x-m22_x)/sqrt(2-m11_x**2-2*m12_x*m21_x-m22_x**2)
#print 'alpha_x',alpha_x
#beta_y = (2*m12_y)/sqrt(2-m11_y**2-2*m12_y*m21_y-m22_y**2)
#print 'beta_y',beta_y
#alpha_y = (m11_y-m22_y)/sqrt(2-m11_y**2-2*m12_y*m21_y-m22_y**2)
#print 'alpha_y',alpha_y

## Save




V = np.array([
    [0,0,1/2],
    [0,1,0],
    [1/2,0,0]
    ])

#lh = T_x.transpose()*V*T_x
#print lh
#print V
