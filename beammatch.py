from __future__ import division # needed for 1/2 = 0.5

from sympy import *
import numpy as np

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

## Putting in correct value
kval = 0.8
L_defocus = 0.4
L_drift = 0.4

m11_x = fodof_x[0][0]
m12_x = fodof_x[0][1]
m21_x = fodof_x[1][0]
m22_x = fodof_x[1][1]

m11_y = fodof_y[0][0]
m12_y = fodof_y[0][1]
m21_y = fodof_y[1][0]
m22_y = fodof_y[1][1]


x = -6.16013714551e-05 # this is just for one particle
xprime = -6.3139281937e-05
y = 1.62143430768e-05
yprime = 4.98176950809e-05
#x = -2.71243417537e-06 # mean
#xprime = -1.93168145531e-06
#y = -1.09113505722e-06
#yprime = 1.83854775315e-08
#x = 9.86402683437e-05 # std
#xprime = 9.8092029984e-05
#y = 0.00010098891309
#yprime = 0.000100293178212


print 'product',(2-m11_x**2-2*m12_x*m21_x-m22_x**2) # bigger than zero -> works

#print 's94 i wille'
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