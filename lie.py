import numpy as np
#from sympy import *
from sympy import Function, Symbol, symbols, summation
from sympy.mpmath import *

## Sum test
#x = np.array([1, 2, 4, 7, 11, 16], dtype=np.float)
#dx = np.gradient(x)
#
#print 'x=',x
#print 'dx=',dx
#
#nN = 3 # Number of dimensions
#result = np.zeros(nN)
#for i in range(1,nN+1):
#    result[i-1] = i
#resultsum = sum(result)
#print 'result=',result
#print 'resultsum=',resultsum


## Operator time!

# old dumb functions without sympy
#def lieop(nN,dfdq,dgdp,dfdp,dgdq):
#    result = np.zeros(nN)
#    for i in range(0,nN):
#        result[i] = dfdq[i]*dgdp[i]-dfdp[i]*dgdq[i]
#    resultsum = sum(result)
#    return resultsum
#
#def parder(fun,dir):
#    # something clever with sympy here
#    return fun/dir


## sympy test
#print('Sympy test...\n')
#x = symbols('x')
#print x
#a = Integral(cos(x)*exp(x), x)
#print a
#print Eq(a, a.doit())


print('- Lie calculations...\n')
q = Symbol('q')
p = Symbol('p')
f = Function('f')(q,p)
g = Function('g')(q,p)

dfdq = f.diff(q)
dfdp = f.diff(p)
dgdq = g.diff(q)
dgdp = g.diff(p)

def lieop(f,g):
    dfdq = f.diff(q)
    dfdp = f.diff(p)
    dgdq = g.diff(q)
    dgdp = g.diff(p)

    sumterm = dfdq*dgdp-dfdp*dgdq

    i, a, b = symbols('i a b', integer=True)
    colfcolg = summation(sumterm, (i,a,b))

    return colfcolg

colfcolg = lieop(f,g)

print colfcolg
#print colfcolg(1,2)

#nN2 = 6
#dfdq = np.array([1, 2, 4, 7, 11, 16], dtype=np.float)
#dgdp = np.array([1, 2, 4, 7, 11, 16], dtype=np.float)
#dfdp = np.array([1, 2, 4, 7, 11, 16], dtype=np.float)
#dgdq = np.array([1, 2, 4, 7, 11, 16], dtype=np.float)
#
#dfdp = np.multiply(2,dfdp)
#dgdq = np.multiply(2,dgdq)
#
#print lieop(nN2,dfdq,dgdp,dfdp,dgdq) 
