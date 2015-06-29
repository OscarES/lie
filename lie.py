import numpy as np
from scipy.misc import *
from sympy.parsing.sympy_parser import parse_expr
from sympy import *
#printing.init_printing(use_latex='mathjax') # latex output in ipython notebook
#from sympy import Function, Symbol, symbols, summation
#from sympy.mpmath import *

qx = Symbol('qx')
qy = Symbol('qy')
qz = Symbol('qz')
px = Symbol('px')
py = Symbol('py')
pz = Symbol('pz')

qi = [qx, qy, qz] 
pi = [px, py, pz] 

v = qi + pi

x = Symbol('x')
xprime = Symbol('xprime')
y = Symbol('y')
yprime = Symbol('yprime')
z = Symbol('z')
zprime = Symbol('zprime')

x0 = Symbol('x0')
xprime0 = Symbol('xprime0')

#m = 1.67262178*10**-27 # mass of proton
m = Symbol('m') # arbitrary mass

q = Symbol('q')
g = Symbol('g')
l = Symbol('l') # arbitrary length
k = Symbol('k')

print('- Lie calculations...\n')
# variable definitions
q = Symbol('q')
p = Symbol('p')

# arbitrary funs here
#f = Function('f')(q,p)
#g = Function('g')(q,p)


def lieop(f,g):

    # The old way only one dim
#    dfdq = f.diff(q) # These four rows can be further developed to support diff for all three qs and three ps
#    dfdp = f.diff(p)
#    dgdq = g.diff(q)
#    dgdp = g.diff(p)
#    sumterm = dfdq*dgdp-dfdp*dgdq
#
#    i, a, b = symbols('i a b', integer=True)
    #colfcolg = summation(sumterm, (i,a,b))

    # New way all three dims
#    gradf = gradient(f,qi) # gradient isn't available
    dfdqx = f.diff(qx) 
    dfdpx = f.diff(px)
    dgdqx = g.diff(qx)
    dgdpx = g.diff(px)
    sumtermx = dfdqx*dgdpx-dfdpx*dgdqx

    dfdqy = f.diff(qy) 
    dfdpy = f.diff(py)
    dgdqy = g.diff(qy)
    dgdpy = g.diff(py)
    sumtermy = dfdqy*dgdpy-dfdpy*dgdqy

    dfdqz = f.diff(qz) 
    dfdpz = f.diff(pz)
    dgdqz = g.diff(qz)
    dgdpz = g.diff(pz)
    sumtermz = dfdqz*dgdpz-dfdpz*dgdqz

    colfcolg = sumtermx + sumtermy + sumtermz

    return colfcolg


## Lie transformation
# old version
def lietransform(ham, vof0, t, order):
    voft = vof0
    for i in range(1,order+1):
        lieterm = simplify(lieop(ham,vof0))
        for j in range(0,i-1):
            lieterm = simplify(lieop(ham,lieterm))

        #print "lieterm:",lieterm
        #voft = voft + t**i / factorial(i) * lieterm # for my formalism
        voft = voft + lieterm / factorial(i) # for Ems formalism

    return voft

# new version
#def lietransform(ham, vof0, t, order):
#    return

# x = mm, v = mum
#ham = p**2 / (2*m)
#vof0 = sin(q) + cos(p)
#t = 10 # arbitrary time
#order = 2
#
#transresult = lietransform(ham, vof0, t, order)
#print transresult

## elementary elements
# drift
print "Drift..."
#driftham = px**2 / (2*m) + py**2 / (2*m) + pz**2 / (2*m) # for my formalism
driftham = -l/2*(px**2 + py**2 + pz**2)

#t = 10
t = Symbol('t')
order = 5

transresultqx = lietransform(driftham, v[0], t, order)
#print "transresultqx:", transresultqx
transresultpx = lietransform(driftham, v[3], t, order)
#print "transresultpx:", transresultpx

# substitution
def substitution(expr):
    expr = expr.subs(t, l*m/pz)
    expr = expr.subs(qx, x0)
    expr = expr.subs(px/pz, xprime)
    expr = expr.subs(q*g/pz,k**2)
    return expr

xofl = substitution(transresultqx)
print "xofl:", xofl

#xprimeofl = substitution(transresultpx/pz) # for my formalism
xprimeofl = substitution(transresultpx) # for Ems formalism
print "xprimeofl:", xprimeofl


# select only terms with a certain var in expr, assumes naively that var is the last symbol in each term
def selectterm(expr, var):
    strexpr = str(expr)
    strvar = str(var)
    idx = strexpr.find(strvar)
    if idx < 0:
        return 0
    if strexpr[idx-1] == " ": # takes care when there is no coefficients
        return 1
    if len(strexpr) == len(strvar): # If the only thing in the expression was the var
        return 1
    strexpr = strexpr[:idx] # gets rid of the variable and terms to the right (risky way to select the desired term)
    idx2 = strexpr.rfind(" ")
    strexpr = strexpr[idx2+1:] 
    strexpr = strexpr.strip()
    if len(strexpr) == 0:
        strexpr = "0"
    if strexpr.endswith("+") or strexpr.endswith("-") or strexpr.endswith("*") or strexpr.endswith("/"):
        strexpr = strexpr[:-1] 
    newexpr = parse_expr(strexpr)
    return newexpr

def removeprefixones(expr): # still in an early form, BUG: it removes 1.0 if there is a space behind it. BUG: Does not work with several 1.0, especially if there is an early 1.0* then the other 1.0 won't be removed and we're stuck in a loop
    strexpr = str(expr)
    idx = 0
    while True:
        idx = strexpr.find("1.0")
        if idx < 0:
            break
        if not strexpr[idx+3].isdigit(): # makes sure that for example 1.00002 doesn't get removed
            strexpr = strexpr[0:idx] + strexpr[idx+4:]
    newexpr = parse_expr(strexpr)
    return newexpr

# Matrix construction, only works with my formalism
print "Matrix construction..."
m11 = removeprefixones(selectterm(xofl,x0))
m12 = removeprefixones(selectterm(xofl,xprime))
m21 = removeprefixones(selectterm(xprimeofl,x0))
print xprimeofl
print selectterm(xprimeofl,xprime)
m22 = removeprefixones(selectterm(xprimeofl,xprime))

driftmatrix = [[m11, m12], [m21, m22]]
print driftmatrix


# quad
print "Quad..."
#quadham = px**2 / (2*m) + py**2 / (2*m) + pz**2 / (2*m) + q*g*pz*qx**2 / (2*m) - q*g*pz*qy**2 / (2*m) # for some reason the denominator always has to be last in the term. for my formalism
quadham = -l/2*(k**2*(qx**2+qy**2+qz**2)+px**2+py**2+pz**2) # for Ems formalism, replace k with -k for defocus

#t = 10
t = Symbol('t')
order = 5

transresultqx = lietransform(quadham, v[0], t, order)
#print "transresultqx:", transresultqx
transresultpx = lietransform(quadham, v[3], t, order)
#print "transresultpx:", transresultpx

# substitution
xofl = substitution(transresultqx)
#fourier() here to go from taylor to sin & cos ?
print "xofl:", xofl

#xprimeofl = simplify(transresultpx/pz) # for my formalism
xprimeofl = simplify(transresultpx) # for Ems formalism
xprimeofl = substitution(xprimeofl)
print "xprimeofl:", xprimeofl

#for the other dimension replace x with y and z, will make functions for this later
