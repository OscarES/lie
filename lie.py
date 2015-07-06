import numpy as np
from scipy.misc import *
#from scipy.linalg import *
from scipy import linalg
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
t = Symbol('t')

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
#quadham = -l/2*(k**2*(qx**2+qy**2+qz**2)+px**2+py**2+pz**2) # for Ems formalism, replace k with -k for defocus
quadham = -l/2*(k**2*(qx**2+qy**2)+px**2+py**2+pz**2) # for Ems formalism, replace k with -k for defocus. Without quad term in z dir

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


########### How many terms in the taylor series...
#I = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
#minusI = np.zeros([3,3])-I
##print minusI
##print np.dot(I,I)
#zero = np.zeros([3,3])
#Supper = np.concatenate((zero,I), axis=1)
#Sunder = np.concatenate((minusI,zero), axis=1)
#S = np.concatenate((Supper,Sunder), axis=0)
##print S
##print np.dot(S,S) # correct! (-I)
#
#sasexp = linalg.expm(S)
#left = np.dot(S,sasexp)
#right = np.dot(sasexp,S)
#
##print left
##print right
##print np.around(left-right) # proves (1.33)
#
## jacobian
#print ''
#print 'Jacobian...'
#dxdx = xofl.diff(x0)
#dxdp = xofl.diff(px)
#dpdx = xprimeofl.diff(x0)
#dpdp = xprimeofl.diff(px)
#
#print 'dxdx:',dxdx
#print 'dxdp:',dxdp
#print 'dpdx:',dpdx
#print 'dpdp:',dpdp
#print ''
#
## Either det(J) == 1
#dxdxtimdpdp = (dxdx*dpdp)
#dxdptimdpdx = (dxdp*dpdx)
##print 'dxdx*dpdp:',dxdxtimdpdp
##print 'dxdp*dpdx:',dxdptimdpdx
##print ''
#
#
## determinant by 'hand'
#detJ = dxdx*dpdp - dxdp*dpdx
#detJ = simplify(detJ)
#print 'detJ =', detJ # when this prints out 1 the simplecity is proved and no further terms are needed
##J = np.mat('[1 2;3 4]')
##J = np.mat('[dxdx dxdp; dpdx dpdp]')
##J = Matrix([[dxdx, dxdp], [dpdx, dpdp]])
#
##print J
#
## or np.dot(J.transpose(),np.dot(S,J)) == S. Jt*S*J
## -> dxdx*dpdp - dxdp*dpdx
#m12 = dxdx*dpdp - dxdp*dpdx
#m12 = simplify(m12)
#print 'm12:',m12
#
## As long as detJ and m12 gives: 1+ordo(order+1) it should be fine right?
## But 1+ordo(2) is the result for order 1, can that really work? Shouldn't the order be higher in order to describe the quadrupole?
#
### arbitrary functions and order into lietrans, one dim for simplicity
#print ''
#print 'arbitrary functions into lietransformation...'
#def simplelieop(f,g):
#    # The old way only one dim
#    dfdq = f.diff(q) # These four rows can be further developed to support diff for all three qs and three ps
#    dfdp = f.diff(p)
#    dgdq = g.diff(q)
#    dgdp = g.diff(p)
#    result = dfdq*dgdp-dfdp*dgdq
#    return result
#
### Lie transformation
## old version
#def simplelietransform(ham, vof0, t, order):
#    voft = vof0
#    for i in range(1,order+1):
#        lieterm = simplify(simplelieop(ham,vof0))
#        for j in range(0,i-1):
#            lieterm = simplify(simplelieop(ham,lieterm))
#
#        #print "lieterm:",lieterm
#        #voft = voft + t**i / factorial(i) * lieterm # for my formalism
#        voft = voft + lieterm / factorial(i) # for Ems formalism
#
#    return voft
#
#
#f = Function('f')(q,p)
#g = Function('g')(q,p)
##f = Function('f')(qx,qy,qz,px,py,pz)
##g = Function('g')(qx,qy,qz,px,py,pz)
##j = Symbol('j') # order
#j = 2
#result = simplelietransform(f, g, t, j) # v[0] is for q
#print result







##################### Classes and elements
print ''
print 'Classes and elements...'

def xfunFromHam(ham, order):
    transresultqx = lietransform(ham, v[0], t, order)
    xofl = substitution(transresultqx)
    return xofl

def xprimefunFromHam(ham, order):
    transresultpx = lietransform(ham, v[3], t, order)
    #xprimeofl = simplify(transresultpx/pz) # for my formalism
    xprimeofl = simplify(transresultpx) # for Ems formalism
    xprimeofl = substitution(xprimeofl)
    return xprimeofl

class Element:
    


    def __init__(self, name, ham, k, l, order):
        self.name = name
        self.ham = ham
        self.k = k
        self.l = l
        self.xfun = xfunFromHam(ham, order)
        self.xprimefun = xprimefunFromHam(ham, order)

    def printInfo(self):
        return self.name

order = 5
elem = Element('quad', quadham, k, l, order)
print elem.name
print elem.printInfo()
print 'xfun:',elem.xfun
print 'xprimefun:',elem.xprimefun

def straight():
    xp = np.arange(0,1)        
    x = np.arange(-0.05,0.06,0.01) #x and xp for when xp is 0
    return x,xp

def scanned():
    xp = np.arange(-0.0001,0.0001,0.00001)
    x = np.arange(-0.05,0.06,0.01)     # x and xp for when xp is scanned
    return x,xp
                            
def randomed():
    x = []
    for a in range (0, 10):
        x.append(random.uniform(-0.05, 0.05))
    
    xp = []
    for a in range (0, 10):
        xp.append(random.uniform(-0.00001, 0.00001)) #x and xp dfor random values
    return x,xp
                                                                                
def gaussian():
    x = []
    for a in range (0, 10):
        x.append(random.gauss(0, 0.01))
    xp = []
    for a in range (0, 10):
        xp.append(random.gauss(0, 0.00001))
    return x,xp
