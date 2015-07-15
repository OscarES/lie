import numpy as np
from scipy.misc import *
#from scipy.linalg import *
from scipy import linalg
from sympy.parsing.sympy_parser import parse_expr
from sympy import *
import math
import random
#printing.init_printing(use_latex='mathjax') # latex output in ipython notebook
#from sympy import Function, Symbol, symbols, summation
#from sympy.mpmath import *

###################### Symbol setup, define all symbols here
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
p = Symbol('p')
g = Symbol('g')
l = Symbol('l') # arbitrary length
k = Symbol('k')
t = Symbol('t')

## Hamiltonians
#ems way
driftham = -l/2*(px**2 + py**2 + pz**2)
quadham = -l/2*(k**2*(qx**2+qy**2)+px**2+py**2+pz**2) # (Ems formalism), replace k with -k for defocus. Without quad term in z dir
quadhamdefocus = -l/2*(-k**2*(qx**2+qy**2)+px**2+py**2+pz**2) # (Ems formalism), replace k with -k for defocus. Without quad term in z dir

#old way
#driftham = px**2 / (2*m) + py**2 / (2*m) + pz**2 / (2*m) # for my formalism
#quadham = px**2 / (2*m) + py**2 / (2*m) + pz**2 / (2*m) + q*g*pz*qx**2 / (2*m) - q*g*pz*qy**2 / (2*m) # for some reason the denominator always has to be last in the term. for my formalism

## One can transform arbitrary functions if one would like so
# arbitrary funs here
#f = Function('f')(q,p)
#g = Function('g')(q,p)

###################################### Operators and transforms here
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
def lietransform(ham, vof0, t, order):
    voft = vof0
    for i in range(1,order+1):
        lieterm = simplify(lieop(ham,vof0))
        for j in range(0,i-1):
            lieterm = simplify(lieop(ham,lieterm))

        #voft = voft + t**i / factorial(i) * lieterm # for my formalism
        voft = voft + lieterm / factorial(i) # for Ems formalism

    return voft


# substitution
def substitution(expr):
    expr = expr.subs(t, l*m/pz)
    expr = expr.subs(qx, x0)
    expr = expr.subs(px/pz, xprime)
    expr = expr.subs(q*g/pz,k**2)
    return expr




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

## Matrix construction, only works with my formalism
#xofl = substitution(transresultqx)
#print "xofl:", xofl
#xprimeofl = substitution(transresultpx/pz) # for my formalism
#xprimeofl = substitution(transresultpx) # for Ems formalism
#print "xprimeofl:", xprimeofl
#print "Matrix construction..."
#m11 = removeprefixones(selectterm(xofl,x0))
#m12 = removeprefixones(selectterm(xofl,xprime))
#m21 = removeprefixones(selectterm(xprimeofl,x0))
#print xprimeofl
#print selectterm(xprimeofl,xprime)
#m22 = removeprefixones(selectterm(xprimeofl,xprime))
#
#driftmatrix = [[m11, m12], [m21, m22]]
#print driftmatrix




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

    def evaluate(self, (mulxin, mulxpin)): # sending in xin and xpin as a vector (same for return) allows "recursive" calls and a lattice can be constructed
        xoutmul = list()
        xpoutmul = list()
        for xin, xpin in zip(mulxin, mulxpin):
            expr = self.xfun.subs([(x0,xin), (px,xpin), (k,self.k), (l,self.l)])
            xout = expr.evalf()
            xoutmul.append(xout)
            expr = self.xprimefun.subs([(x0,xin), (px,xpin), (k,self.k), (l,self.l)])
            xpout = expr.evalf()
            xpoutmul.append(xpout)
        return (xoutmul,xpoutmul)


#print cos(3.14/2).evalf() # python uses radians

############# Randoms
## Old dust
#print ''
#print 'Randoms...'

# with arange
#def straight():
#    xp = np.arange(0,1)        
#    x = np.arange(-0.05,0.06,0.01) #x and xp for when xp is 0
#    return x,xp
#
#def scanned():
#    xp = np.arange(-0.0001,0.0001,0.00001)
#    x = np.arange(-0.05,0.06,0.01)     # x and xp for when xp is scanned
#    return x,xp
#                            
#def randomed():
#    x = []
#    for a in range (0, 10):
#        x.append(random.uniform(-0.05, 0.05))
#    
#    xp = []
#    for a in range (0, 10):
#        xp.append(random.uniform(-0.00001, 0.00001)) #x and xp dfor random values
#    return x,xp
#                                                                                
#def gaussian():
#    x = []
#    for a in range (0, 10):
#        x.append(random.gauss(0, 0.01))
#    xp = []
#    for a in range (0, 10):
#        xp.append(random.gauss(0, 0.00001))
#    return x,xp

## New stuff
def straight(particles):
    xp = np.linspace(0,0,particles)        
    x = np.linspace(-1,1,particles) #x and xp for when xp is 0
    return x,xp

def scanned(particles):
    xp = np.linspace(-0.0001,0.0001,particles)
    x = np.linspace(-0.05,0.05,particles)     # x and xp for when xp is scanned
    return x,xp
 
 # old should not be used                           
def randomed(particles):
    # new
    x = [random.uniform(-0.05, 0.05) for _ in xrange(particles)]
    xp = [random.uniform(-0.00001, 0.00001) for _ in xrange(particles)]
    # old
    #x = []
    #for a in range (0, nbrofparticles):
    #    x.append(random.uniform(-0.05, 0.05))
    #
    #xp = []
    #for a in range (0, nbrofparticles):
    #    xp.append(random.uniform(-0.00001, 0.00001)) #x and xp dfor random values
    return x,xp
 
# slow                                                                               
#def gaussian():
    # new
    #x = [random.gauss(0, 0.01) for _ in xrange(nbrofparticles)]
    #xp = [random.gauss(0, 0.00001) for _ in xrange(nbrofparticles)]
    # old
    #x = []
    #for a in range (0, nbrofparticles):
    #    x.append(random.gauss(0, 0.01))
    #xp = []
    #for a in range (0, nbrofparticles):
    #    xp.append(random.gauss(0, 0.00001))
    #return x,xp

def gaussian(particles):
    x = np.random.normal(0,0.001,particles)
    xp = np.random.normal(0,0.000001,particles)
    y = np.random.normal(0,0.001,particles)
    yp = np.random.normal(0,0.000001,particles)
    return x,xp,y,yp

##print 'straight:', straight(nbrofparticles)
##print 'scanned:', scanned(nbrofparticles)
##print 'randomed:', randomed(nbrofparticles)
#print 'gaussian:', gaussian(nbrofparticles)

## Multiple particles and lattice construction
print 'Multiple particles and lattice construction...'
nbrofparticles = raw_input('Enter number of particles:')
nbrofparticles = int(nbrofparticles)

order = 5
myK = 2
myQuadL = 0.5
myDriftL = 1.83

myQuad = Element('quad', quadham, myK, myQuadL, order)
myDrift = Element('drift', driftham, 0, myDriftL, order)

lattice = list()
lattice.append(myQuad)
lattice.append(myDrift)
#print lattice[0].name
#print lattice[1].name

def evalLattice(lattice,(xin,xpin)):
    xout, xpout = xin,xpin
    for elem in lattice:
        xout, xpout = elem.evaluate((xout,xpout))
    return xout, xpout

#print 'evalLattice:',evalLattice(lattice,([1],[0]))
#print 'evalLattice with mul particles',evalLattice(lattice,(mulx,mulxp))

##### Saved and loaded data
datamode = raw_input('Save or load data (S/l):')
if datamode != 's' and datamode != 'l':
    datamode = 's'
datafile = raw_input('Enter file name:')
if len(datafile) < 1 : datafile = "test.txt"
if datamode == 's':
    x, xp, y, yp = gaussian(nbrofparticles)
    dt = np.dtype([('x', 'd'), ('xp', 'd'), ('y', 'd'), ('yp', 'd'), ('alpha', 'd'), ('beta', 'd'), ('epsilon', 'd')])
    a = np.zeros(nbrofparticles, dt)
    a['x'] = x
    a['xp'] = xp
    a['y'] = y
    a['yp'] = yp
    a['alpha'] = 0
    a['beta'] = 0
    a['epsilon'] = 0

    np.savetxt(datafile, a, '%10s')

elif datamode == 'l':
    try:
        x, xp, y, yp, alpha, beta, epsilon = np.loadtxt(datafile,unpack = True)
    except:
        print 'Bad datafile!'
        quit()
    #print stx, stxp, scx, scxp, rax, raxp, gax, gaxp

###### eval data, lattice defined above
#print 'Straight...'
#print 'Input x and xp:', stx, stxp
#stxo, stxpo = evalLattice(lattice,(stx,stxp))
#print 'Output x and xp:', stxo, stxpo
#
#print 'Scanned...'
#print 'Input x and xp:', scx, scxp
#scxo, scxpo = evalLattice(lattice,(scx,scxp))
#print 'Output x and xp:', scxo, scxpo
#
#print 'Random...'
#print 'Input x and xp:', rax, raxp
#raxo, raxpo = evalLattice(lattice,(rax,raxp))
#print 'Output x and xp:', raxo, raxpo

#print 'Gaussian...'
#print 'Input x and xp:', x, xp
#gaxo, gaxpo = evalLattice(lattice,(x,xp))
#print 'Output x and xp:', gaxo, gaxpo


##### FODO
print 'FODO...'
fF = Element('quad', quadham, myK, myQuadL, order)
oO1 = Element('drift', driftham, 0, myDriftL, order)
dD = Element('quaddefocus', quadhamdefocus, myK, myQuadL, order)
oO2 = Element('drift', driftham, 0, myDriftL, order)

fodoLattice = list()
fodoLattice.append(fF)
fodoLattice.append(oO1)
fodoLattice.append(dD)
fodoLattice.append(oO2)

# input from randoms and loads above

#stxoFODO, stxpoFODO = evalLattice(fodoLattice,(stx,stxp))
stxoFODO, stxpoFODO = evalLattice(fodoLattice,(x,xp))

#print 'Output x and xp:', stxoFODO, stxpoFODO

outputfile = raw_input('Enter file for output data:')
if len(outputfile) < 1 : outputfile = "out.txt"
if len(outputfile) > 0:
    dt = np.dtype([('x', 'd'), ('xp', 'd')]) #('x', 'd'), ('xp', 'd'), ('alpha', 'd'), ('beta', 'd'), ('epsilon', 'd')])
    a = np.zeros(nbrofparticles, dt)
    a['x'] = stxoFODO
    a['xp'] = stxpoFODO
#    a['y'] = styoFODO
#    a['yp'] = stypoFODO
#    a['alpha'] = 0
#    a['beta'] = 0
#    a['epsilon'] = 0

    np.savetxt(outputfile, a, '%10s')

##### Focal length
#print ''
#print 'Focal length...'
#singleQuad = list()
#singleQuad.append(fF)
#stxoQ, stxpoQ = evalLattice(singleQuad,(stx,stxp))
#print 'beforeQ:',stx, stxp 
#print 'afterQ:',stxoQ, stxpoQ 
#scxoQ, scxpoQ = evalLattice(singleQuad,(scx,scxp))
#raxoQ, raxpoQ = evalLattice(singleQuad,(rax,raxp))
#gaxoQ, gaxpoQ = evalLattice(singleQuad,(gax,gaxp))
#
#stFlength = -stxoQ[0]/stxpoQ[0]
#print 'Straight focal length:', stFlength
#scFlength = -scxoQ[0]/scxpoQ[0]
#print 'Scanned focal length:', scFlength
#raFlength = -raxoQ[0]/raxpoQ[0]
#print 'Random focal length:', raFlength
#gaFlength = -gaxoQ[0]/gaxpoQ[0]
#print 'Gaussian focal length:', gaFlength

##### Bugs to fix and TODOs
# Fix k substitution (remove the substitution and just calculate for k's value) -> seems like the expansion doesn't hold for k > 1...

# Finish writing the output "spool" so that all required data is there

# Start simulating 2D particles