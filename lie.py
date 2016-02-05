from __future__ import division # needed for 1/2 = 0.5
import numpy as np
from scipy.misc import *
#from scipy.linalg import *
from scipy import linalg
from sympy.parsing.sympy_parser import parse_expr
from sympy import *
import math
import random
import matplotlib.pyplot as plt

#printing.init_printing(use_latex='mathjax') # latex output in ipython notebook
#from sympy import Function, Symbol, symbols, summation
#from sympy.mpmath import *

###################### Symbol setup, define all symbols here
plt.rcParams.update({'font.size': 22})
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


l = Symbol('l') # arbitrary length
k = Symbol('k')


## Hamiltonians
driftham = -l/2*(px**2 + py**2 + pz**2)
quadham = -l/2*(k**2*(qx**2-qy**2)+px**2+py**2+pz**2) # replace k with -k for defocus. Without quad term in z dir
quadhamdefocus = -l/2*(-k**2*(qx**2-qy**2)+px**2+py**2+pz**2) # replace k with -k for defocus. Without quad term in z dir
sextupoleham = -l/2*(2/3*k**2*(qx**3-3*qx*qy**2)+(px**2+py**2)) # should the ps' perhaps be divided by 2 as in nonlinear2013_3.pdf? That division is assumed to be the l/2 in the beginning, . k is actually k**2
octupoleham = -l/2*(2/4*k*(qx**4-6*qx**2*qy**2+qy**4)+(px**2+py**2)) # same decision as above


###################################### Operators and transforms here
def lieop(f,g):
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
def lietransform(ham, vof0, order):#,t):, #zzz
    voft = vof0
    for i in range(1,order+1):
        lieterm = simplify(lieop(ham,vof0))
        for j in range(0,i-1):
            lieterm = simplify(lieop(ham,lieterm))

        #voft = voft + t**i / factorial(i) * lieterm # for my formalism, #zzz
        voft = voft + lieterm / factorial(i) # for Ems formalism, shouldn't each term also 

    return voft


########### How many terms in the taylor series... Remake this code so that it finds the optimal order for the run (determinant of the Jacobian differ by less than 10^-5 or 10^-6 from 1) Use and print this order. The suggested changes should be inserted into the element class so that each element has a unique order. Or maybe The optimal order for each pole should be calculated and then just inserted, i.e. 1 for drift and ...
def findOrder(ham,K,L,acceptableError,approxOffsetQ,approxOffsetP): # Since this also find the functions maybe it should also return them. Since why do the same thing twice?
    if K*L > 1:
        print 'K*L larger than 1!'
        quit()
    order = 1
    while True:

        xfun = funFromHam(ham, order, qx)
        xprimefun = funFromHam(ham, order, px)
        yfun = funFromHam(ham, order, qy)
        yprimefun = funFromHam(ham, order, py)

        xf = xfun.subs([(k,K),(l,L)])
        xpf = xprimefun.subs([(k,K), (l,L)])
        yf = yfun.subs([(k,K),(l,L)])
        ypf = yprimefun.subs([(k,K), (l,L)])

        J = np.array([[xf.diff(qx),xf.diff(qy),xf.diff(px),xf.diff(py)],[yf.diff(qx),yf.diff(qy),yf.diff(px),yf.diff(py)],[xpf.diff(qx),xpf.diff(qy),xpf.diff(px),xpf.diff(py)],[ypf.diff(qx),ypf.diff(qy),ypf.diff(px),ypf.diff(py)]]) # Note that here I deviate from the normal order x,xp,y,yp

        ## Approximative qx and qy set to 0.1 (px and py always disapper in diffs but if needed the same could be done to them)
        ## Maybe the approxOffset should be set to the mean of the Gaussian
        #approxOffsetQ = 0.01
        #approxOffsetP = 0.00001
        for i in range(0,4):
            for j in range(0,4):
                J[i,j] = J[i,j].subs([(qx,approxOffsetQ),(qy,approxOffsetQ),(px,approxOffsetP),(py,approxOffsetP)])

        detJ = np.linalg.det(J)

        error = detJ - 1
        if abs(error) < acceptableError:
            break
        order += 1
        #print 'order',order
        #print 'error',error

    return order

###### Class for elements
def funFromHam(ham, order, vof0):
    transresult = lietransform(ham, vof0, order)
    fun = simplify(transresult)
    return fun

class Element:
    def __init__(self, name, ham, kval, lval, order):
        self.name = name
        self.ham = ham
        self.L = lval

        # Same this done 4 times, time for a new function?
        self.xfun = funFromHam(ham, order, qx)
        self.xprimefun = funFromHam(ham, order, px)
        self.yfun = funFromHam(ham, order, qy)
        self.yprimefun = funFromHam(ham, order, py)

        self.xf = self.xfun.subs([(k,kval),(l,lval)])
        self.xpf = self.xprimefun.subs([(k,kval), (l,lval)])
        self.yf = self.yfun.subs([(k,kval),(l,lval)])
        self.ypf = self.yprimefun.subs([(k,kval), (l,lval)])

        self.xf = lambdify((qx,px,qy,py),self.xf, "numpy")
        self.xpf = lambdify((qx,px,qy,py),self.xpf, "numpy")
        self.yf = lambdify((qx,px,qy,py),self.yf, "numpy")
        self.ypf = lambdify((qx,px,qy,py),self.ypf, "numpy")

    def printInfo(self):
        return self.name

    def evaluate(self, (mulxin, mulxpin, mulyin, mulypin)): # sending in xin and xpin as a vector (same for return) allows "recursive" calls and a lattice can be constructed
        xout = self.xf(mulxin, mulxpin, mulyin, mulypin)
        xpout = self.xpf(mulxin, mulxpin, mulyin, mulypin)
        yout = self.yf(mulxin, mulxpin, mulyin, mulypin)
        ypout = self.ypf(mulxin, mulxpin, mulyin, mulypin)

        return (xout,xpout,yout,ypout)


#print cos(3.14/2).evalf() # python uses radians

## Plot def
def plotPhaseSpace(x,xp,y,yp):
    plt.subplot(121)
    plt.plot(x,xp,'ro')
    plt.xlabel('x')
    plt.ylabel('xp')
    
    plt.subplot(122)
    plt.plot(y,yp,'ro')
    plt.xlabel('y')
    plt.ylabel('yp')

    plt.show()

############# Randoms
def straight(particles):
    x = np.linspace(-1,1,particles) #x and xp for when xp is 0
    xp = np.linspace(0,0,particles)
    y = np.linspace(-1,1,particles) #y and yp for when yp is 0
    yp = np.linspace(0,0,particles)
    return x,xp,y,yp

def scanned(particles):
    x = np.linspace(-0.05,0.05,particles)     # x and xp for when xp is scanned
    xp = np.linspace(-0.0001,0.0001,particles)
    y = np.linspace(-0.05,0.05,particles)     # y and yp for when yp is scanned
    yp = np.linspace(-0.0001,0.0001,particles)
    return x,xp,y,yp
                           
def randomed(particles):
    x = [random.uniform(-0.05, 0.05) for _ in xrange(particles)]
    xp = [random.uniform(-0.00001, 0.00001) for _ in xrange(particles)]
    y = [random.uniform(-0.05, 0.05) for _ in xrange(particles)]
    yp = [random.uniform(-0.00001, 0.00001) for _ in xrange(particles)]
    return x,xp,y,yp

def gaussian(particles):
    x = np.random.normal(0,0.001,particles)
    xp = np.random.normal(0,0.000001,particles)
    y = np.random.normal(0,0.001,particles)
    yp = np.random.normal(0,0.000001,particles)
    return x,xp,y,yp

def gaussiantwiss(particles, alpha, beta, epsilon):
    xi = np.random.normal(0,1,particles)
    xip = np.random.normal(0,1,particles)

    M = np.array([
        [1/sqrt(beta*epsilon), 0],
        [alpha/sqrt(beta*epsilon), sqrt(beta/epsilon)]
        ])

    Minv = np.linalg.inv(M)

    x = np.zeros(particles)
    xp = np.zeros(particles)

    for i in range(particles):
        x[i] = Minv[0,0]*xi[i] + Minv[0,1]*xip[i]
        xp[i] = Minv[1,0]*xi[i] + Minv[1,1]*xip[i]

    return x,xp

################## Multiple particles and lattice construction
print 'Multiple particles and lattice construction...'

#def rms(a):
#    sq = np.power(a,2)
#    me = np.mean(sq)
#    root = np.sqrt(me)
#    return root

def evalLattice(lattice,(xin,xpin,yin,ypin)):
    xout, xpout, yout, ypout = xin,xpin,yin,ypin
    envx = np.zeros((len(lattice)+1,2))
    envy = np.zeros((len(lattice)+1,2))

    envx[0][1] = np.std(xin)
    envy[0][1] = np.std(yin)

    elementnbr = 1

    for elem in lattice:
        xout, xpout, yout, ypout = elem.evaluate((xout,xpout,yout,ypout))
        envx[elementnbr][0] = envx[elementnbr-1][0] + elem.L
        envx[elementnbr][1] = np.std(xout)

        envy[elementnbr][0] = envy[elementnbr-1][0] + elem.L
        envy[elementnbr][1] = np.std(yout)

        elementnbr += 1

    return xout, xpout, yout, ypout, envx, envy

##### Saved and loaded data
datamodepart = raw_input('Save or load particles (S/l):')
if datamodepart != 's' and datamodepart != 'l':
    datamodepart = 's'
datamodetwiss = raw_input('Save or load twiss (S/l):')
if datamodepart != 's' and datamodepart != 'l':
    datamodepart = 's'
datafilepart = raw_input('Enter particle file name:')
if len(datafilepart) < 1 : datafilepart = "inpart.txt"

datafiletwiss = raw_input('Enter twiss file name:')
if len(datafiletwiss) < 1 : datafiletwiss = "intwiss.txt"

dta = np.dtype([('x', 'd'), ('xp', 'd'), ('y', 'd'), ('yp', 'd')])
dtb = np.dtype([('alpha_x', 'd'), ('beta_x', 'd'), ('epsilon_x', 'd'), ('alpha_y', 'd'), ('beta_y', 'd'), ('epsilon_y', 'd')])

if datamodepart == 's':
    nbrofparticles = raw_input('Enter number of particles:')
    if len(nbrofparticles) < 1 : nbrofparticles = "1000"
    nbrofparticles = int(nbrofparticles)
    
    # can only save twiss when saving particles
    if datamodetwiss == 's':
        alpha_x = raw_input('Enter alpha_x:') # set this to 0
        if len(alpha_x) < 1 : alpha_x = "0"
        alpha_x = float(alpha_x)
        beta_x = raw_input('Enter beta_x:') # set this to 0.001
        if len(beta_x) < 1 : beta_x = "5"
        beta_x = float(beta_x)
        epsilon_x = raw_input('Enter epsilon_x:') # set this to 0.000001
        if len(epsilon_x) < 1 : epsilon_x = "0.00001"
        epsilon_x = float(epsilon_x)
        alpha_y = raw_input('Enter alpha_y:') # set this to 0
        if len(alpha_y) < 1 : alpha_y = "0"
        alpha_y = float(alpha_y)
        beta_y = raw_input('Enter beta_y:') # set this to 0.001
        if len(beta_y) < 1 : beta_y = "5"
        beta_y = float(beta_y)
        epsilon_y = raw_input('Enter epsilon_y:') # set this to 0.000001
        if len(epsilon_y) < 1 : epsilon_y = "0.00001"
        epsilon_y = float(epsilon_y)
    elif datamodetwiss == 'l':
        try:
            alpha_x, beta_x, epsilon_x, alpha_y, beta_y, epsilon_y = np.loadtxt(datafiletwiss,unpack = True)
        except:
            print 'Bad datafile!'
            quit()
    
    #x, xp, y, yp = gaussian(nbrofparticles) # uncoupled x and xp, y and yp
    x, xp = gaussiantwiss(nbrofparticles,alpha_x,beta_x,epsilon_x) # coupled x and xp
    y, yp = gaussiantwiss(nbrofparticles,alpha_y,beta_y,epsilon_y) # coupled y and yp
    
    #plotPhaseSpace(x,xp,y,yp)

    a = np.zeros(nbrofparticles, dta)
    b = np.zeros(1,dtb)

    a['x'] = x
    a['xp'] = xp
    a['y'] = y
    a['yp'] = yp
    b['alpha_x'] = alpha_x
    b['beta_x'] = beta_x
    b['epsilon_x'] = epsilon_x
    b['alpha_y'] = alpha_y
    b['beta_y'] = beta_y
    b['epsilon_y'] = epsilon_y


    np.savetxt(datafilepart, a, '%10s')
    np.savetxt(datafiletwiss, b, '%10s')

elif datamodepart == 'l':
    try:
        x, xp, y, yp= np.loadtxt(datafilepart,unpack = True)
        alpha_x, beta_x, epsilon_x, alpha_y, beta_y, epsilon_y = np.loadtxt(datafiletwiss,unpack = True)

        nbrofparticles = len(x)
    except:
        print 'Bad datafile!'
        quit()



##### FODO
print 'FODO...'
#myKsquared = 0.64 # Wille and lie formalism, with Ems this is myK, zzzzz
myK = 0.8
#myK = sqrt(myKsquared)  # Wille and lie formalism, with Ems this is sqrt(myK), zzzzz
myKfocus = -myK # this should be -myK for focus in x if myK is positive
myKdefocus = myK # take this times 100 to get an interesting envelope # this should be myK for defocus in x if myK is positive
myfQuadL = 0.4 # If FODOF cells set this length to half of mydQuadL
mydQuadL = 0.4

myDriftL = 1.0

mySextuK = 0.6
mySextuL = 0.3

myOctuL = 0.05

approxOffsetQ = max(abs(x))+3*np.std(x)+max(abs(y))+3*np.std(y)
approxOffsetP = max(abs(xp))+3*np.std(xp)+max(abs(yp))+3*np.std(yp)

print 'approxOffsetQ:', approxOffsetQ
print 'approxOffsetP:', approxOffsetP

acceptableError = 1e-9
driftOrder = findOrder(driftham,0,myDriftL,acceptableError,approxOffsetQ,approxOffsetP)
quadOrderf = findOrder(quadham,myKfocus,myfQuadL,acceptableError,approxOffsetQ,approxOffsetP)
quadOrderd = findOrder(quadhamdefocus,myKdefocus,mydQuadL,acceptableError,approxOffsetQ,approxOffsetP)

print 'driftOrder:',driftOrder
print 'quadOrderf:',quadOrderf
print 'quadOrderd:',quadOrderd

sextuOrder = findOrder(sextupoleham,mySextuK,mySextuL,acceptableError,approxOffsetQ,approxOffsetP)

print 'sextuOrder:',sextuOrder


fF = Element('quadfocus', quadham, myKfocus, myfQuadL, quadOrderf)
oO1 = Element('drift', driftham, 0, myDriftL, driftOrder)
dD = Element('quaddefocus', quadhamdefocus, myKdefocus, mydQuadL, quadOrderd)
oO2 = Element('drift', driftham, 0, myDriftL, driftOrder)


## higher order elements
sextupole = Element('sextupole', sextupoleham, mySextuK, mySextuL,sextuOrder)
#octupole = Element('octupole', octupoleham, myK, myOctuL, order)


fodoLattice = list()
nbroffodos = raw_input('Enter number of FODO cells:')
if len(nbroffodos) < 1 : nbroffodos = "50"
nbroffodos = int(nbroffodos)
for i in range(nbroffodos):
    fodoLattice.append(fF)
    #fodoLattice.append(sextupole)
    fodoLattice.append(oO1)
    #fodoLattice.append(oO1)
    #fodoLattice.append(oO1)
    #fodoLattice.append(oO1)
    #fodoLattice.append(dD)
    #fodoLattice.append(sextupole)
    #fodoLattice.append(oO2)
    #fodoLattice.append(oO2)
    #fodoLattice.append(oO2)
    #fodoLattice.append(oO2)

    #fodoLattice.append(fF)

    #if i == 5:
        #fodoLattice.append(sextupole)
    #fodoLattice.append(octupole)

# input from randoms and loads above

#stxoFODO, stxpoFODO = evalLattice(fodoLattice,(stx,stxp))
xoFODO, xpoFODO, yoFODO, ypoFODO, envx, envy = evalLattice(fodoLattice,(x,xp,y,yp)) # Calculate the output values

######## Print phase space
print 'Plotting...'

def plotEverything(xin,xpin,yin,ypin,alpha_x,beta_x,epsilon_x,alpha_y,beta_y,epsilon_y,xo,xpo,yo,ypo,envx,envy):
    plt.figure(0)
    ax1 = plt.subplot2grid((4,3), (0,0))
    plt.plot(xin,xpin,'ro')
    plt.title('Initial value in x')
    plt.xlabel('x')
    plt.ylabel('xp')
    ax2 = plt.subplot2grid((4,3), (0,1))
    plt.plot(yin,ypin,'bo')
    plt.title('Initial value in y')
    plt.xlabel('y')
    plt.ylabel('yp')
    ax3 = plt.subplot2grid((4,3), (0,2))
    plt.plot(xin,yin,'go')
    plt.title('Initial values x and y')
    plt.xlabel('x')
    plt.ylabel('y')
    ax4 = plt.subplot2grid((4,3), (1, 0), colspan=3)
    plt.plot(envx[:,0],envx[:,1],'ro')
    plt.title('Envelope in x by z')
    plt.xlabel('z')
    plt.ylabel('Envelope in x')
    ax5 = plt.subplot2grid((4,3), (2, 0), colspan=3)
    plt.plot(envy[:,0],envy[:,1],'bo')
    plt.title('Envelope in y by z')
    plt.xlabel('z')
    plt.ylabel('Envelope in y')
    ax6 = plt.subplot2grid((4,3), (3, 0))
    plt.plot(xo,xpo,'ro')
    plt.title('Values after FODO in x')
    plt.xlabel('x')
    plt.ylabel('xp')
    ax7 = plt.subplot2grid((4,3), (3, 1))
    plt.plot(yo,ypo,'bo')
    plt.title('Values after FODO in y')
    plt.xlabel('y')
    plt.ylabel('yp')
    ax8 = plt.subplot2grid((4,3), (3, 2))
    plt.plot(xo,yo,'go')
    plt.title('Values after FODO in y and x')
    plt.xlabel('x')
    plt.ylabel('y')

    plt.suptitle("Plots")
    plt.show()

def plotEnvelope(envx,envy):
    plt.figure(0)
    ax4 = plt.subplot2grid((2,3), (0, 0), colspan=3)
    plt.plot(envx[:,0],envx[:,1],'ro')
    plt.title('Envelope in x by z')
    plt.xlabel('z')
    plt.ylabel('Envelope in x')
    ax5 = plt.subplot2grid((2,3), (1, 0), colspan=3)
    plt.plot(envy[:,0],envy[:,1],'bo')
    plt.title('Envelope in y by z')
    plt.xlabel('z')
    plt.ylabel('Envelope in y')

    plt.show()

#plotPhaseSpace(xoFODO, xpoFODO, yoFODO, ypoFODO)
plotEnvelope(envx,envy)

plotEverything(x,xp,y,yp,alpha_x,beta_x,epsilon_x,alpha_y,beta_y,epsilon_y,xoFODO,xpoFODO,yoFODO,ypoFODO,envx,envy)

######## Output
def saveOutput(x,xp,y,yp,alpha_x,beta_x,epsilon_x,alpha_y,beta_y,epsilon_y):
    outputfile = raw_input('Enter file for output particle data:')
    if len(outputfile) < 1 : outputfile = "out.txt"
    outputtwiss = raw_input('Enter file for output twiss data:')
    if len(outputtwiss) < 1 : outputtwiss = "outtwiss.txt"
    
    dt = np.dtype([('x', 'd'), ('xp', 'd'), ('y', 'd'), ('yp', 'd')])
    dtb = np.dtype([('alpha_x', 'd'), ('beta_x', 'd'), ('epsilon_x', 'd'), ('alpha_y', 'd'), ('beta_y', 'd'), ('epsilon_y', 'd')])
    a = np.zeros(len(x), dt)
    b = np.zeros(1,dtb)
    a['x'] = x
    a['xp'] = xp
    a['y'] = y
    a['yp'] = yp

    b['alpha_x'] = alpha_x
    b['beta_x'] = beta_x
    b['epsilon_x'] = epsilon_x
    b['alpha_y'] = alpha_y
    b['beta_y'] = beta_y
    b['epsilon_y'] = epsilon_y
    
    np.savetxt(outputfile, a, '%10s')
    np.savetxt(outputtwiss, b, '%10s')

saveOutput(xoFODO, xpoFODO, yoFODO, ypoFODO,alpha_x,beta_x,epsilon_x,alpha_y,beta_y,epsilon_y)

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

# Finish writing the output "spool" so that all required data is there

# Particles should be stored in the output file when they are done to limit memory usage. This might not be optimal for when space charge forces are introduced

# Merge with Anton's code so that I can easily double check that results are correct. New name is needed as well; lie + Accelerator-?code? = ??????? particle tracker?   

# Make a lattice class with funs: add(), evaluate(x,xp,y,yp), __init__()

##### Comments and stuff not to forget
# 1. Seems like the expansion doesn't hold for k > 1...

# 2. Changed to a minus sign in front of the qy**2 in the Ham to get proper focus/defocus

# 3. Arguments should always be ordered x,xp,y,yp!

# 4. Python uses radians, i.e # print cos(3.14/2).evalf()