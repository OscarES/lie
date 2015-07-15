import numpy as np

##### Diff mellan mina o antons resultat, datan ska va sparad som x column xp column
antonsdata = raw_input('Enter anton\'s datafile name:')
if len(antonsdata) < 1 : antonsdata = "out.txt"
oscarsdata = raw_input('Enter oscar\'s datafile name:')
if len(oscarsdata) < 1 : oscarsdata = "out.txt"

antonx, antonxp = np.loadtxt(antonsdata,unpack = True)
oscarx, oscarxp = np.loadtxt(oscarsdata,unpack = True)

diffx = antonx - oscarx
diffxp = antonxp - oscarxp

print 'diffx:',diffx
print 'diffxp:',diffxp

stdx = sum(diffx**2) - sum(diffx)**2
stdxp = sum(diffxp**2) - sum(diffxp)**2

print 'stdx:',stdx
print 'stdxp:',stdxp

## TODO: make the program do it for 2D and follow the format (column vectors):
# x, xp, y, yp, alpha, beta, epsilon