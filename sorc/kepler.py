from math import *
import numpy as np

"""
M = E - e*sinE
M = mean anomaly, = 2*pi*t/Period
E = eccentric anomaly = angle w.r.t origin (not focus) of actual planet

Solve Kepler's equation for the earth's distance to the sun through time
"""

def kepler(M,e):
#solve for E given M,e
# 2 iterations usually enough to get residual to < 1e-13
  residual = 10.
  iter = 0
  E = M #first guess for a newton's method

  while (iter < 20 and abs(residual) > 1.e-14):
    residual = ( E - M - e*sin(E))
    #newton increment 
    E -= residual / (1. - e*cos(E)) 
    #debug: print("iter ",iter, abs(residual),flush=True)
    iter += 1
    
  return E

def radius(a, e, E):
  r = a*(1-e*cos(E))
  return r


#--------------------------------
# Astronomical numbers:
period = 365.2422
semimajor = 149.598e6 #km
e = 0.01671022         #J2000
semimajor = 1.00000011 #J2000 AU

dt = 0.1
length = int(0.5 + period*10/dt + 1)
time = np.zeros((length))
M = np.zeros((length))
E = np.zeros((length))
de = np.zeros((length))
dr = np.zeros((length))

alpha = 2.*pi/period*dt
for t in range(0,length):
  time[t] = np.double(t)*dt
  M[t] = np.double(t)*alpha
  E[t] = kepler(M[t],e)
  de[t] = M[t] - E[t]
  dr[t] = radius(1., e, E[t])-1.
  #debug: print(t, M[t], E[t], de[t], dr[t])

dr[t] -= dr.mean()

# Begin simple time series analysis
import matplotlib
import matplotlib.pyplot as plt
import scipy

tmp = scipy.fft.fft(dr)
fftdr_amp = np.absolute(tmp)
# Scaling of the fft:
fftdr_amp *= 2./len(dr)

#debug: print("len dr fftdr_amp ",len(dr), len(fftdr_amp))
#debug: print(fftdr_amp.dtype, fftdr_amp[0] )
#debug: print("period = ",period)

#CPD:
fftdr_fre = scipy.fft.fftfreq(len(fftdr_amp), dt) #last = delta_t
#CPY:
#fftdr_fre = scipy.fft.fftfreq(len(fftdr_amp), 1./period) #last = delta_t

fftde_amp = np.absolute(scipy.fft.fft(de))
fftde_fre = scipy.fft.fftfreq(len(fftde_amp), dt) #last = delta_t

#fig,ax = plt.subplots()
#ax.set(xlabel = "cpd")
#ax.plot(fftdr_fre, fftdr_amp)
##amp = np.absolute(fftde_amp)
#ax.plot(fftde_fre, amp)
#ax.grid()
#plt.savefig("ktest.png")

#i=0
#print(i,dr[i], fftdr_amp[i], fftdr_fre[i], 0.)
#for i in range(1, len(fftdr_amp) ):
#  print(i,dr[i], fftdr_amp[i], fftdr_fre[i], 1./fftdr_fre[i] )
#
#debug: print("mean, dr.sum ",dr.mean(), dr.sum() )

#--------------------------------------
for t in range(0,length):
  dr[t] += e*cos(alpha*t)
  #debug: print(t, dr[t])

# Convert to micro-AU
dr *= 1.e6

#CPD:
tmp = scipy.fft.fft(dr)
fftdr_amp = np.absolute(tmp)
fftdr_amp *= 2./len(dr)
fftdr_fre = scipy.fft.fftfreq(len(fftdr_amp), dt) #last = delta_t

#fig,ax = plt.subplots()
#ax.set(title = "after 1 cpy removed")
#ax.plot(fftdr_fre, fftdr_amp)
#ax.grid()
#plt.savefig("ktest2.png")

i=0
print(i,dr[i], fftdr_amp[i], fftdr_fre[i], 0.)
for i in range(1, len(fftdr_amp) ):
  print(i,dr[i], fftdr_amp[i], fftdr_fre[i], 1./fftdr_fre[i] )
