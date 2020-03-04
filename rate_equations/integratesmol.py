import numpy as np
from scipy.integrate import odeint as ode
import matplotlib.pyplot as plt
import math
from smoluchowski import smoluchowskiModel as smol
from smoluchowski import InitialConditions as ics

# Moments
def calculateMass(s, nc):
    M=0
    for r, c in enumerate(s[1:],nc):
        M=M+c*r
    return M

# Initial Conditions
c0 = 5.0
a=0.0231
b=0.03
kn=3*math.pow(10,-3)
rsc=2.5
r1=2.5
rc=1.9
phi=0
Nmax = 100
nc=3

c1, gamma, alpha = ics(c0, Nmax, nc, phi, r1=r1, rsc=rsc, rc=rc)
dcdt=smol(a, b, nc, alpha, gamma, kn=kn)
t=np.linspace(0,2000,100)
sol1 = ode(dcdt, c1, t)
mass1=[]
for s in sol1:
    mass1.append(calculateMass(s,nc))
plt.plot(t,mass1)
plt.plot(t, sol1[:,0])
plt.show()