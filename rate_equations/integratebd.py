import numpy as np
from scipy.integrate import odeint as ode
import matplotlib.pyplot as plt
import math
from beckerdoring import beckerDoringModel as bdm
from beckerdoring import BDinitialConditions as ics

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
kn=3*math.pow(10,-8)
rsc=2.5
r1=2.5
rc=1.9
phi=0
Nmax = 1000
nc=3

c1, gamma, alpha = ics(c0, Nmax, nc, phi, r1=r1, rsc=rsc, rc=rc)
dcdt=bdm(a, b, nc, alpha, gamma, kn=kn)
t=np.linspace(0,20000,100)
sol1 = ode(dcdt, c1, t)
mass1=[]
for s in sol1:
    mass1.append(calculateMass(s,nc))

# phi=0.0375
# c2, gamma, alpha = ics(c0, Nmax, nc, phi, r1=r1, rsc=rsc, rc=rc)
# dcdt=bdm(a, b, nc, alpha, gamma, kn=kn)
# t=np.linspace(0,20000,100)
# sol2 = ode(dcdt, c2, t)
# mass2=[]
# for s in sol2:
#     mass2.append(calculateMass(s,nc))

# phi=0.075
# rc=1.95
# c3, gamma, alpha = ics(c0, Nmax, nc, phi, r1=r1, rsc=rsc, rc=rc)
# dcdt=bdm(a, b, nc, alpha, gamma, kn=kn)
# t=np.linspace(0,20000,100)
# sol3 = ode(dcdt, c3, t)
# mass3=[]
# for s in sol3:
#     mass3.append(calculateMass(s,nc))

# phi=0.15
# rc=2.3
# c4, gamma, alpha = ics(c0, Nmax, nc, phi, r1=r1, rsc=rsc, rc=rc)
# dcdt=bdm(a, b, nc, alpha, gamma, kn=kn)
# t=np.linspace(0,20000,100)
# sol4 = ode(dcdt, c4, t)
# mass4=[]
# for s in sol4:
#     mass4.append(calculateMass(s,nc))
plt.plot(t, [m/c0 for m in mass1], 'b')
# plt.plot(t, mass2, 'p')
# plt.plot(t, mass3, 'r')
# plt.plot(t, mass4, 'k')
plt.show()
