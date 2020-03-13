from smoluchowskiMPL import smolMP
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math, numpy

def gamma(phi, r1, rc, rsc):
    if phi != 0:
        z = phi/(1-phi)
        R = rsc/rc
        A1 = R*R*R + 3*R*R + +3*R
        A2 = 3*R*R*R + 4.5*R*R
        A3 = 3*R*R*R
        lng = -1.0*math.log(1-phi) + A1*z + A2*z*z + A3*z*z*z
        return math.exp(lng)
    else:
        return 1

def alpha(phi, r1, rc, rsc):
    if phi != 0:
        R = rsc/rc
        R1 = r1/rsc
        z = phi/(1-phi)
        lna = (2.0/3.0)*R1*R1*R1*(1.5*(R*R+R+1)*z + 4.5*(R*R+R)*z*z + 4.5*R*R*z*z*z)
        return math.exp(lna)
    else:
        return 1

# Initial conditions
M_0 = 0.0
P_0 = 0.0
c_0 = 10

# Parameters
nc = 2 # Critical nucleus
n2 = 2 # Secondary nucleus
phi = 0.0 # Volume fraction of crowders
r1 = 1.0 # Monomer radius
rc = 1.0 # Crowder radius
rsc = 1.0 # Sphero-cylinder radius

alph = alpha(phi, r1, rc, rsc)
gamm = gamma(phi, r1, rc, rsc)

# Rate constants
kp_0 = 2200 * math.pow(10,-6) # Crowderless addition
km = 0.001 # Subtraction
kn_0 =  0.00005 * math.pow(10, -6*(nc-1)) # Crowderless primary nucleation
kfp_0 = 1000 * math.pow(10,-6) # Crowderless coagulation
kfm = 0.001 # Fragmentation
k2_0 = .001 * math.pow(10,-6*n2) # Crowderless secondary nucleation (one-step)

kp = (gamm/alph)*kp_0 # Crowded addition
kfp = (gamm/alph)*kfp_0 # Crowded coagulation
kn = math.pow(gamm/alph, nc-1) * kn_0 # Crowded primary nucleation
k2 = math.pow(gamm, n2)*k2_0

c = [c_0, 0, 0]

dcdt = smolMP(kp, km, kfp, kfm, kn, k2, nc, n2)
t = numpy.linspace(0, 10000, 200)

sol = odeint(dcdt, c, t)

L=[0]
L.extend([a/b for a,b in zip(sol[1:,1],sol[1:,2])])

fig, ax1 = plt.subplots()
ax1.set_xlabel('time')
ax1.plot(t, sol[:,0], 'b')
ax1.plot(t, sol[:,1], 'r')
ax1.plot(t, sol[:,2], 'g')

# ax2 = ax1.twinx()
# ax2.plot(t, L)
fig.tight_layout()
plt.show()
