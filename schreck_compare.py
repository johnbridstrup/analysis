from rate_equations.smoluchowskiMPL import smolMP
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math
import numpy as np
from utilities import getData

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
c_0 = 5

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
kp_0 = 455.5 # Crowderless addition
km = 50 # Subtraction
kn_0 = .0000102 # Crowderless primary nucleation
kfp_0 = 0 # Crowderless coagulation
kfm = 0 # Fragmentation
k2_0 = 0 # Crowderless secondary nucleation (one-step)

kp = (gamm/alph)*kp_0 # Crowded addition
kfp = (gamm/alph)*kfp_0 # Crowded coagulation
kn = math.pow(gamm/alph, nc-1) * kn_0 # Crowded primary nucleation
k2 = math.pow(gamm, n2)*k2_0

c = [c_0, 0, 0]

dcdt = smolMP(kp, km, kfp, kfm, kn, k2, nc, n2)
t = np.linspace(0, 10, 400)

sol = odeint(dcdt, c, t)

L=[0]
L.extend([a/b for a,b in zip(sol[1:,1],sol[1:,2])])

fig, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
ax1.set_xlabel('time')
# ax1.plot(t, sol[:,0], 'b')
ax1.plot(t, [i/c_0 for i in sol[:,1]], color='k',linestyle='--')
ax2.plot(t,[i/j for i,j in zip(sol[:,1],sol[:,2])], color='k', linestyle='--')
# plt.figure()
# plt.plot(t, [i/j for i,j in zip(sol[:,1],sol[:,2])], 'g')

# ax2 = ax1.twinx()
# ax2.plot(t, L)


folder_path = 'data/schreck_comp'
data = getData(folder_path, 'N')

M = [i/30000 for i in data['30000']['M']]
M60 = [i/60000 for i in data['60000']['M']]
M50 = [i/50000 for i in data['50000']['M']]
M40 = [i/40000 for i in data['40000']['M']]
M20 = [i/20000 for i in data['20000']['M']]
M10 = [i/10000 for i in data['10000']['M']]
L = data['30000']['L']
L10 = data['10000']['L']
L20 = data['20000']['L']
L40 = data['40000']['L']
L50 = data['50000']['L']
L60 = data['60000']['L']
st = data['30000']['t']
# ax1.plot(st, M, 'b')
# ax1.plot(st, M10, 'g')
# ax1.plot(st, M20, 'k')
ax1.plot(st, M60, 'b')
# ax1.plot(st, M10, 'g')
# ax2.plot(st, L, 'b')
ax2.plot(st, L10, 'g')
# ax2.plot(st, L20, 'k')
# ax2.plot(st, L40, 'm')
ax2.plot(st, L60, 'b')
fig.tight_layout()
plt.show()
