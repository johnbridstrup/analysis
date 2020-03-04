import numpy as np
import math

def BDinitialConditions(c0, Nmax, nc=2, phi=0, **kwargs):
    if (Nmax<2):
        # some kind of error here
        return 0
    if (nc<2):
        return 0
    if (Nmax<nc):
        return 0
    size = Nmax + 2 - nc
    c=[0 for i in range(size)]
    c[0] = c0
    if phi==0:
        alpha=1
        gamma=1
    else:
        r1=kwargs['r1']
        rsc=kwargs['rsc']
        rc=kwargs['rc']
        R=rsc/rc
        R1 = r1/rsc 
        z=phi/(1-phi)
        A1=R*R*R+3*R*R+3*R
        A2=3*R*R*R+4.5*R*R
        A3=3*R*R*R
        lng = math.log(1-phi) + A1*z + A2*z*z + A3*z*z*z
        lna = (2.0/3.0)*R1*R1*R1*(1.5*(R*R+R+1)*z + 4.5*(R*R+R)*z*z + 4.5*R*R*z*z*z)
        gamma = math.exp(lng)
        alpha = math.exp(lna)
    return (c, gamma, alpha)

def J(c, r, a, b, gamma=1.0, alpha=1.0):
    goa=gamma/alpha
    try:
        j = goa*a*c[0]*(c[r-1]-c[r]) + b*(c[r+1]-c[r])
    except IndexError:
        j = goa*a*c[0]*(c[r-1]-c[r]) - b*c[r]
    return j

def beckerDoringModel(a, b, nc=2, alpha=1, gamma=1, **kwargs ):
    if 'kn' not in kwargs:
        kn = a/nc
    else:
        kn = kwargs['kn']
    goa=gamma/alpha
    goanc = math.pow(goa, nc-1)
    def DCDT(c, t):
        cp = sum(c[1:-1])
        crGTnc = sum(c[2:])
        dc1dt = b*(crGTnc + nc*c[1]) - nc*goanc*kn*pow(c[0],nc) - goa*a*c[0]*cp
        dcndt = goanc*kn*pow(c[0],nc) + b*(c[2]-c[1]) - goa*a*c[0]*c[1]
        dcdt = [dc1dt, dcndt]
        dcrdt = [J(c, r, a, b, gamma, alpha) for r, _ in enumerate(c[2:], 2)]
        dcdt.extend(dcrdt)
        return dcdt  
    return DCDT