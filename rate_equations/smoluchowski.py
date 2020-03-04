import math

def InitialConditions(c0, Nmax, nc=2, phi=0, **kwargs):
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

def delta(r, s):
    if r==s:
        return 1
    else:
        return 0

def bdJ(c, r, a, b, gamma=1.0, alpha=1.0):
    goa=gamma/alpha
    try:
        j = goa*a*c[0]*(c[r-1]-c[r]) + b*(c[r+1]-c[r])
    except IndexError:
        j = goa*a*c[0]*(c[r-1]-c[r]) - b*c[r]
    return j

def merJ(c, rps_idx, ka, nc=2, gamma=1.0, alpha=1.0):
    goa=gamma/alpha
    j=0
    rps = rps_idx + nc - 1
    # loss due to merging
    ## POSSIBLE WEIRDNESS WITH INDEXING HERE
    for s,cs in enumerate(c[1:-rps],1):
        j = j - goa*ka*c[rps_idx]*cs*(1+delta(rps_idx,s))
    # gain due to merging
    rps_idx_half = (rps_idx+nc-1)/2.0 - nc + 1
    if rps_idx_half > 0:
        for s_idx,cs in enumerate(c[1:math.floor(rps_idx_half)],1):
            s = s_idx + nc - 1
            r = rps - s
            r_idx = r - nc + 1
            j = j + goa*ka*c[r_idx]*cs*(1.0-0.5*delta(r,s))
    return j

def breJ(c, r_idx, kb, nc=2):
    j=0
    r = r_idx + nc - 1
    # Loss due to breaking
    if r > nc:
        # break into any two mers bigger than monomer
        j = j - (r - 3)*kb*c[r_idx]
    # gain due to breaking
    for rps_idx, crps in enumerate(c[r_idx+2:],r_idx+2):
        rps = rps_idx + nc + 1
        if r > rps/2.0:
            j = j + kb*crps
        else:
            j = j + 2*kb*crps
    return j

def c1j(c, a, b, kb, kn, nc, gamma, alpha):
    goa=gamma/alpha
    goanc=math.pow(goa, nc-1)
    cp = sum(c[1:-1])
    crGTnc = sum(c[2:])
    dc1dtbd = b*(crGTnc + nc*c[1]) - nc*goanc*kn*pow(c[0],nc) - goa*a*c[0]*cp
    j = 0
    for rps_idx, crps in enumerate(c[1:],1):
        rps = rps_idx + nc - 1
        if rps > 3:
            if rps < nc+2:
                j = j + rps*(rps-3)*kb*crps
            else:
                for s in range(2, nc):
                    j = j + 2*s*kb*crps
    dc1dt = dc1dtbd + j
    return dc1dt

def cnj(c, a, b, ka, kb, kn, nc, gamma, alpha):
    goa = gamma/alpha
    goanc = math.pow(goa, nc-1)
    j=0
    # Loss due to breaking
    if nc > 3:
        # break into any two mers bigger than monomer
        j = j - (nc - 3)*kb*c[1]
    # loss due to merging
    ## POSSIBLE WEIRDNESS WITH INDEXING HERE
    for s_idx,cs in enumerate(c[1:-nc],1):
        s = s_idx + nc - 1
        j = j - goa*ka*c[1]*cs*(1+delta(nc,s))
    
    # gain due to breaking
    for rps_idx, crps in enumerate(c[3:],3):
        rps = rps_idx + nc + 1
        if nc > rps/2.0:
            j = j + kb*crps
        else:
            j = j + 2*kb*crps
    
    # gain from regular nucleation
    j = j + goanc*kn*math.pow(c[0],nc)
    # gain from monomer subtraction
    j = j + b*c[2]
    # loss from addition
    j = j - goa*a*c[0]*c[1]
    # loss from subtraction
    j = j - b*c[1]
    return j

def J(c, a, b, ka, kb, nc, idx, gamma, alpha):
    return bdJ(c, idx, a, b, gamma, alpha) + merJ(c, idx, ka, nc, gamma, alpha) + breJ(c, idx, kb, nc)

def smoluchowskiModel(a, b, nc=2, gamma=1, alpha=1, **kwargs):
    if 'kn' not in kwargs:
        kn = a/nc
    else:
        kn = kwargs['kn']
    if 'ka' not in kwargs:
        ka = a
    else:
        ka = kwargs['ka']
    if 'kb' not in kwargs:
        kb = b
    else:
        kb = kwargs['kb']
    def DCDT(c, t):
        dc1dt = c1j(c, a, b, kb, kn, nc, gamma, alpha)
        dcndt = cnj(c, a, b, ka, kb, kn, nc, gamma, alpha)
        dcdt=[dc1dt,dcndt]
        dcrdt=[J(c, a, b, ka, kb, nc, r, gamma, alpha) for r, _ in enumerate(c[2:],2)]
        dcdt.extend(dcrdt)
        return dcdt
    return DCDT