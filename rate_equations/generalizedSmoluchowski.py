import math

def ICs(c0, Nmax):
    c = [0 for i in range(Nmax)]
    c[0] = c0
    return c

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

def merJ(c, r, s, ka, nc=2, gamma=1.0, alpha=1.0):
    goa=gamma/alpha
    j=0
    # check that it is possible just in case
    if r+s<len(c):
        j = goa*ka*c[r-1]*c[s-1]
    return j

def breJ(c, r, s, kb, nc=2):
    j=0
    if r>=nc and s>=nc:
        j = kb*c[r+s-1]
    return j