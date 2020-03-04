import math

def smolMP(kp, km, kpb, kmb, kn, k2, nc, n2):
    def dcdt(c, t):
        P = c[2]
        M = c[1]
        c1 = c[0]
        dP = -kpb*P*P + kmb*(M - (2*nc-1)*P) + kn*math.pow(c1,nc) + k2*math.pow(c1,n2)*M
        dM = 2*(c1*kp - km - 0.5*kmb*(nc*(nc-1)))*P + kn*nc*math.pow(c1,nc) + k2*n2*math.pow(c1,n2)*M
        dc1 = -dM
        return [dc1, dM, dP]
    return dcdt