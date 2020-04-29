import matplotlib.pyplot as plt
import numpy as np
import math

def Gamma(phi, R):
    # R = rsc/rc
    if phi != 0:
        z = phi/(1-phi)
        A1 = R*R*R + 3*R*R + +3*R
        A2 = 3*R*R*R + 4.5*R*R
        A3 = 3*R*R*R
        lng = -math.log(1-phi) + A1*z + A2*z*z + A3*z*z*z
        return math.exp(lng)
    else:
        return 1

def Alpha(phi, R, R1):
    # R = rsc/rc
    # R1 = r1/rsc
    if phi != 0:
        z = phi/(1-phi)
        lna = (2.0/3.0)*R1*R1*R1*(1.5*(R*R+R+1)*z + 4.5*(R*R+R)*z*z + 4.5*R*R*z*z*z)
        return math.exp(lna)
    else:
        return 1

Rr = [i for i in np.linspace(0.2,5,100)]
R1r = [i for i in np.linspace(0.5,2,100)]
phir = [i for i in np.linspace(0, 0.5, 100)]
R = 1.2
R1 = 1
phi = 0.1

goa_R = [math.log(Gamma(phi, i)/Alpha(phi,i,R1)) for i in Rr]
goa_R1 = [math.log(Gamma(phi, R)/Alpha(phi, R, i)) for i in R1r]
goa_phi = [math.log(Gamma(i, R)/Alpha(i, R, R1)) for i in phir]
gamma_R = [math.log(Gamma(phi, i)) for i in Rr]
gamma_phi = [math.log(Gamma(i, R)) for i in phir]

print('goa\n2.35')
print(Gamma(0.2, 2.5/2.35)/Alpha(0.1, 2.5/2.35, 2.5/2.35))
print('1.9')
print(Gamma(0.2, 2.5/1.9)/Alpha(0.1, 2.5/1.9, 2.5/1.9))
print('squared\n2.35')
print((Gamma(0.2, 2.5/2.35)/Alpha(0.1, 2.5/2.35, 2.5/2.35))**2)
print('1.9')
print((Gamma(0.2, 2.5/1.9)/Alpha(0.1, 2.5/1.9, 2.5/1.9))**2)
print('gamma cubed\n2.35')
print((Gamma(0.2, 2.5/2.35))**3)
print('1.9')
print((Gamma(0.2, 2.5/1.9))**3)

fig, axs = plt.subplots(3,2)
axs[0, 0].plot(Rr, goa_R)
axs[0, 0].set_title('log g/a vs rsc/rc')
axs[0, 0].text(0.4,0.65, 'phi = {}\nr1/rsc = {}'.format(phi, R1), transform=axs[0,0].transAxes)
axs[0, 1].plot(R1r, goa_R1, 'tab:orange')
axs[0, 1].set_title('log g/a vs r1/rsc')
axs[0, 1].text(0.4,0.65, 'phi = {}\nrsc/rc = {}'.format(phi, R), transform=axs[0,1].transAxes)
axs[1, 0].plot(phir, gamma_phi, 'tab:green')
axs[1, 0].set_title('log g/a vs phi')
axs[1, 0].text(0.4,0.65, 'rsc/rc = {}\nr1/rsc = {}'.format(R, R1), transform=axs[1,0].transAxes)
axs[1, 1].plot(Rr, gamma_R, 'tab:red')
axs[1, 1].set_title('log gamma1 vs rsc/rc')
axs[1, 1].text(0.4,0.8, 'phi = {}'.format(phi), transform=axs[1,1].transAxes)
axs[2, 0].plot(phir, gamma_phi, 'tab:purple')
axs[2, 0].set_title('log gamma1 vs phi')
axs[2, 0].text(0.4,0.8, 'rsc/rc = {}'.format(R), transform=axs[2,0].transAxes)
fig.tight_layout(pad=1.0)
plt.show()

# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# ax1.set_title('log gamma vs r1')
# ax1.text(0.4,0.85, 'phi = {}\nrc = {}\nrsc = {}'.format(phif, rcf, rscf), transform=ax1.transAxes)
# ax1.plot(r1r, gamma_r1r)
# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
# ax2.set_title('log gamma vs rc')
# ax2.text(0.4,0.85, 'phi = {}\nr1 = {}\nrsc = {}'.format(phif, r1f, rscf), transform=ax2.transAxes)
# ax2.plot(rcr, gamma_rcr)
# fig3 = plt.figure()
# ax3 = fig3.add_subplot(111)
# ax3.set_title('log gamma vs rsc')
# ax3.text(0.4,0.85, 'phi = {}\nr1 = {}\nrc = {}'.format(phif, r1f, rcf), transform=ax3.transAxes)
# ax3.plot(rscr, gamma_rscr)
# fig4 = plt.figure()
# ax4 = fig4.add_subplot(111)
# ax4.set_title('log gamma vs phi')
# ax4.text(0.4,0.85, 'r1 = {}\nrc = {}\nrsc = {}'.format(r1f, rcf, rscf), transform=ax4.transAxes)
# ax4.plot(phir, gamma_phir)
# plt.show()


