
from wlqppy import WLController
import numpy as np
import matplotlib.pyplot as plt

wl = WLController()

pdotdes = [0, 0, 10, 0, 0, 0]
mb = 100
g = 9.81e-3
h0 = [0, 0, mb * g, 0, 0, 0]
u0 = [140.0, 0., 0., 0.]

Nt = 500
U = np.zeros((Nt,4))
U[0,:] = u0
for i in range(1,Nt):
    if i > 100:
        pdotdes[2] = 20
    if i > 200:
        pdotdes[2] = 30
    if i > 300:
        pdotdes[2] = 40
    if i > 400:
        pdotdes[2] = 50
    U[i,:] = wl.update(U[i-1,:], h0, pdotdes)
print(U[0,:])

# convert
Vmean = U[:,0]
uoffs = U[:,1]
udiff = U[:,2]

vleft = Vmean * (1 + udiff) * (1 + uoffs)
vright = Vmean * (1 - udiff) * (1 + uoffs)

drv_pch = Vmean * uoffs
drv_amp = np.abs((vleft-vright))/2 + np.minimum(vleft,vright)
drv_roll = (vleft-vright)/4; 
drv_bias = np.maximum(vleft,vright) + 2* np.abs(drv_pch) # voltage

# print(udiff, uoffs, drv_pch)

fig, ax = plt.subplots(2)
ax[0].plot(Vmean, label='Vmean')
ax[0].plot(drv_bias, label='drv_bias')
ax[0].plot(drv_amp, label='drv_amp')
ax[0].set_ylim((100,220))
ax[0].legend()

ax[1].plot(drv_roll, label='drv_roll')
ax[1].plot(drv_pch, label='drv_pch')
ax[1].axhline(0, alpha=0.3)
ax[1].set_ylim((-30,30))
ax[1].legend()

plt.show()

