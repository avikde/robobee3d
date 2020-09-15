
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
for i in range(Nt):
    U[i,:] = wl.update(u0, h0, pdotdes)
print(U[0,:])

# convert
Vmean = U[:,0]
uoffs = U[:,1]
udiff = U[:,2]

drv_pch = Vmean * (1 + udiff) * (1 + uoffs)
drv_roll = Vmean * (1 - udiff) * (1 + uoffs)

fig, ax = plt.subplots(2)
ax[0].plot(Vmean, label='Vmean')
ax[0].legend()

ax[1].plot(drv_roll, label='drv_roll')
ax[1].plot(drv_pch, label='drv_pch')
# ax[1].set_ylim((-100,100))
ax[1].legend()

plt.show()

