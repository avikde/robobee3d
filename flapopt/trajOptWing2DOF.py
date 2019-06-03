import autograd.numpy as np
from autograd import jacobian
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
sys.path.append('..')
import planar.FlappingModels as FlappingModels
np.set_printoptions(precision=4, suppress=True, linewidth=200)


# FIXME: need a new one with no body coords
m = FlappingModels.Wing2DOF()

# discrete => do not need solve_ivp

yi = m.y0
dt = 0.01
u = np.array([1,1])
tf = 0.1
tvec = np.arange(0, tf, dt)
yi = np.zeros((len(tvec), len(yi)))
yi[0,:]

for ti in range(1, len(tvec)):
    yi[ti,:] = yi[ti-1,:] + m.dydt(yi[ti-1,:], u) * dt

# display

fig, ax = plt.subplots(2)

ax[0].plot(tvec, yi[:,0:2], '.-')

plt.show()
