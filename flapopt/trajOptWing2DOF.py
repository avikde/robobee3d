import autograd.numpy as np
from autograd import jacobian
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
sys.path.append('..')
import planar.FlappingModels as FlappingModels
from scipy.integrate import solve_ivp
np.set_printoptions(precision=4, suppress=True, linewidth=200)


# FIXME: need a new one with no body coords
m = FlappingModels.Wing2DOF()

# discrete => do not need solve_ivp

dt = 0.001
tf = 0.1
tvec = np.arange(0, tf, dt)
yi = np.zeros((m.nx, len(tvec)))
yi[:,0] = np.array([1e-2, 0, 0, 0])

# params
params = []

for ti in range(1, len(tvec)):
    u0 = 0 * 1e0 * np.sin(100 * 2 * np.pi * tvec[ti])
    u = np.array([u0])
    yi[:,ti] = m.dynamics(yi[:,ti-1], u, dt, params)


# compare to continuous
sol = solve_ivp(lambda t, y: m.dydt(y, u), [0,tf], yi[:,0], dense_output=True, t_eval=tvec)

# display

fig, ax = plt.subplots(2)

ax[0].plot(tvec, yi[0,:], '.-')
ax[0].plot(tvec, sol.y[0,:], '.-')

ax[1].plot(tvec, yi[1,:], '.-')
ax[1].plot(tvec, sol.y[1,:], '.-')

# print(sol.y.shape)

plt.show()
