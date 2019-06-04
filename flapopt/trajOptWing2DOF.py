import autograd.numpy as np
from autograd import jacobian
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
sys.path.append('..')
import planar.FlappingModels as FlappingModels
from scipy.integrate import solve_ivp
from matplotlib import animation
from matplotlib.collections import PatchCollection
np.set_printoptions(precision=4, suppress=True, linewidth=200)

m = FlappingModels.Wing2DOF()

# discrete => do not need solve_ivp

dt = 0.00001
tf = 0.1
tvec = np.arange(0, tf, dt)
yi = np.zeros((m.nx, len(tvec)))
yi[:,0] = np.array([1e-2, 0, 0, 0])

# params
params = []

def sigmades(t):
    return 10e-3 * np.sin(100 * 2 * np.pi * t)
def closedLoop(t, y):
    # u0 = 1e-3 * np.sin(100 * 2 * np.pi * tvec[ti])
    # pos servoing
    u0 = 1e2 * (sigmades(t) - y[0]) - 1e-2 * y[2]
    # print(u0)
    return m.dydt(y, [u0], params)

for ti in range(1, len(tvec)):
    yi[:,ti] = yi[:,ti-1] + dt * closedLoop(tvec[ti], yi[:,ti-1])


# compare to continuous
sol = solve_ivp(closedLoop, [0,tf], yi[:,0], dense_output=True, t_eval=tvec)

# plots

fig, ax = plt.subplots(3)

# display


# ax[0].plot(tvec, yi[0,:])
ax[0].plot(tvec, sol.y[0,:])
ax[0].plot(tvec, sigmades(tvec))

# ax[1].plot(tvec, yi[1,:])
ax[1].plot(tvec, sol.y[1,:])

# def makeAnim(_ax, _t, _y):
_t = sol.t
_y = sol.y
_ax = ax[2]

p1, = _ax.plot([], [], '.-')
_ax.grid(True)
_ax.set_aspect(1)
_ax.set_ylim(m.cbar * np.array([-2, 2]))

def _init():
    return p1,

def _animate(i):
    wing1 = np.array([_y[0,i], 0])
    c, s = np.cos(_y[1,i]), np.sin(_y[1,i])
    wing2 = wing1 + np.array([[c, -s], [s, c]]) @ np.array([0, -m.cbar])
    p1.set_xdata([wing1[0], wing2[0]])
    p1.set_ydata([wing1[1], wing2[1]])
    return p1,

anim = animation.FuncAnimation(fig, _animate, init_func=_init, frames=len(_t), interval=5000*dt, blit=True)
# makeAnim(ax[2], sol.t, sol.y)

plt.tight_layout()
# print(sol.y.shape)

plt.show()
