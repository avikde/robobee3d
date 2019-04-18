'''
Test template MPC ideas on:
- actuated IP as template
- double pendulum, acrobot as anchors
'''
import autograd.numpy as np
from autograd import jacobian
import sys
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits import mplot3d

sys.path.append('..')
from controlutils.py import lqr, mpc
from controlutils.py.models.pendulums import Pendulum, DoublePendulum

np.set_printoptions(precision=4, suppress=True, linewidth=200)

'''
Simulate an IP
'''
pendulum = Pendulum()
doublePendulum = DoublePendulum()

# auto linearization
yup = np.array([np.pi, 0.0])
uup = 0.0
Afun = jacobian(lambda y: pendulum.dynamics(y, uup))
Bfun = jacobian(lambda u: pendulum.dynamics(yup, u))
B = Bfun(uup)[:,np.newaxis]
Aup = Afun(yup)
# LQR
Q = np.eye(2)
R = np.eye(1)
K, P = lqr.lqr(Aup, B, Q, R)

# Simulation
tf = 5.0
dt = 0.01
t_eval = np.arange(0, tf, dt)
y0 = np.array([2,0.0])
sol = solve_ivp(lambda t, y: pendulum.dynamics(y, K @ (yup - y)), [0, tf], y0, dense_output=True, t_eval=t_eval)

# double pendulum with MPC controller
N = 5
wx = np.full(4, 1)
wu = np.full(4, 0.01)
dpmpc = mpc.LTVMPC(doublePendulum, N, wx, wu, verbose=False, polish=False, scaling=0, eps_rel=1e-2, eps_abs=1e-2, kdamping=0)

def mpcDoublePendulum(t,y):
    """Closed-loop dynamics with MPC controller"""
    return doublePendulum.dynamics(y, np.zeros(2))

# simulate
y02 = np.array([1,-1, 0, 0])
sol2 = solve_ivp(mpcDoublePendulum, [0, tf], y02, dense_output=True, t_eval=t_eval)

# visualize value function
def lqrValueFunc(x1,x2):
    # quadratic
    val = P[0,0] * (x1 - yup[0])**2 + P[1,1] * x2**2 + (P[0,1] + P[1,0]) * (x1 - yup[0])*x2
    return val

xx, yy = np.meshgrid(np.linspace(0, 2*np.pi, 30), np.linspace(-10, 10, 30))

# Display -------------------

fig, ax = plt.subplots(3)

ax[0].plot(sol.t, sol.y[0,:], label='sp')
ax[0].plot(sol2.t, sol2.y[0,:], label='dp0')
ax[0].plot(sol2.t, sol2.y[1,:], label='dp1')
ax[0].legend()

ax[1].contourf(xx, yy, lqrValueFunc(xx, yy), cmap='gray_r')
ax[1].plot(yup[0], yup[1], 'r*')

# animation
line1, = ax[2].plot([], [], '.-', lw=2)
line2, = ax[2].plot([], [], '.-', lw=2)
patches = [line1, line2]
ax[2].set_aspect(1)
ax[2].set_xlim((-2,2))
ax[2].set_ylim((-2,2))
ax[2].grid(True)
plt.tight_layout()

def _init():
    line1.set_data([], [])
    line2.set_data([], [])
    return patches

def _animate(i):
    # get the vertices of the pendulum
    p1 = pendulum.kinematics(sol.y[0:1, i])
    line1.set_data([0, p1[0]], [0, p1[1]])

    p1, p2 = doublePendulum.kinematics(sol2.y[0:2, i])
    line2.set_data([0, p1[0], p2[0]], [0, p1[1], p2[1]])
    return patches

anim = animation.FuncAnimation(fig, _animate, init_func=_init, frames=len(sol2.t), interval=1000*dt, blit=True)

plt.show()
