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

# Models
pendulum = {'model': Pendulum()}
pendulum2 = {'model': DoublePendulum()}
acrobot = {'model': DoublePendulum()}

"""Single ----------------------------
"""
yup = np.array([np.pi, 0.0])
uup = np.array([0.0])
A, B, c = pendulum['model'].autoLin(yup, uup)
# LQR
K1, pendulum['P'] = lqr.lqr(A, B, Q=np.eye(2), R=0.01*np.eye(1))

# Simulation
tf = 5.0
dt = 0.01
t_eval = np.arange(0, tf, dt)
y0 = np.array([2,0.0])
pendulum['sol'] = solve_ivp(lambda t, y: pendulum['model'].dynamics(y, K1 @ (yup - y)), [0, tf], y0, dense_output=True, t_eval=t_eval)

"""Double pendulum ----------------------------
"""
yup2 = np.array([np.pi, 0, 0, 0])
uup2 = np.zeros(2)
# LQR
A, B, c = pendulum2['model'].autoLin(yup2, uup2)
K2, pendulum2['P'] = lqr.lqr(A, B, Q=np.eye(4), R=0.01*np.eye(2))
# Simulate
y02 = np.array([0.2, 0.3, 0, 0])
pendulum2['sol'] = solve_ivp(lambda t, y: pendulum2['model'].dynamics(y, K2 @ (yup2 - y)), [0, tf], y02, dense_output=True, t_eval=t_eval)

# MPC
N = 10
wx = np.full(4, 100000)
wu = np.full(2, 0.01)
# add model limits here by munging
umax = np.full(2, np.inf)
xmax = np.full(4, np.inf)
pendulum2['model'].limits = -umax, umax, -xmax, xmax
# Instantiate MPC
pendulum2['model'].dt = 0.2
dpmpc = mpc.LTVMPC(pendulum2['model'], N, wx, wu, verbose=False, polish=False, scaling=0, eps_rel=1e-2, eps_abs=1e-2, kdamping=0)

# simulate
# need to have discrete changes to the MPC
def discretizationEvent(t, y):
    return t - tev[-1] - pendulum2['model'].dt
discretizationEvent.direction = 1
discretizationEvent.terminal = True

# initial states
tev = np.zeros(1)
uev = np.zeros((2, 1))
yev = np.zeros((4, 1))
yev[:, 0] = y02
uev[:, 0] = dpmpc.update(y02, yup2)
pendulum2['tt'] = np.zeros(0)
pendulum2['yy'] = np.zeros((4, 0))

while True:
    sol = solve_ivp(lambda t, y: pendulum2['model'].dynamics(y, uev[:, -1]), [tev[-1], tf], yev[:, -1], dense_output=True, events=discretizationEvent, t_eval=np.arange(tev[-1], tf - 1e-3, dt))

    if sol.status == 1:
        # events
        tev = np.hstack((tev, sol.t_events[0][0]))
        yev = np.hstack((yev, sol.sol(tev[-1])[:, np.newaxis]))
        uev = np.hstack((uev, dpmpc.update(yev[:, -1], yup2)[:, np.newaxis]))
        # continuous for plotting
        pendulum2['tt'] = np.hstack((pendulum2['tt'], sol.t))
        pendulum2['yy'] = np.hstack((pendulum2['yy'], sol.y))
        # should we go again?
        if tf - tev[-1] < dt:
            break

    if sol.status == 0:
        break

"""Acrobot ----------------------------
"""
# LQR
A, B, c = acrobot['model'].autoLin(yup2, np.zeros(1)) # only 1 input
K3, acrobot['P'] = lqr.lqr(A, B, Q=np.eye(4), R=0.01*np.eye(1))
# simulate
y03 = np.array([3,0,0,0])
# acrobot['sol'] = solve_ivp(lambda t, y: acrobot['model'].dynamics(y, K3 @ (yup2 - y)), [0, tf], y03, dense_output=True, t_eval=t_eval)
print(K3, K3 @ (yup2 - y03))
acrobot['sol'] = solve_ivp(lambda t, y: acrobot['model'].dynamics(y, np.zeros(2)), [0, tf], y03, dense_output=True, t_eval=t_eval)

# ---

# visualize value function
def lqrValueFunc(x1, x2, P):
    # quadratic
    val = P[0,0] * (x1 - yup[0])**2 + P[1,1] * x2**2 + (P[0,1] + P[1,0]) * (x1 - yup[0])*x2
    return val

xx, yy = np.meshgrid(np.linspace(0, 2*np.pi, 30), np.linspace(-10, 10, 30))

# Display -------------------

fig, ax = plt.subplots(3)

print(pendulum2['yy'].shape)

ax[0].plot(pendulum['sol'].t, pendulum['sol'].y[0, :], label='sp')
ax[0].plot(pendulum2['sol'].t, pendulum2['sol'].y[:2, :].T, label='dp')
ax[0].plot(pendulum2['tt'], pendulum2['yy'][:2, :].T, label='dpmpc')
ax[0].legend()

ax[1].contourf(xx, yy, lqrValueFunc(xx, yy, pendulum['P']), cmap='gray_r')
ax[1].plot(yup[0], yup[1], 'r*')

# animation
line1, = ax[2].plot([], [], '.-', lw=2, label='sp')
line2, = ax[2].plot([], [], '.-', lw=2, label='dp')
line3, = ax[2].plot([], [], '.-', lw=2, label='dpmpc')
line4, = ax[2].plot([], [], '.-', lw=2, label='acro')
patches = [line1, line2, line3]
ax[2].set_aspect(1)
ax[2].set_xlim((-2,2))
ax[2].set_ylim((-2,2))
ax[2].grid(True)
ax[2].legend(bbox_to_anchor=(2,1))
plt.tight_layout()

def _init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    return patches

def _animate(i):
    # get the vertices of the pendulum
    if i < len(pendulum['sol'].t):
        p1 = pendulum['model'].kinematics(pendulum['sol'].y[0:1, i])
        line1.set_data([0, p1[0]], [0, p1[1]])

    if i < len(pendulum2['sol'].t):
        p1, p2 = pendulum2['model'].kinematics(pendulum2['sol'].y[0:2, i])
        line2.set_data([0, p1[0], p2[0]], [0, p1[1], p2[1]])

    if i < len(pendulum2['tt']):
        p1, p2 = pendulum2['model'].kinematics(pendulum2['yy'][0:2, i])
        line3.set_data([0, p1[0], p2[0]], [0, p1[1], p2[1]])

    if i < len(acrobot['sol'].t):
        p1, p2 = acrobot['model'].kinematics(acrobot['sol'].y[0:2, i])
        line4.set_data([0, p1[0], p2[0]], [0, p1[1], p2[1]])
    return patches

anim = animation.FuncAnimation(fig, _animate, init_func=_init, frames=len(pendulum2['sol'].t), interval=1000*dt, blit=True)

plt.show()
