import autograd.numpy as np
from autograd import jacobian, hessian
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from scipy.integrate import solve_ivp
from matplotlib import animation
from matplotlib.collections import PatchCollection
np.set_printoptions(precision=4, suppress=True, linewidth=200)
import osqp
import wingopt

# ---------------------------------------------------
m = wingopt.m

# PARAMETERS -------------
dt = 1e-4 # need for discretization

# Functions ---
_, _, _, xmax = m.limits
def sigmades(t):
    return 0.75 * xmax[0] * np.sin(150 * 2 * np.pi * t)
def controller(t, y):
    return 1e2 * (sigmades(t) - y[0]) - 1e-2 * y[2]
def strokePosControlVF(t, y):
    # u0 = 1e-3 * np.sin(100 * 2 * np.pi * tvec[ti])
    # pos servoing
    u0 = controller(t, y)
    # print(u0)
    return m.dydt(y, [u0], wingopt.params)


def Jcostsol(solt, soly, params):
    Nt = len(solt)
    Jcosti = 0
    for i in range(Nt):
        Jcosti += wingopt.Jobjinst(soly[:,i], controller(solt[i], soly[:,i]), params)
    return Jcosti

"""Get initial OL trajectory with a small-timestep simulation ----------------------
"""
# Sim params
tf = 0.1
tvec = np.arange(0, tf, dt)
y0 = np.array([1e-2, 0, 0, 0])
# Sim of an openloop controller
sol = solve_ivp(strokePosControlVF, [0,tf], y0, dense_output=True, t_eval=tvec)

# Decimate the rate to get a starting traj for optimization with fewer knot points
# At this point there is a "nominal trajectory" in sol.y
yu0 = sol.y.copy()
# yu0 = yi.copy()
utest = np.zeros(yu0.shape[1])
for i in range(yu0.shape[1]):
    utest[i] = controller(sol.t[i], sol.y[:,i])
yu0 = np.vstack((yu0, utest))
# This range and decimation depends on the OL traj. TODO: Should automate finding 1 cycle
olTraj = (yu0.T)[170:238:3,:]
olTrajt = sol.t[170:238:3]
olTrajdt = np.mean(np.diff(olTrajt))
m.dt = olTrajdt  # for the discretized dynamics

# Get the decimated trajectory to be feasible by iterating the linearized dynamics with the inputs
yi2 = olTraj.copy()
for ti in range(1, len(olTrajt)):
    ui = olTraj[ti-1, 4:] # OL input from previous traj
    # Linearized
    A, B, c = m.getLinearDynamics(olTraj[ti-1, :4], ui, wingopt.params)
    olTraj[ti, :4] = A @ olTraj[ti-1, :4] + B @ ui + c
    # Nonlinear
    yi2[ti, :4] = yi2[ti-1, :4] + olTrajdt * m.dydt(yi2[ti-1, :4], ui, wingopt.params)
# plotTrajs(olTraj, yi2, yilin)


# Use the class above to step the QP
Nknot = olTraj.shape[0]  # number of knot points in this case
nx = 4
nu = 1
wx = np.ones(nx) * 1e-6
wu = np.ones(nu) * 1e3
kdampx = 1e-5 * np.ones(4)
kdampu = np.zeros(1)
# Must be 1 smaller to have the correct number of xi
wqp = wingopt.WingQP(m, Nknot-1, wx, wu, kdampx, kdampu, verbose=False, eps_rel=1e-4, eps_abs=1e-4, max_iter=10000)
wqp.trajt = olTrajt
# Test warm start
# wqp.ltvsys.prob.warm_start(x=dirTranForm(olTraj, Nknot, 4, 1))
traj2 = wqp.update(olTraj)
traj3 = wqp.update(traj2)
# TEST: reduce wu 
# wqp.ltvsys.qof.wu = np.ones(nu) * 1e2
traj4 = wqp.update(traj3)
# wqp.debugConstraintViol(olTraj, wqp.dirtranx)

# Create shorter timestep sims
ctstrajs = []
tvec = np.arange(0, olTrajt[-1] - olTrajt[0], dt)
tf = tvec[-1]
y0 = olTraj[0,:4]

def knotPointControl(t, y, traj):
    # dumb TODO: interpolate
    u0 = 0 # which knot point input to use
    for i in range(len(olTrajt)):
        if olTrajt[0] + t > olTrajt[i]:
            u0 = traj[i,4]
    return m.dydt(y, [u0], wingopt.params)

for opttraj in [olTraj, traj2, traj3, traj4]:
    # Sim of an openloop controller
    sol = solve_ivp(lambda t, y: knotPointControl(t, y, opttraj), [0, tf], y0, dense_output=True, t_eval=tvec)
    ctstrajs.append(sol.y.T)

# Optim wrt params ----

JT = lambda T : wingopt.Jcost_dirtran(wqp.dirtranx, Nknot, [wingopt.params[0], T])
DJT = jacobian(JT)

Ttest = np.arange(0.5, 1.5, 0.1)
plt.plot(Ttest, [JT(Ti) for Ti in Ttest])

# Display -------------

# wqp.plotTrajs(olTraj, traj2, traj3, traj4)
# wingopt.trajAnim(tvec, ctstrajs)
plt.show()
sys.exit(0)

# ---------------------- OLD gradient descent ------------

Jgrad = jacobian(lambda yu : wingopt.Jcosttraj(yu, params))
# gradient wrt params
Jgradp = jacobian(lambda p : wingopt.Jcosttraj(yu0, p))

# Gradient descent
print(Jcosttraj(yu0, params))
yu1 = yu0.copy()
for i in range(5):
    g1 = Jgrad(yu1)
    # print(g1[:20])
    # # Stupid "line search" for step size
    # for i in range(-10,10):
    #     print(i, Jcosttraj(yu0 - np.power(10.0,i) * g1, params))
    # sys.exit(0)
    yu1 -= 1e1 * g1
    print(i, Jcosttraj(yu1, params))

params1 = params.copy()
for i in range(5):
    g1 = Jgradp(params1)
    # # Stupid "line search" for step size
    # for i in range(-10,10):
    #     print(i, Jcosttraj(yu0, params - np.power(10.0,i) * g1))
    # sys.exit(0)
    params1 -= 1e-6 * g1
    print(i, Jcosttraj(yu0, params1))

# --------------------------------------------------------

# plots

fig, ax = plt.subplots(2)

# display

# # ax[0].plot(tvec, yi[0,:])
# ax[0].plot(tvec, yu0[0,:], label='0')
# ax[0].plot(tvec, yu1[0,:], label='1')
# ax[0].plot(tvec, sigmades(tvec), 'k--')
# ax[0].legend()

# ax[0].plot(tvec, yi[1,:])
ax[0].plot(tvec, yu0[1,:], label='0')
ax[0].plot(tvec, yu1[1,:], label='1')
ax[0].legend()

# def makeAnim(_ax, _t, _y):

_ax = ax[-1]

p1, = _ax.plot([], [], 'b.-', linewidth=4)
paero, = _ax.plot([], [], 'r', linewidth=2)
p2, = _ax.plot([], [], 'b.-', linewidth=4)
paero2, = _ax.plot([], [], 'r', linewidth=2)
p3, = _ax.plot([], [], 'b.-', linewidth=4)
paero3, = _ax.plot([], [], 'r', linewidth=2)
_ax.grid(True)
_ax.set_aspect(1)
_ax.set_ylim(m.cbar * np.array([-4, 2]))

def _init():
    return p1, paero, p2, paero2, p3, paero3, 

def _animate(i):
    wingopt.flapVisUpdate(yu0[:,i], [-0.02,0], params, p1, paero)
    wingopt.flapVisUpdate(yu1[:,i], [0,0], params, p2, paero2)
    wingopt.flapVisUpdate(yu0[:,i], [0.02,0], params1, p3, paero3)
    return p1, paero, p2, paero2, p3, paero3, 

anim = animation.FuncAnimation(fig, _animate, init_func=_init, frames=yu0.shape[1], interval=2e5*dt, blit=True)
if False:
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
    import time
    timestamp = time.strftime('%Y%m%d%H%M%S', time.localtime())
    anim.save('trajOptWing2DOF'+timestamp+'.mp4', writer=writer)

plt.tight_layout()
# print(sol.y.shape)

plt.show()
