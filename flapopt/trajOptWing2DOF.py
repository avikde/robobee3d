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
import osqp
import scipy.sparse as sparse
import controlutils.py.mpc as mpc
import controlutils.py.ltvsystem as ltvsystem

# ---------------------------------------------------

m = FlappingModels.Wing2DOF()

# PARAMETERS -------------
dt = 1e-4 # need for discretization
# params: [cbar, ]
params = np.array([5e-3])

# Functions ---

def sigmades(t):
    return 15e-3 * np.sin(150 * 2 * np.pi * t)
def controller(t, y):
    return 1e2 * (sigmades(t) - y[0]) - 1e-2 * y[2]
def closedLoop(t, y):
    # u0 = 1e-3 * np.sin(100 * 2 * np.pi * tvec[ti])
    # pos servoing
    u0 = controller(t, y)
    # print(u0)
    return m.dydt(y, [u0], params)

def Jobjinst(y, u, params):
    _, Faero = m.aero(y, u, params)
    return -Faero[1] # minimization

def Jcostinst_dynpenalty(ynext, y, u, params):
    '''error on dynamics for penalty method'''
    dynErr = ynext - (y + m.dydt(y, u, params) * dt)
    return 1/2 * dynErr.T @ dynErr

# FIXME: need to find y that make one cycle
# as a first pass just average over the whole time

def Jcosttraj(yu, params):
    '''this is over a traj. yu = (nx+nu,Nt)-shaped'''
    Nt = yu.shape[1]
    c = 0
    PENALTY = 1e-6
    for i in range(Nt-1):
        c += Jobjinst(yu[:m.nx,i], yu[m.nx:,i], params) + PENALTY * Jcostinst_dynpenalty(yu[:m.nx,i+1], yu[:m.nx,i], yu[m.nx:,i], params)
    # TODO: any final cost?
    c += Jobjinst(yu[:m.nx,-1], yu[m.nx:,-1], params)
    return c

def Jcostsol(solt, soly, params):
    Nt = len(solt)
    Jcosti = 0
    for i in range(Nt):
        Jcosti += Jobjinst(soly[:,i], controller(solt[i], soly[:,i]), params)
    return Jcosti

Jgrad = jacobian(lambda yu : Jcosttraj(yu, params))

# --- SIMULATIONS DISCRETIZED AND CTS ---
# Sim params
tf = 0.1
tvec = np.arange(0, tf, dt)
yi = np.zeros((m.nx, len(tvec)))
yi[:,0] = np.array([1e-2, 0, 0, 0])
yilin = yi.copy() # linearized for comparison TODO: remove

for ti in range(1, len(tvec)):
    # Nonlinear
    yi[:,ti] = yi[:,ti-1] + dt * closedLoop(tvec[ti], yi[:,ti-1])
    # Linearized
    ui = np.array([controller(tvec[ti], yilin[:,ti-1])])
    A, B, c = m.getLinearDynamics(yilin[:,ti-1], ui)
    # yilin[:,ti] = A @ yilin[:,ti-1] + B @ ui + c

# compare to continuous
sol = solve_ivp(closedLoop, [0,tf], yi[:,0], dense_output=True, t_eval=tvec)

print('Avg cost =', Jcostsol(sol.t, sol.y, params))

# --------
# At this point there is a "nominal trajectory" in sol.y
# FIXME: update this traj and have one cycle
nominalTraj = (sol.y.copy().T)[0:100,:]

# Create "MPC" object which will be used for SQP
wx = np.array([1,1,1,1])
wu = np.array([1])
peps = 1e-2
Nknot = nominalTraj.shape[0]  # number of knot points in this case
ltvqp = mpc.LTVMPC(m, Nknot, wx, wu, verbose=True, polish=False, scaling=0, eps_rel=peps, eps_abs=peps, max_iter=100, kdamping=0)
# x0 must be a (N,nx) trajectory
xr = np.zeros(m.nx)  # FIXME: does not make sense
ctrl = ltvqp.update(x0=nominalTraj, xr=xr, trajMode=ltvsystem.GIVEN_POINT_OR_TRAJ)

# # Test the linearized dynamics
# strokeEnd = 1e-3
# ytest = np.array([-strokeEnd, 0, 0, 0])
# umax = 1e-3
# utest = np.array([umax])
# A, B, c = m.getLinearDynamics(ytest, utest)
# print(A, B, c)

# TODO: simulate forward with nonlin and lin dynamics for a short time and compare the trajectories


sys.exit(0)

# Trajectory to begin gradient descent from ------------
yu0 = sol.y.copy()
# yu0 = yi.copy()
utest = np.zeros(yu0.shape[1])
for i in range(yu0.shape[1]):
    utest[i] = controller(sol.t[i], sol.y[:,i])
yu0 = np.vstack((yu0, utest))
# gradient wrt params
Jgradp = jacobian(lambda p : Jcosttraj(yu0, p))

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

def flapkin(yui, xyoff, _params):
    # wing extents
    wing1 = np.array([yui[0], 0]) + np.asarray(xyoff)
    c, s = np.cos(yui[1]), np.sin(yui[1])
    wing2 = wing1 + np.array([[c, -s], [s, c]]) @ np.array([0, -2*_params[0]])
    # aero arrow extents
    _, Faero = m.aero(yui[:m.nx], yui[m.nx:], _params)
    pcop = (wing1 + wing2)/2
    aeroEnd = pcop + 0.3 * Faero
    return wing1, wing2, pcop, aeroEnd

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
    wing1, wing2, pcop, aeroEnd = flapkin(yu0[:,i], [-0.02,0], params)
    p1.set_xdata([wing1[0], wing2[0]])
    p1.set_ydata([wing1[1], wing2[1]])
    paero.set_xdata([pcop[0], aeroEnd[0]])
    paero.set_ydata([pcop[1], aeroEnd[1]])
    wing1, wing2, pcop, aeroEnd = flapkin(yu1[:,i], [0.0,0], params)
    p2.set_xdata([wing1[0], wing2[0]])
    p2.set_ydata([wing1[1], wing2[1]])
    paero2.set_xdata([pcop[0], aeroEnd[0]])
    paero2.set_ydata([pcop[1], aeroEnd[1]])
    wing1, wing2, pcop, aeroEnd = flapkin(yu0[:,i], [0.02,0], params1)
    p3.set_xdata([wing1[0], wing2[0]])
    p3.set_ydata([wing1[1], wing2[1]])
    paero3.set_xdata([pcop[0], aeroEnd[0]])
    paero3.set_ydata([pcop[1], aeroEnd[1]])
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
