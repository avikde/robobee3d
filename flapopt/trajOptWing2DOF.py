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

for ti in range(1, len(tvec)):
    # Nonlinear
    yi[:,ti] = yi[:,ti-1] + dt * closedLoop(tvec[ti], yi[:,ti-1])

# compare to continuous
sol = solve_ivp(closedLoop, [0,tf], yi[:,0], dense_output=True, t_eval=tvec)

print('Avg cost =', Jcostsol(sol.t, sol.y, params))

# -------- Open loop trajectory from the simulation above ------------------
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

def plotTrajs(*args):
    """Helper function to plot a bunch of trajectories superimposed"""
    umin, umax, xmin, xmax = m.limits
    _, ax = plt.subplots(3)
    for arg in args:
        ax[0].plot(olTrajt, arg[:,0], '.-')
    ax[0].axhline(y=xmin[0], color='k', alpha=0.3)
    ax[0].axhline(y=xmax[0], color='k', alpha=0.3)
    ax[0].set_ylabel('stroke (m)')
    for arg in args:
        ax[1].plot(olTrajt, arg[:,1], '.-')
    ax[1].axhline(y=xmin[1], color='k', alpha=0.3)
    ax[1].axhline(y=xmax[1], color='k', alpha=0.3)
    ax[1].set_ylabel('hinge angle (rad)')
    for arg in args:
        ax[2].plot(olTrajt, arg[:,4], '.-')
    ax[2].axhline(y=umin[0], color='k', alpha=0.3)
    ax[2].axhline(y=umax[0], color='k', alpha=0.3)
    ax[2].set_ylabel('stroke force (N)')
    plt.show()
    sys.exit(0)
# --------------------------------------

# Test the linearized dynamics 
yi2 = olTraj.copy()
yilin = olTraj.copy()
for ti in range(1, len(olTrajt)):
    ui = yilin[ti-1, 4:]
    # Nonlinear
    yi2[ti, :4] = yi2[ti-1, :4] + olTrajdt * m.dydt(yi2[ti-1, :4], ui, params)
    # Linearized
    # ui = np.array([controller(tvec[ti], yilin[:,ti-1])])
    A, B, c = m.getLinearDynamics(yilin[ti-1, :4], ui)
    yilin[ti, :4] = A @ yilin[ti-1, :4] + B @ ui + c
# plotTrajs(olTraj, yi2, yilin)

# Wing traj opt using QP -------------------------------------------------
def dirTranForm(xtraj, N, nx, nu):
    # convert from the (N-1,nx+nu) array to the dirtran form with x0,..,x[N-1],u0,...,u[N-2] i.e. shorter by 1
    return np.hstack((
            np.reshape(xtraj[:,:nx], (N+1)*nx, 'C'),
            np.reshape(xtraj[:-1,nx:], (N)*nu, 'C')
        ))
def xuMatForm(dirtranx, N, nx):
    Xmat = np.reshape(dirtranx[:(N)*nx], (N,nx), 'C')
    Umat = dirtranx[(N)*nx:][:,np.newaxis]
    # Need to add another u since dirtran form has one fewer u. Just copy the last one
    Umat = np.vstack((Umat, Umat[-1,:]))
    return np.hstack((Xmat, Umat))
# Create "MPC" object which will be used for SQP
class QOFAvgLift:
    def __init__(self, N, wx, wu, kdampx, kdampu):
        self.N = N
        self.nx = len(wx)
        self.nu = len(wu)
        self.wx = wx
        self.wu = wu
        self.kdampx = kdampx
        self.kdampu = kdampu

    def getPq(self, xtraj):
        nX = (self.N+1) * self.nx + self.N*self.nu
        # vector of weights for the whole dirtran x
        w = np.hstack((np.tile(self.wx, self.N+1), np.tile(self.wu, self.N)))
        kdamp = np.hstack((np.tile(self.kdampx, self.N+1), np.tile(self.kdampu, self.N)))

        self.P = sparse.diags(w + kdamp).tocsc()
        dirtranx = dirTranForm(xtraj, self.N, self.nx, self.nu)
        self.q = -np.multiply(w, dirtranx)
        return self.P, self.q

class WingQP:
    def __init__(self, model, N, wx, wu, kdampx, kdampu, **settings):
        self.ltvsys = ltvsystem.LTVSolver(model)
        # Dynamics and constraints
        self.ltvsys.initConstraints(model.nx, model.nu, N, polyBlocks=None)
        self.ltvsys.initObjective(QOFAvgLift(N, wx, wu, kdampx, kdampu))
        self.ltvsys.initSolver(**settings)

    def update(self, xtraj):
        N = self.ltvsys.N
        nx = self.ltvsys.nx
        # TODO: check which traj mode
        u0 = xtraj[:,4][:,np.newaxis]
        # NOTE: confirmed that updateTrajectory correctly updates the traj, and that updateDynamics is updating the A, B
        xtraj = self.ltvsys.updateTrajectory(xtraj[:,:4], u0, trajMode=ltvsystem.GIVEN_POINT_OR_TRAJ)
        self.ltvsys.updateObjective()
        self.dirtranx, res = self.ltvsys.solve(throwOnError=False)
        if res.info.status not in ['solved', 'solved inaccurate', 'maximum iterations reached']:
            self.ltvsys.debugResult(res)
            # dirtranx = dirTranForm(xtraj, N, nx, self.ltvsys.nu)
            # print(self.ltvsys.u - self.ltvsys.A @ dirtranx, self.ltvsys.A @ dirtranx - self.ltvsys.l)

            raise ValueError(res.info.status)
        # debug
        # print(self.ltvsys.u - self.ltvsys.A @ dirtranx, self.ltvsys.A @ dirtranx - self.ltvsys.l)
        # reshape into (N,nx+nu)
        return xuMatForm(self.dirtranx, N+1, nx)

# Use the class above to step the QP
Nknot = olTraj.shape[0]  # number of knot points in this case
nx = 4
nu = 1
wx = np.ones(nx) * 10
wu = np.ones(nu) * 0.1
kdampx = np.ones(4)
kdampu = np.ones(1)
# Must be 1 smaller to have the correct number of xi
wqp = WingQP(m, Nknot-1, wx, wu, kdampx, kdampu, verbose=False, eps_rel=1e-2, eps_abs=1e-2)
# Test warm start
# wqp.ltvsys.prob.warm_start(x=dirTranForm(olTraj, Nknot, 4, 1))
traj2 = wqp.update(olTraj)
# print(olTraj - wqp.ltvsys.xtraj) # <these are identical: OK; traj update worked

# # Debug the solution
# print(olTraj.shape, wqp.dirtranx.shape)
# print(wqp.ltvsys.A.shape, olTraj.shape)
# olTrajDirTran = dirTranForm(olTraj, wqp.ltvsys.N, wqp.ltvsys.nx,  wqp.ltvsys.nu)
# traj2DirTran = dirTranForm(traj2, wqp.ltvsys.N, wqp.ltvsys.nx,  wqp.ltvsys.nu)
# fig, ax = plt.subplots(2)
# ax[0].plot(wqp.ltvsys.A @ olTrajDirTran - wqp.ltvsys.l, label='1')
# ax[0].plot(wqp.ltvsys.A @ traj2DirTran - wqp.ltvsys.l, label='3')
# ax[0].axhline(0, color='k', alpha=0.3)
# ax[0].legend()
# ax[1].plot(wqp.ltvsys.u - wqp.ltvsys.A @ olTrajDirTran, label='1')
# ax[1].plot(wqp.ltvsys.u - wqp.ltvsys.A @ traj2DirTran, label='3')
# ax[1].axhline(0, color='k', alpha=0.3)
# ax[1].legend()
# plt.show()
# sys.exit(0)

# print(olTraj.shape, traj2.shape, olTrajt.shape)
plotTrajs(olTraj, traj2)# debug the 1-step solution

# wx = np.array([1,1,1,1])
# wu = np.array([1])
# peps = 1e-2
# ltvqp = mpc.LTVMPC(m, Nknot, wx, wu, verbose=True, polish=False, scaling=0, eps_rel=peps, eps_abs=peps, max_iter=100, kdamping=0)
# # x0 must be a (N,nx) trajectory
# xr = np.zeros(m.nx)  # FIXME: does not make sense
# ctrl = ltvqp.update(x0=olTraj, xr=xr, trajMode=ltvsystem.GIVEN_POINT_OR_TRAJ)

sys.exit(0)

# OLD gradient descent ------------
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
