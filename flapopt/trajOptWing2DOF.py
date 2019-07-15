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

"""QP based ----------------------------------------
"""
# # Use the class above to step the QP
Nknot = olTraj.shape[0]  # number of knot points in this case
# nx = 4
# nu = 1
# wx = np.ones(nx) * 1e-6
# wu = np.ones(nu) * 1e3
# kdampx = 1e-5 * np.ones(4)
# kdampu = np.zeros(1)
# # Must be 1 smaller to have the correct number of xi
# wqp = wingopt.WingQP(m, Nknot-1, wx, wu, kdampx, kdampu, verbose=False, eps_rel=1e-4, eps_abs=1e-4, max_iter=10000)
# wqp.trajt = olTrajt
# # Test warm start
# # wqp.ltvsys.prob.warm_start(x=dirTranForm(olTraj, Nknot, 4, 1))
# traj2 = wqp.update(olTraj)
# traj3 = wqp.update(traj2)
# # TEST: reduce wu 
# # wqp.ltvsys.qof.wu = np.ones(nu) * 1e2
# traj4 = wqp.update(traj3)
# # wqp.debugConstraintViol(olTraj, wqp.dirtranx)

# # # Optim wrt params ----

# # JT = lambda T : wingopt.Jcost_dirtran(wqp.dirtranx, Nknot, [wingopt.params[0], T])
# # DJT = jacobian(JT)

# # Ttest = np.arange(0.5, 1.5, 0.1)
# # plt.plot(Ttest, [JT(Ti) for Ti in Ttest])

# # Display -------------

# # tvec, ctstrajs = wingopt.createCtsTraj(dt, olTrajt, [olTraj, traj2, traj3, traj4])
# # wqp.plotTrajs(olTraj, traj2, traj3, traj4)
# # wingopt.trajAnim(tvec, ctstrajs)
# # plt.show()
# # sys.exit(0)

"""
Penalty-based NL optim ----------------------------------------
"""

wpo = wingopt.WingPenaltyOptimizer(Nknot-1)
# Initial trajectory
traj0 = wingopt.dirTranForm(olTraj, Nknot-1, m.nx, m.nu)
params0 = wingopt.params
# with params as well
# traj0 = np.hstack((traj0, wingopt.params))

optavglift = {'dynamics':1e-3, 'periodic':0, 'input':1e4, 'state': 1e0}
optavgliftparams = {'dynamics':1e3, 'periodic':0, 'input':1e4, 'state': 1e0}

print("hi 0")
trajs = [traj0]
params = [params0]
# test
# params = [[0.2, 2.0]]
for ii in range(2):
    trajs.append(wpo.update(trajs[-1], params[-1], mode=wpo.WRT_TRAJ, opt=optavglift)[0])
    # traj2 = wpo.update(traj0, opt={'dynamics':1e-3, 'periodic':0, 'input':1e4, 'state': 1e0, 'odrag': 1})
    pnew, J, _ = wpo.update(trajs[-1], params[-1], mode=wpo.WRT_PARAMS, opt=optavglift)
    # pnew, J, _ = wpo.update(trajs[-1], params[-1], mode=wpo.WRT_TRAJ_PARAMS, opt=optavgliftparams)
    params.append(pnew[-len(params0):])

print(params)
wpo.plotTrajs(*trajs)

# # Try optimizing traj wrt given params separately
# cbars = np.linspace(0.002, 0.02, 5)
# Ts = np.linspace(0.7, 2, 5)
# cbarsM, TsM = np.meshgrid(cbars, Ts)
# JbestsM = np.zeros_like(cbarsM)
# for ii in range(JbestsM.shape[0]):
#     for jj in range(JbestsM.shape[1]):
#         print("Progress", JbestsM.shape, ii, jj)
#         # no opt, just change params
#         c, r = wingopt.Jcosttraj_penalty(trajs[-1], wpo.N, [cbarsM[ii,jj], TsM[ii,jj]], opt=optavgliftparams)
#         JbestsM[ii,jj] = c + r.T @ r
#         # # opt traj with these params
#         # _, _, JbestsM[ii,jj] = wpo.update(trajs[-1], [cbarsM[ii,jj], TsM[ii,jj]], mode=wpo.WRT_TRAJ, opt=optavgliftparams)

# print(JbestsM)
# # np.savez('res', cbars, Ts, JbestsM)
# # display this
# from mpl_toolkits.mplot3d import axes3d
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(cbarsM, TsM, JbestsM, cmap=plt.get_cmap('gist_earth'))
# ax.set_xlabel('cbar')
# ax.set_ylabel('T')

# tvec, ctstrajs = wingopt.createCtsTraj(dt, olTrajt, [traj0, traj])
# wingopt.trajAnim(tvec, ctstrajs)
plt.show()
sys.exit(0)

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
