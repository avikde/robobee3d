import autograd.numpy as np
from autograd import jacobian
import sys
from scipy.integrate import solve_ivp
import scipy.sparse as sparse
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits import mplot3d
from matplotlib.collections import PatchCollection
import osqp
# import control

sys.path.append('..')
from controlutils.py import lqr, solver
import controlutils.py.models.pendulums as pendulums
import controlutils.py.models.aerial as aerial
import controlutils.py.misc as misc

np.set_printoptions(precision=4, suppress=True, linewidth=200)

# Create all the models
lip = pendulums.LIP()
ip = pendulums.PendulumFlyWheel()

# LQR sols here --------------------------------
# LIP
yup = np.zeros(4)
uup = np.array([0.0])
A, B, c = lip.autoLin(yup, uup)
Q1 = np.eye(4)
R1 = np.eye(1)
# NOTE: the one step cost g(y,u) = y.Q.y + u.R.u
# LQR
K1, S1 = lqr.lqr(A, B, Q=Q1, R=R1)
# print(K1, P1)

# IP
yupip = np.array([0., lip.z0, 0., 0., 0., 0.])
uupip = np.array([ip.m * pendulums.g, 0.0])
# check that is an equilibrium
assert np.allclose(ip.dynamics(yupip, uupip), np.zeros(6))
Aip, Bip, cip = ip.autoLin(yupip, uupip)
Qip = np.eye(6)
Rip = np.eye(2)
Kip, Sip = lqr.lqr(Aip, Bip, Q=Qip, R=Rip)
# print(Kip)


# ---
prob = osqp.OSQP()

# # the QP decision var, x = [y, u]
# nx = 5
# # col-major 1..nx^2 elements
# Pqp = np.reshape(np.arange(1,nx*nx+1), (nx,nx), order='F')
# Pqp = sparse.csc_matrix(Pqp)
# # now Pqp.data is just 1..nx^2
# qqp = np.zeros(nx)

# the QP decision var is u
nx = 1
# col-major 1..nx^2 elementss
Pqp = sparse.csc_matrix(R1)
# now Pqp.data is just 1..nx^2
qqp = np.full(nx, 1)

Aqp = sparse.csc_matrix(np.zeros((0,nx)))
lqp = np.zeros(0)
uqp = np.zeros(0)
prob.setup(Pqp, qqp, Aqp, lqp, uqp, verbose=False)

# ---
probIP = osqp.OSQP()
# QP decision var is u for IP
nx = 2
# add a ganch regularization of forces = u.Rip.u
Pqp = sparse.csc_matrix([[R1[0,0] + Rip[0,0], 0], [0, Rip[1,1]]])
qqp = np.full(nx, 1)
# no constraints, as before
probIP.setup(Pqp, qqp, Aqp, lqp, uqp, verbose=False)
# Proj to template
PiT = np.array([[1, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 1]])

# ---
def valFuncQP(t, y):
    global prob, probIP
    # update
    if len(y) == 4:
        qqp = B.T @ S1 @ y
        # prob.update(Px=np.reshape(Pqp, nx*nx, order='F'))
        prob.update(q=qqp)
        res = prob.solve()
        return lip.dynamics(y, res.x)
    elif len(y) == 6:
        # Linearize at the current state
        # TODO: store the jacobians and evaluate
        # Aip, Bip, cip = ip.autoLin(y, uupip)
        # Use value function from the template (with the appropriate Jacobians)--compare to above
        Vanch = np.array([0, 0 * (y[1] - lip.z0), 0, 0, -100 * y[4], 0])
        qqp = Bip.T @ (Vanch)# + PiT.T @ S1 @ PiT @ y)
        probIP.update(q=qqp)
        res = probIP.solve()
        return ip.dynamics(y, res.x)
    else:
        raise ValueError("Should be LIP or IP")


# Simulation setup
tf = 0.1
dt = 0.01
t_eval = np.arange(0, tf, dt)
y0 = np.array([2., 0., 0., 0.])  # for 4dim LIP
y0ip = np.array([0.1, 0.11, y0[1], y0[2], 0.0, y0[3]])  # for 6dim IP
lipsol = solve_ivp(lambda t, y: lip.dynamics(y, K1 @ (yup - y)), [0, tf], y0, dense_output=True, t_eval=t_eval)
# LIP but with QP controller
lipqpsol = solve_ivp(valFuncQP, [0, tf], y0, dense_output=True, t_eval=t_eval)

# IP with feedback controller
def ipClosedLoop(t, y):
    # l = np.sqrt(y[0]**2 + y[1]**2)
    # dl = 2 * (y[0] * y[3] + y[1] * y[4])
    # fk = 10000 * (lip.z0 - y[1]) - 100 * y[4]  # attract to z0
    # tauh = 100 * y[0] + 1 * y[3]
    # return ip.dynamics(y, np.array([fk, tauh]))

    return ip.dynamics(y, Kip @ (yupip - y))
iplqrsol = solve_ivp(ipClosedLoop, [0, tf], y0ip, dense_output=True, t_eval=t_eval)
# IP with QP controller
ipqpsol = solve_ivp(valFuncQP, [0, tf], y0ip, dense_output=True, t_eval=t_eval)

# Quadrotor 2d ---------------------------------------

# planar quadrotor
q2d = aerial.Quadrotor2D()

# LQR --
yhover = np.zeros(6)
uhover = np.full(2, q2d.m * aerial.g / 2.0)
# check that is an equilibrium
assert np.allclose(q2d.dynamics(yhover, uhover), np.zeros(6))
# q2d.fakeDamping = True
Aq2d, Bq2d, cq2d = q2d.autoLin(yhover, uhover)
Qq2d = np.diag([10, 10, 1, 1, 1, 0.1])
Rq2d = 0.001 * np.eye(2)
# print(Aq2d, Bq2d)  # , Kq2d)
Kq2d, Sq2d = lqr.lqr(Aq2d, Bq2d, Qq2d, Rq2d)
# print(Kq2d)

# QP for control with LQR VF --
probQ1 = osqp.OSQP()
# the QP decision var is u
nx = 2
# col-major 1..nx^2 elementss
Pqp = sparse.csc_matrix(Rq2d)
# now Pqp.data is just 1..nx^2
qqp = np.full(nx, 1)

Aqp = sparse.csc_matrix(np.zeros((0, nx)))
lqp = np.zeros(0)
uqp = np.zeros(0)
probQ1.setup(Pqp, qqp, Aqp, lqp, uqp, verbose=False)
uprev = uhover
def valFuncQuadQP(t, y):
    global probQ1, uprev
    # update
    if len(y) == 6:
        # Use value function from the template (with the appropriate Jacobians)--compare to above

        # Linearize at the current state
        # Bnow = jacobian(lambda u: q2d.dynamics(y, u))
        # qqp = (Bnow(uprev)).T @ Sq2d @ y

        qqp = (Bq2d).T @ Sq2d @ y
        # prob.update(Px=np.reshape(Pqp, nx*nx, order='F'))
        probQ1.update(q=qqp)
        res = probQ1.solve()
        uprev = res.x
        return q2d.dynamics(y, uprev)
    else:
        raise ValueError("Should be LIP or IP")

# Simulations --

tf = 3
dt = 0.05
t_eval = np.arange(0, tf, dt)
y0 = np.array([2, -1, 0, 0, 0, 0])
qlqrsol = solve_ivp(lambda t, y: q2d.dynamics(y, Kq2d @ (yhover - y)), [0, tf], y0, dense_output=True, t_eval=t_eval)
# ---
qsols = [
    {'name': 'lqr', 'col': 'b', 'sol': qlqrsol, 'model': q2d},
    {'name': 'qp', 'col': 'r--', 'sol': solve_ivp(valFuncQuadQP, [0, tf], np.array([2, -1, 0, 0, 0, 0]), dense_output=True, t_eval=t_eval), 'model': q2d},
]


# ------------ Display -----------------------

# Quadrotor display ---
fig, ax = plt.subplots(4)

for qsol in qsols:
    ax[0].plot(qsol['sol'].y[0, :], qsol['sol'].y[1, :], qsol['col'])
ax[0].set_aspect(1)
ax[0].set_ylabel('xz')

for qsol in qsols:
    ax[1].plot(qsol['sol'].t, qsol['sol'].y[2, :], qsol['col'])
ax[1].set_ylabel('phi')
plt.tight_layout()

# Plot CTG along trajectory
for qsol in qsols:
    q2dval = np.zeros_like(qsol['sol'].t)
    for ti in range(len(q2dval)):
        yi = qsol['sol'].y[:, ti]
        q2dval[ti] = yi @ Sq2d @ yi
    ax[2].plot(qsol['sol'].t, q2dval, qsol['col'])
ax[2].set_ylabel('CTG')

# Animation --
patches = [misc.rectangle(y0[0:2], y0[2], 2*q2d.r, 0.1*q2d.r, color=qsol['col'][0]) for qsol in qsols]
ax[3].set_aspect(1)
ax[3].set_xlim((-2,2))
ax[3].set_ylim((-1,1))
ax[3].grid(True)
ax[3].add_patch(patches[0])
ax[3].add_patch(patches[1])
plt.tight_layout()

def _init():
    return patches[0], patches[1], 

def _animate(i):
    for mi in range(len(qsols)):
        qsol = qsols[mi]['sol']
        if i < len(qsol.t):
            rawxy = misc.rectangle(qsol.y[0:2, i], qsol.y[2, i], 2*q2d.r, 0.5*q2d.r, rawxy=True)
            patches[mi].set_xy(rawxy)

    return patches[0], patches[1],

anim = animation.FuncAnimation(fig, _animate, init_func=_init, frames=len(qsols[0]['sol'].t), interval=1000*dt, blit=True)
# --


# # --- IP display ---

# dispsols = [
#     {'name': 'liplqr', 'col': 'b', 'sol': lipsol, 'model': lip},
#     {'name': 'lipqp', 'col': 'r', 'sol': lipqpsol, 'model': lip},
#     {'name': 'iplqr', 'col': 'g', 'sol': iplqrsol, 'model': ip},
#     # {'name': 'ipqp', 'col': 'g', 'sol': ipqpsol, 'model': ip}
# ]

# # visualize value function
# def lqrValueFunc(x1, x2, P):
#     # TODO: will need slices
#     val = P[0,0] * (x1 - yup[0])**2 + P[1,1] * x2**2 + (P[0,1] + P[1,0]) * (x1 - yup[0])*x2
#     return val

xx, yy = np.meshgrid(np.linspace(0, 2*np.pi, 30), np.linspace(-10, 10, 30))

# fig, ax = plt.subplots(3)

# # for dispsol in dispsols:
# #     ax[0].plot(dispsol['sol'].t, dispsol['sol'].y[0, :], label=dispsol['name'])
# ax[0].plot(ipqpsol.t, ipqpsol.y[0:2,:].T)
# ax[0].legend()

# # ax[1].contourf(xx, yy, lqrValueFunc(xx, yy, pendulum['P']), cmap='gray_r')
# # ax[1].plot(yup[0], yup[1], 'r*')

# # Plot value along trajectory
# lipval = np.zeros_like(lipsol.t)
# for ti in range(len(lipval)):
#     yi = lipsol.y[:, ti]
#     lipval[ti] = yi @ S1 @ yi
# ax[1].plot(lipsol.t, lipval)
# ax[1].set_ylabel('CTG')

# # animation
# patches = []
# for dispsol in dispsols:
#     linei, = ax[2].plot([], [], lw=2, label=dispsol['name'])
#     patches.append(linei)
# ax[2].set_aspect(1)
# ax[2].set_xlim((-0.3,0.3))
# ax[2].set_ylim((-0.3,0.3))
# ax[2].grid(True)
# ax[2].legend(bbox_to_anchor=(2,1))
# plt.tight_layout()

# def _init():
#     for pp in patches:
#         pp.set_data([], [])
#     return patches

# def _animate(i):
#     # get the vertices of the pendulum
#     for mi in range(len(dispsols)):
#         dispsol = dispsols[mi]
#         if i < len(dispsol['sol'].t):
#             p1 = dispsol['model'].kinematics(dispsol['sol'].y[:, i])
#             patches[mi].set_data([0, p1[0]], [0, p1[1]])

#     return patches

# anim = animation.FuncAnimation(fig, _animate, init_func=_init, frames=len(lipsol.t), interval=1000*dt, blit=True)

plt.show()
