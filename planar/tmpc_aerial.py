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

sys.path.append('..')
from controlutils.py import lqr, solver
import controlutils.py.models.aerial as aerial
import controlutils.py.misc as misc
from FlappingModels import PlanarThrustStrokeDev2

np.set_printoptions(precision=4, suppress=True, linewidth=200)

# Quadrotor 2d ---------------------------------------

# planar quadrotor
q2d = aerial.Quadrotor2D()
ptsd = PlanarThrustStrokeDev2()

# LQR --
yhover = np.zeros(6)
uhover = np.full(2, q2d.m * aerial.g / 2.0)
uhoverPTSD = np.array([ptsd.m * aerial.g, 0])
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
probQ1.setup(Pqp, qqp, Aqp, lqp, uqp, verbose=False, eps_rel=1e-4, eps_abs=1e-4, max_iter=20)
uprev = uhover


def valFuncQuadQP(t, y, anch):
    global probQ1, uprev
    # update
    # Use value function from the template (with the appropriate Jacobians)--compare to above

    if anch == 'q2d':
        # Linearize at the current state
        # Bnow = jacobian(lambda u: q2d.dynamics(y, u))
        # qqp = (Bnow(uprev)).T @ Sq2d @ y

        # B from eq
        qqp = (Bq2d).T @ Sq2d @ y
    elif anch == 'ptsd':
        # FIXME: autograd not working on this
        Bnow = jacobian(lambda u: ptsd.dynamics(y, u))
        # qqp = (Bnow(uhoverPTSD)).T @ Sq2d @ y

        # u1 = ptsd.m * aerial.g
        # u2 = 0
        u1 = uprev[0]
        u2 = uprev[1]
        # print(uhoverPTSD, u1, u2)
        sphi = np.sin(y[2])
        cphi = np.cos(y[2])
        Bsym = np.vstack((np.zeros((3, 2)),
                          np.array([[-sphi / ptsd.m, 0], [cphi/ptsd.m, 0], [u2/ptsd.ib, u1/ptsd.ib]])))
        # print(Bnow(uhoverPTSD))
        # print(Bnow(uhoverPTSD), Bsym)
        qqp = Bq2d.T @ Sq2d @ y
    else:
        raise ValueError("specify anch")

    # prob.update(Px=np.reshape(Pqp, nx*nx, order='F'))
    probQ1.update(q=qqp)
    res = probQ1.solve()
    uprev = res.x
    if anch == 'q2d':
        return q2d.dynamics(y, uprev)
    elif anch == 'ptsd':
        uprev[0] = 5
        uprev[1] = 0
        print(uprev)
        return ptsd.dynamics(y, uprev)

# Simulations --


tf = 3
dt = 0.05
t_eval = np.arange(0, tf, dt)
y0 = np.array([2, -1, 0, 0, 0, 0])
qlqrsol = solve_ivp(lambda t, y: q2d.dynamics(y, Kq2d @ (yhover - y)), [0, tf], y0, dense_output=True, t_eval=t_eval)

qqpsol = solve_ivp(lambda t, y: valFuncQuadQP(t, y, 'q2d'), [0, tf], np.array([2, 1, 0, 0, 0, 0]), dense_output=True, t_eval=t_eval)

uprev = uhoverPTSD
probQ1.update(Px=np.array([1, 10]))  # more input weight
qptsdsol = solve_ivp(lambda t, y: valFuncQuadQP(t, y, 'ptsd'), [0, tf], np.array([1, 0, 0, 0, 0, 0]), dense_output=True, t_eval=t_eval)
# ---
qsols = [
    {'name': 'lqr', 'col': 'b', 'sol': qlqrsol, 'model': q2d},
    {'name': 'qp', 'col': 'r', 'sol': qqpsol, 'model': q2d},
    {'name': 'ptsd', 'col': 'g', 'sol': qptsdsol, 'model': ptsd},
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
ax[3].set_xlim((-2, 2))
ax[3].set_ylim((-1, 1))
ax[3].grid(True)
for patch in patches:
    ax[3].add_patch(patch)
plt.tight_layout()


def _init():
    return patches[0], patches[1], patches[2],


def _animate(i):
    for mi in range(len(qsols)):
        qsol = qsols[mi]['sol']
        if i < len(qsol.t):
            rawxy = misc.rectangle(qsol.y[0:2, i], qsol.y[2, i], 2*q2d.r, 0.5*q2d.r, rawxy=True)
            patches[mi].set_xy(rawxy)

    return patches[0], patches[1], patches[2],


anim = animation.FuncAnimation(fig, _animate, init_func=_init, frames=len(qsols[0]['sol'].t), interval=1000*dt, blit=True)
# --

plt.show()
