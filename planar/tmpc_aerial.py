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
q2d = {'m': aerial.Quadrotor2D()}
ptsd = {'m': PlanarThrustStrokeDev2()}

# LQR --
# Hover conditions
q2d['y0'] = np.zeros(6)
ptsd['y0'] = np.zeros(6)
q2d['u0'] = np.full(2, q2d['m'].m * aerial.g / 2.0)
ptsd['u0'] = np.array([ptsd['m'].m * aerial.g, 0])

# Costs
q2d['Q'] = np.diag([10, 10, 1, 1, 1, 0.1])
q2d['R'] = 0.001 * np.eye(2)
ptsd['Q'] = np.diag([10, 10, 1, 1, 1, 0.1])
ptsd['R'] = np.diag([0.001, 0.1])

# Some computation for all the systems
for S in [q2d, ptsd]:
    # check that is an equilibrium
    assert np.allclose(S['m'].dynamics(S['y0'], S['u0']), np.zeros(6))
    # q2d.fakeDamping = True
    S['A'], S['B'], S['c'] = S['m'].autoLin(S['y0'], S['u0'])
    # print(Aq2d, Bq2d)  # , Kq2d)
    S['K'], S['S'] = lqr.lqr(S['A'], S['B'], S['Q'], S['R'])
    print(S['K'])
    # LQR solution u* = - K y = - R^{-1} B^T S y
    assert np.allclose(np.linalg.inv(S['R']) @ S['B'].T @ S['S'], S['K'])
# sys.exit(0)
# QP for control with LQR VF --
probQ1 = osqp.OSQP()
# the QP decision var is u
nx = 2
# col-major 1..nx^2 elementss
Pqp = sparse.csc_matrix(q2d['R'])
# now Pqp.data is just 1..nx^2
qqp = np.full(nx, 1)

Aqp = sparse.csc_matrix(np.zeros((0, nx)))
lqp = np.zeros(0)
uqp = np.zeros(0)
probQ1.setup(Pqp, qqp, Aqp, lqp, uqp, verbose=False, eps_rel=1e-4, eps_abs=1e-4, max_iter=20)
uprev = q2d['u0']

def valFuncQuadQP(t, y, anch):
    global probQ1, uprev
    # update
    # Use value function from the template (with the appropriate Jacobians)--compare to above

    if anch == 'q2d':
        # Linearize at the current state
        # Bnow = jacobian(lambda u: q2d.dynamics(y, u))
        # qqp = (Bnow(uprev)).T @ Sq2d @ y

        # B from eq
        qqp = (q2d['B']).T @ q2d['S'] @ y
    elif anch == 'ptsd':
        # FIXME: autograd not working on this
        Bnow = jacobian(lambda u: ptsd['m'].dynamics(y, u))
        # qqp = (Bnow(uhoverPTSD)).T @ Sq2d @ y

        # u1 = ptsd.m * aerial.g
        # u2 = 0
        u1 = uprev[0]
        u2 = uprev[1]
        # print(uhoverPTSD, u1, u2)
        sphi = np.sin(y[2])
        cphi = np.cos(y[2])
        mb = ptsd['m'].m
        ib = ptsd['m'].ib
        Bsym = np.vstack((np.zeros((3, 2)),
                          np.array([[-sphi/mb, 0], [cphi/mb, 0], [u2/ib, u1/ib]])))
        # print(Bnow(uhoverPTSD))
        # print(Bnow(uhoverPTSD), Bsym)
        qqp = q2d['B'].T @ q2d['S'] @ y
    else:
        raise ValueError("specify anch")

    # prob.update(Px=np.reshape(Pqp, nx*nx, order='F'))
    probQ1.update(q=qqp)
    res = probQ1.solve()
    uprev = res.x
    if anch == 'q2d':
        return q2d['m'].dynamics(y, uprev)
    elif anch == 'ptsd':
        uprev[0] = 5
        uprev[1] = 0
        print(uprev)
        return ptsd['m'].dynamics(y, uprev)

# Simulations --


tf = 2
dt = 0.05
t_eval = np.arange(0, tf, dt)
y0 = np.array([2, -1, 0, 0, 0, 0])
qlqrsol = solve_ivp(lambda t, y: q2d['m'].dynamics(y, q2d['K'] @ (q2d['y0'] - y)), [0, tf], y0, dense_output=True, t_eval=t_eval)

qqpsol = solve_ivp(lambda t, y: valFuncQuadQP(t, y, 'q2d'), [0, tf], np.array([2, 1, 0, 0, 0, 0]), dense_output=True, t_eval=t_eval)

# NOTE: using S from PTSD
Ktest = np.linalg.inv(ptsd['R']) @ ptsd['B'].T @ q2d['S']  # ptsd['K']
# print(Ktest)
# print(q2d['S'] - ptsd['S'])
# print(np.linalg.inv(ptsd['R']) @ ptsd['B'].T)
# sys.exit(0)
ptsdlqrsol = solve_ivp(lambda t, y: ptsd['m'].dynamics(y, Ktest @ (ptsd['y0'] - y)), [0, tf], y0, dense_output=True, t_eval=t_eval)

uprev = ptsd['u0']
probQ1.update(Px=np.array([1, 10]))  # more input weight
qptsdsol = solve_ivp(lambda t, y: valFuncQuadQP(t, y, 'ptsd'), [0, tf], np.array([1, 0, 0, 0, 0, 0]), dense_output=True, t_eval=t_eval)

# ---
qsols = [
    {'name': 'lqr', 'col': 'b', 'sol': qlqrsol, 'model': q2d['m']},
    {'name': 'qp', 'col': 'r', 'sol': qqpsol, 'model': q2d['m']},
    {'name': 'ptsdlqr', 'col': 'g', 'sol': ptsdlqrsol, 'model': ptsd['m']},
    {'name': 'ptsd', 'col': 'c', 'sol': qptsdsol, 'model': ptsd['m']},
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
        q2dval[ti] = yi @ q2d['S'] @ yi
    ax[2].plot(qsol['sol'].t, q2dval, qsol['col'])
ax[2].set_ylabel('CTG')

# Animation --
r = q2d['m'].r
patches = [misc.rectangle(y0[0:2], y0[2], 2*r, 0.1*r, color=qsol['col'][0]) for qsol in qsols]
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
            rawxy = misc.rectangle(qsol.y[0:2, i], qsol.y[2, i], 2*r, 0.5*r, rawxy=True)
            patches[mi].set_xy(rawxy)

    return patches[0], patches[1], patches[2],


anim = animation.FuncAnimation(fig, _animate, init_func=_init, frames=len(qsols[0]['sol'].t), interval=1000*dt, blit=True)
# --

plt.show()
