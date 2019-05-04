import autograd.numpy as np
from autograd import jacobian
import sys
from scipy.integrate import solve_ivp
import scipy.sparse as sparse
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits import mplot3d
import osqp

sys.path.append('..')
from controlutils.py import lqr, solver
from controlutils.py.models.pendulums import PendulumFlyWheel, LIP

np.set_printoptions(precision=4, suppress=True, linewidth=200)

lip = LIP()
ip = PendulumFlyWheel()

# LQR solutions here
yup = np.zeros(4)
uup = np.array([0.0])
A, B, c = lip.autoLin(yup, uup)
Q1 = np.eye(4)
R1 = np.eye(1)
# NOTE: the one step cost g(y,u) = y.Q.y + u.R.u
# LQR
K1, S1 = lqr.lqr(A, B, Q=Q1, R=R1)
# print(K1, P1)

yupip = np.array([0., 0.1, 0., 0., 0., 0.])
uupip = np.array([ip.m * ip.g, 0.0])
Aip, Bip, cip = ip.autoLin(yupip, uupip)
Qip = np.eye(6)
Rip = np.eye(2)
Kip, Sip = lqr.lqr(Aip, Bip, Q=Qip, R=Rip)
# print(Sip)

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
# probIP = osqp.OSQP()
# # QP decision var is u for IP
# nx = 2
# Pqp = sparse.csc_matrix([[R1, 0], [0, 0]])
# qqp = np.full(nx, 1)
# # no constraints, as before
# probIP.setup(Pqp, qqp, Aqp, lqp, uqp, verbose=False)

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
        qqp = B.T @ S1 @ y
        probIP.update(q=qqp)
        res = probIP.solve()
        return ip.dynamics(y, res.x)
    else:
        raise ValueError("Should be LIP or IP")


# Simulation
tf = 2
dt = 0.01
t_eval = np.arange(0, tf, dt)
y0 = np.array([2.,0.,0.,0.])
lipsol = solve_ivp(lambda t, y: lip.dynamics(y, K1 @ (yup - y)), [0, tf], y0, dense_output=True, t_eval=t_eval)
# LIP but with QP controller
lipqpsol = solve_ivp(valFuncQP, [0, tf], y0, dense_output=True, t_eval=t_eval)

# IP with LQR controller
def ipClosedLoop(t, y):
    # l = np.sqrt(y[0]**2 + y[1]**2)
    # dl = 2 * (y[0] * y[3] + y[1] * y[4])
    fk = 10000 * (lip.z0 - y[1]) - 100 * y[4]  # attract to z0
    tauh = 100 * y[0] + 1 * y[3]
    return ip.dynamics(y, np.array([fk, tauh]))
y0ip = np.array([0.1, 0.11, y0[1], y0[2], 0.0, y0[3]])
iplqrsol = solve_ivp(ipClosedLoop, [0, tf], y0ip, dense_output=True, t_eval=t_eval)
# IP with QP controller
# ipsol = solve_ivp(valFuncQP, [0, tf], y0ip, dense_output=True, t_eval=t_eval)

# make a list for display
dispsols = [
    {'name': 'liplqr', 'col': 'b', 'sol': lipsol, 'model': lip},
    {'name': 'lipqp', 'col': 'r', 'sol': lipqpsol, 'model': lip},
    {'name': 'iplqr', 'col': 'g', 'sol': iplqrsol, 'model': ip}
]

# ---

# visualize value function
def lqrValueFunc(x1, x2, P):
    # TODO: will need slices
    val = P[0,0] * (x1 - yup[0])**2 + P[1,1] * x2**2 + (P[0,1] + P[1,0]) * (x1 - yup[0])*x2
    return val

xx, yy = np.meshgrid(np.linspace(0, 2*np.pi, 30), np.linspace(-10, 10, 30))

# Display -------------------

fig, ax = plt.subplots(3)

for dispsol in dispsols:
    ax[0].plot(dispsol['sol'].t, dispsol['sol'].y[0, :], label=dispsol['name'])
# ax[0].plot(iplqrsol.t, iplqrsol.y[1,:])
ax[0].legend()

# ax[1].contourf(xx, yy, lqrValueFunc(xx, yy, pendulum['P']), cmap='gray_r')
# ax[1].plot(yup[0], yup[1], 'r*')

# Plot value along trajectory
lipval = np.zeros_like(lipsol.t)
for ti in range(len(lipval)):
    yi = lipsol.y[:, ti]
    lipval[ti] = yi @ S1 @ yi
ax[1].plot(lipsol.t, lipval)
ax[1].set_ylabel('CTG')

# animation
patches = []
for dispsol in dispsols:
    linei, = ax[2].plot([], [], lw=2, label=dispsol['name'])
    patches.append(linei)
ax[2].set_aspect(1)
ax[2].set_xlim((-0.3,0.3))
ax[2].set_ylim((-0.3,0.3))
ax[2].grid(True)
ax[2].legend(bbox_to_anchor=(2,1))
plt.tight_layout()

def _init():
    for pp in patches:
        pp.set_data([], [])
    return patches

def _animate(i):
    # get the vertices of the pendulum
    for mi in range(len(dispsols)):
        dispsol = dispsols[mi]
        if i < len(dispsol['sol'].t):
            p1 = dispsol['model'].kinematics(dispsol['sol'].y[:, i])
            patches[mi].set_data([0, p1[0]], [0, p1[1]])

    return patches

anim = animation.FuncAnimation(fig, _animate, init_func=_init, frames=len(lipsol.t), interval=1000*dt, blit=True)

plt.show()
