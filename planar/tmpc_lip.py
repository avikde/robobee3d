import autograd.numpy as np
from autograd import jacobian
import sys
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits import mplot3d

sys.path.append('..')
from controlutils.py import lqr, solver
from controlutils.py.models.pendulums import PendulumFlyWheel, LIP

np.set_printoptions(precision=4, suppress=True, linewidth=200)

lip = LIP()

yup = np.array([0.0, 0.0])
uup = np.array([0.0])
A, B, c = lip.autoLin(yup, uup)
# LQR
K1, P1 = lqr.lqr(A, B, Q=np.eye(2), R=0.01*np.eye(1))

print(K1, P1)

# Simulation
tf = 0.2
dt = 0.01
t_eval = np.arange(0, tf, dt)
y0 = np.array([2,0.0])
lipsol = solve_ivp(lambda t, y: lip.dynamics(y, K1 @ (yup - y)), [0, tf], y0, dense_output=True, t_eval=t_eval)

# make a list for display
dispsols = [
    {'name': 'lip', 'col': 'b', 'sol': lipsol, 'model': lip}
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
    ax[0].plot(dispsol['sol'].t, dispsol['sol'].y[0, :], '.-', label=dispsol['name'])
ax[0].legend()

# ax[1].contourf(xx, yy, lqrValueFunc(xx, yy, pendulum['P']), cmap='gray_r')
# ax[1].plot(yup[0], yup[1], 'r*')

# Plot value along trajectory
lipval = np.zeros_like(lipsol.t)
for ti in range(len(lipval)):
    yi = lipsol.y[:, ti]
    lipval[ti] = yi @ P1 @ yi
ax[1].plot(lipsol.t, lipval, '.-')
ax[1].set_ylabel('value')

# animation
patches = []
for dispsol in dispsols:
    linei, = ax[2].plot([], [], '.-', lw=2, label=dispsol['name'])
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
            p1 = dispsol['model'].kinematics(dispsol['sol'].y[0:1, i])
            patches[mi].set_data([0, p1[0]], [0, p1[1]])

    return patches

anim = animation.FuncAnimation(fig, _animate, init_func=_init, frames=len(lipsol.t), interval=1000*dt, blit=True)

plt.show()
