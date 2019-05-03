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

print(A)

# # Simulation
# tf = 5.0
# dt = 0.01
# t_eval = np.arange(0, tf, dt)
# y0 = np.array([2,0.0])
# pendulum['sol'] = solve_ivp(lambda t, y: pendulum['model'].dynamics(y, K1 @ (yup - y)), [0, tf], y0, dense_output=True, t_eval=t_eval)

# # ---

# # visualize value function
# def lqrValueFunc(x1, x2, P):
#     # quadratic
#     val = P[0,0] * (x1 - yup[0])**2 + P[1,1] * x2**2 + (P[0,1] + P[1,0]) * (x1 - yup[0])*x2
#     return val

# xx, yy = np.meshgrid(np.linspace(0, 2*np.pi, 30), np.linspace(-10, 10, 30))

# # Display -------------------

# fig, ax = plt.subplots(3)

# ax[0].plot(pendulum['sol'].t, pendulum['sol'].y[0, :], label='sp')
# ax[0].plot(pendulum2['sol'].t, pendulum2['sol'].y[:2, :].T, label='dp')
# ax[0].plot(pendulum2['sol2'].t, pendulum2['sol2'].y[:2, :].T, label='dpmpc')
# ax[0].plot(acrobot['sol2'].t, acrobot['sol2'].y[:2, :].T, label='acro')
# ax[0].legend()

# ax[1].contourf(xx, yy, lqrValueFunc(xx, yy, pendulum['P']), cmap='gray_r')
# ax[1].plot(yup[0], yup[1], 'r*')

# # animation
# line1, = ax[2].plot([], [], '.-', lw=2, label='sp')
# line2, = ax[2].plot([], [], '.-', lw=2, label='dp')
# line3, = ax[2].plot([], [], '.-', lw=2, label='dpmpc')
# line4, = ax[2].plot([], [], '.-', lw=2, label='acro')
# patches = [line1, line2, line3, line4]
# ax[2].set_aspect(1)
# ax[2].set_xlim((-2,2))
# ax[2].set_ylim((-2,2))
# ax[2].grid(True)
# ax[2].legend(bbox_to_anchor=(2,1))
# plt.tight_layout()

# def _init():
#     line1.set_data([], [])
#     line2.set_data([], [])
#     line3.set_data([], [])
#     line4.set_data([], [])
#     return patches

# def _animate(i):
#     # get the vertices of the pendulum
#     if i < len(pendulum['sol'].t):
#         p1 = pendulum['model'].kinematics(pendulum['sol'].y[0:1, i])
#         line1.set_data([0, p1[0]], [0, p1[1]])

#     if i < len(pendulum2['sol'].t):
#         p1, p2 = pendulum2['model'].kinematics(pendulum2['sol'].y[0:2, i])
#         line2.set_data([0, p1[0], p2[0]], [0, p1[1], p2[1]])

#     if i < len(pendulum2['sol2'].t):
#         p1, p2 = pendulum2['model'].kinematics(pendulum2['sol2'].y[0:2, i])
#         line3.set_data([0, p1[0], p2[0]], [0, p1[1], p2[1]])

#     if i < len(acrobot['sol2'].t):
#         p1, p2 = acrobot['model'].kinematics(acrobot['sol2'].y[0:2, i])
#         line4.set_data([0, p1[0], p2[0]], [0, p1[1], p2[1]])
#     return patches

# anim = animation.FuncAnimation(fig, _animate, init_func=_init, frames=len(pendulum2['sol'].t), interval=1000*dt, blit=True)

# plt.show()
