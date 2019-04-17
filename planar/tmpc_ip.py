'''
Test template MPC ideas on:
- actuated IP as template
- double pendulum, acrobot as anchors
'''
import autograd.numpy as np
from autograd import jacobian
import sys
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

sys.path.append('..')
from controlutils.py import lqr

np.set_printoptions(precision=4, suppress=True, linewidth=200)

'''
Simulate an IP
'''

def pendulum(y, u):
    '''Inverted pendulum actuated at the base'''
    g = 9.81
    l = 1.0
    kd = 0.1
    # y = (theta, dtheta)
    return np.array([y[1], -g * np.sin(y[0]) / l - kd * y[1] + u])

# auto linearization
yup = np.array([np.pi, 0.0])
uup = 0.0
Afun = jacobian(lambda y: pendulum(y, uup))
Bfun = jacobian(lambda u: pendulum(yup, u))
B = Bfun(uup)[:,np.newaxis]
Aup = Afun(yup)
# LQR
Q = np.eye(2)
R = np.eye(1)
K, P = lqr.lqr(Aup, B, Q, R)

# Simulation
y0 = np.array([2,0.0]) 
sol = solve_ivp(lambda t, y: pendulum(y, K @ (yup - y)), [0, 10], y0, dense_output=True)
t, y = sol.t, sol.y

# visualize value function
def f(x1,x2):
    # quadratic
    val = P[0,0] * x1**2 + P[1,1] * x2**2 + (P[0,1] + P[1,0]) * x1*x2
    return val

xx, yy = np.meshgrid(np.linspace(-np.pi, np.pi, 30), np.linspace(-10, 10, 30))
zz = f(xx, yy)

# Display -------------------

fig, ax = plt.subplots(2)

ax[0].plot(t, y[0,:])
ax[1].contourf(xx, yy, zz, cmap='gray')

# # plot
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot_surface(xx,yy,zz, rstride=1, cstride=1, cmap='viridis', edgecolor='none')

plt.show()
