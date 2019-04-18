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
g = 9.81

def pendulum(y, u):
    '''Inverted pendulum actuated at the base'''
    l = 1.0
    kd = 0.1
    # y = (theta, dtheta)
    return np.array([y[1], -g * np.sin(y[0]) / l - kd * y[1] + u])

def doublePendulum(y, u):
    l1 = 1.0
    l2 = 1.0
    m1 = 1.0
    m2 = 1.0
    # unpack
    if len(u) > 1:
        # double pendulum
        tau1 = u[0]
        tau2 = u[1]
    else:
        # acrobot
        tau1 = 0
        tau2 = u[0] 
    dq1 = y[2]
    dq2 = y[3]
    s1 = np.sin(y[0])
    c1 = np.cos(y[0])
    s2 = np.sin(y[1])
    c2 = np.cos(y[1])
    s12 = np.sin(y[0] + y[1])
    dq12 = dq1**2
    dq22 = dq2**2
    l22 = l2**2
    l12 = l1**2
    m22 = m2**2

    # from mathematica
    y2dot = np.array([(g*l1*l2*m1*s1 + g*l1*l2*m2*s1 - c2*g*l1*l2*m2*s12 + 
      dq12*l1*l2*(c2*l1 + l2)*m2*s2 + 2*dq1*dq2*l1*l22*m2*s2 + 
      dq22*l1*l22*m2*s2 + l2*tau1 - c2*l1*tau2 - l2*tau2)/
    (l1**2*l2*(m1 + m2 - c2**2*m2)), 
    -((g*l1*l22*m1*m2*s1 + g*l1*l22*m22*s1 - c2**2*g*l1*l22*m22*s1 - 
        c1*g*l12*l2*m1*m2*s2 + 
        dq12*l1*l2*m2*(2*c2*l1*l2*m2 + l22*m2 + l12*(m1 + m2))*s2 - 
        c1*g*l12*l2*m22*s2 - c1*c2*g*l1*l22*m22*s2 + 
        2*dq1*dq2*l1*(c2*l1 + l2)*l22*m22*s2 + 
        dq22*l1*(c2*l1 + l2)*l22*m22*s2 + c2*l1*l2*m2*tau1 + l22*m2*tau1 - 
        l12*m1*tau2 - l12*m2*tau2 - 2*c2*l1*l2*m2*tau2 - l22*m2*tau2)/
      (l1**2*l2**2*m2*(m1 + m2 - c2**2*m2)))
    ])

    return np.hstack((y[2:], y2dot))

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
# double pendulum
y02 = np.array([1,-1, 0, 0])
sol2 = solve_ivp(lambda t, y: doublePendulum(y, np.zeros(2)), [0, 10], y02, dense_output=True)

# visualize value function
def lqrValueFunc(x1,x2):
    # quadratic
    val = P[0,0] * (x1 - yup[0])**2 + P[1,1] * x2**2 + (P[0,1] + P[1,0]) * (x1 - yup[0])*x2
    return val

xx, yy = np.meshgrid(np.linspace(0, 2*np.pi, 30), np.linspace(-10, 10, 30))

# Display -------------------

fig, ax = plt.subplots(2)

ax[0].plot(sol.t, sol.y[0,:], label='sp')
ax[0].plot(sol2.t, sol2.y[0,:], label='dp0')
ax[0].plot(sol2.t, sol2.y[1,:], label='dp1')
ax[0].legend()

ax[1].contourf(xx, yy, lqrValueFunc(xx, yy), cmap='gray_r')
ax[1].plot(yup[0], yup[1], 'r*')

# # plot
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot_surface(xx,yy,zz, rstride=1, cstride=1, cmap='viridis', edgecolor='none')

plt.show()
