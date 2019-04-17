'''
Test template MPC ideas on:
- actuated IP as template
- double pendulum, acrobot as anchors
'''
import numpy as np
import sys
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

sys.path.append('..')
from controlutils.py import lqr

np.set_printoptions(precision=4, suppress=True, linewidth=200)

'''
Simulate an IP
'''

def pendulum(t, y):
    '''Inverted pendulum actuated at the base'''
    g = 9.81
    l = 1.0
    # y = (theta, dtheta)
    return np.array([y[1], -g * y[0] / l])

sol = solve_ivp(pendulum, [0, 10], np.array([1.0,0.0]), dense_output=True)
t, y = sol.t, sol.y

# Display

fig, ax = plt.subplots(2)

ax[0].plot(t, y[0,:])
ax[1].plot(t, y[1,:])

plt.show()
