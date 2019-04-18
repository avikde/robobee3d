import osqp
import numpy as np
import scipy as sp
import scipy.linalg
# import scipy.sparse as sparse
import matplotlib.pyplot as plt
import sys

# # Solves the problem 
# min 1/2 x^T P x + q^T x
# s.t. l <= A x <= u
# prob.setup(P, q, A, l, u, warm_start=True)

# Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
# i.e. all the states, and all the u's
# Kronecker product, a composite array made of blocks of the second array scaled by the first.
# >>> np.kron(np.eye(2), np.ones((2,2)))
# array([[ 1.,  1.,  0.,  0.],
#        [ 1.,  1.,  0.,  0.],
#        [ 0.,  0.,  1.,  1.],
#        [ 0.,  0.,  1.,  1.]])

# Looks like P contains quadratic error for both state and control effort
# q contains the goal error (x - xr)
# the dynamics are contained in the inequality
# Ad,Bd -> Ax, Bu
# Aeq = [Ax, Bu]
# Aineq has state and input ineq constraints [xmin, xmax], [umin, umax]

# ---- Original example https://osqp.org/docs/examples/mpc.html

CDmax = 3.4
CLmax = 1.8
CD0 = 0.4
# Units: mm, mg
khinge0 = 0.1
mb = 100 # mg
ib = 1./12. * mb * 12**2 # 12mm long body
kaero = 1
d = 5
dt = 0.002
tf = 0.5

# Discrete time model
def getLin(u, sigma, phi):
    # q = (x,z,phi), y = (sigma, q, dq)

    # Need to reverse hinge deflection
    if u < 0:
        khinge = -khinge0
    else:
        khinge = khinge0

    Ad = np.matrix([
        [0., 0., 0., 0, 0., 0., 0.],
        [0., 0., 0., 0, 1., 0., 0.],
        [0., 0., 0., 0., 0., 1., 0.],
        [0., 0., 0., 0., 0., 0., 1.],
        [0., 0.,      0.,     -(kaero*u*(2*CLmax*np.cos(phi)*np.sin(2*khinge*u**2) + (CD0 + CDmax + (CD0 - CDmax)*np.cos(2*khinge*u**2))*np.sin(phi)))/(2.*mb) , 0., 0., 0.],
        [0., 0.,      0.,     (kaero*u*((CD0 + CDmax + (CD0 - CDmax)*np.cos(2*khinge*u**2))*np.cos(phi) - 2*CLmax*np.sin(2*khinge*u**2)*np.sin(phi)))/(2.*mb), 0., 0., 0.],
        [0., 0.,      0.,         0.   , 0., 0., 0. ]
    ])

    Bd = np.matrix([
        [1.],
        [0.],
        [0.],
        [0.],
        [(kaero*(np.cos(phi)*(CD0 + CDmax + (CD0 - CDmax)*np.cos(2*khinge*u**2) - 
          4*(CD0 - CDmax)*khinge*u**2*np.sin(2*khinge*u**2)) - 
       2*CLmax*(4*khinge*u**2*np.cos(2*khinge*u**2) + 
          np.sin(2*khinge*u**2))*np.sin(phi)))/(2.*mb)],
        [(kaero*(2*CLmax*np.cos(phi)*np.sin(2*khinge*u**2) + 
       (CD0 + CDmax - 4*(CD0 - CDmax)*khinge*u**2*
           np.sin(2*khinge*u**2))*np.sin(phi) + 
       np.cos(2*khinge*u**2)*
        (8*CLmax*khinge*u**2*np.cos(phi) + (CD0 - CDmax)*np.sin(phi))))/
   (2.*mb)],
        [(kaero*(-(d*(CD0 + CDmax + (CD0 - CDmax)*np.cos(2*khinge*u**2) - 
            4*(CD0 - CDmax)*khinge*u**2*np.sin(2*khinge*u**2))) + 
       2*CLmax*(4*khinge*u**2*np.cos(2*khinge*u**2) + 
          np.sin(2*khinge*u**2))*sigma))/(2.*ib)]
        ])

    return Ad, Bd
# [nx, nu] = Bd.shape
nx = 7
nu = 1

# MPC stuff -----------

# Constraints
u0 = 0
umin = -np.inf * np.ones((1))
umax = np.inf * np.ones((1))
xmin = -np.inf * np.ones((nx))
xmax = np.inf * np.ones((nx))

# Objective function
Q = np.diag([0., 10., 10., 10., 0., 0., 0.])
QN = Q
R = 0.1

# Initial and reference states
y0 = np.zeros(nx)
tvec = np.arange(0, tf, dt)
y = np.zeros((len(tvec), nx))
y[0, :] = y0
yr = np.array([0.,0., 1.,0.,0.,0.,0.])

# Prediction horizon
N = 10

# Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
# - quadratic objective
P = sp.linalg.block_diag([np.kron(np.eye(N), Q), QN,
                       np.kron(np.eye(N), R)])
# - linear objective
q = np.hstack([np.kron(np.ones(N), -Q.dot(yr)), -QN.dot(yr),
               np.zeros(N*nu)])

# - linear dynamics


leq = np.hstack([-y0, np.zeros(N*nx)])
ueq = leq
# - input and state constraints
Aineq = np.eye((N+1)*nx + N*nu)
lineq = np.hstack([np.kron(np.ones(N+1), xmin), np.kron(np.ones(N), umin)])
uineq = np.hstack([np.kron(np.ones(N+1), xmax), np.kron(np.ones(N), umax)])
# - OSQP constraints
l = np.hstack([leq, lineq])
u = np.hstack([ueq, uineq])

# Create an OSQP object
prob = osqp.OSQP()
# Setup workspace
Ad, Bd = getLin(0, 0, 0)
# LTI version: identical Ad, Bd
Ax = np.kron(np.eye(N+1),-np.eye(nx)) + np.kron(np.eye(N+1, k=-1), Ad)
Bu = np.kron(np.vstack([np.zeros((1, N)), np.eye(N)]), Bd)
# print(Ax.shape)
# print(Bu.shape)
Aeq = np.hstack([Ax, Bu])
A = np.vstack([Aeq, Aineq])

prob.setup(P, q, A, l, u, warm_start=True)

# LTV version: need to generate new A,B as time goes on

for ti in range(len(tvec)):
    if ti == 1:
        continue

    # need u(t)
    t = tvec[ti]
    unom = np.array([15 * np.sin(2 * np.pi * 170 * t)])
    # need sigma(t)
    sigma = y[ti - 1, 0]
    phi = y[ti - 1, 3]
    Ad, Bd = getLin(unom, sigma, phi)	
    
    # LTV set different Ad according to a nominal stroke pattern
    # FIXME make 
    Ax = np.kron(np.eye(N+1),-np.eye(nx)) + np.kron(np.eye(N+1, k=-1), Ad)
    Bu = np.kron(np.vstack([np.zeros((1, N)), np.eye(N)]), Bd)
    Aeq = np.hstack([Ax, Bu])
    A = np.vstack([Aeq, Aineq])
    # print(A.shape)

    # Need to update the dynamics
    prob.update(Ax=A)

    u = unom

    y[ti, :] = y[ti - 1, :] + (Ad.dot(y[ti - 1, :]) + Bd * u) * dt


# PLOT ---

plt.figure()
np = 3

plt.subplot(np, 1, 1)
plt.plot(tvec, y[:, 0], label='sigma')
plt.legend()
plt.title('timestep = '+str(dt))

plt.subplot(np, 1, 2)
plt.plot(y[:, 1], y[:, 2], label='xz')
plt.legend()

plt.subplot(np, 1, 3)
plt.plot(tvec, y[:, 3], label='phi')
plt.legend()

plt.show()

# sys.exit(0)
