import osqp
import numpy as np
import scipy as sp
import scipy.linalg
import scipy.sparse as sparse
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

# [nx, nu] = Bd.shape
nx = 7
nu = 2

# ---- Original example https://osqp.org/docs/examples/mpc.html

# Discrete time model
def getLin(u0, tf0, y):
    # q = (x,z,phi), y = (sigma, q, dq)
    sigma0 = y[0]
    phi = y[3]
    dphi = y[6]
    dx = y[4]
    dz = y[5]
    g = 9.81
    d = 2
    
    tf02 = tf0**2
    u03 = u0**3
    u02 = u0**2
    cphi = np.cos(phi)
    sphi = np.sin(phi)
    pi = np.pi
    kaero2 = 0.091875
    su0 = np.sign(u0)

    # Copied from mathematica 

    Ad = np.matrix([
        [1,0,0,0,0,0,0],
        [0,1,0,0,tf0,0,0],
        [0,0,1,0,0,tf0,0],
   [0,0,0,1,0,0,tf0],
     [0,0,0,(tf0*(-(cphi*kaero2*u02) + kaero2*sphi*su0*u02))/100.,1,0,0],
        [0,0,0,(tf0*(-(kaero2*sphi*u02) - cphi*kaero2*su0*u02))/100.,0,1,0],
        [(kaero2*tf0*u02)/7200.,0,0,0,0,0,1]
    ])

    Bd = np.matrix([
        [tf0,u0],
        [0,dx],
        [0,dz],
        [0,dphi],
        [(tf0*(-2*kaero2*sphi*u0 - \
            2*cphi*kaero2*su0*u0))/100.,(-(kaero2*sphi*u02) - \
            cphi*kaero2*su0*u02)/100.],
        [(tf0*(2*cphi*kaero2*u0 - \
            2*kaero2*sphi*su0*u0))/100.,(-0.9800000000000001 + cphi*kaero2*u02 - \
            kaero2*sphi*su0*u02)/100.],
        [(tf0*(2*d*kaero2*su0*u0 + 2*kaero2*u0*(sigma0 + (tf0*u0)/2.) + \
            (kaero2*tf0*u02)/2.))/7200.,(d*kaero2*su0*u02 + kaero2*(sigma0 + \
            (tf0*u0)/2.)*u02)/7200. + (kaero2*tf0*u03)/14400.]
        ])

    return Ad, Bd

def getCondensed(u0, tf0, y0, l, u):
    # LTV set different Ad according to a nominal stroke pattern

    # Setup matrices
    Ax = np.kron(np.eye(N+1),-np.eye(nx))# + np.kron(np.eye(N+1, k=-1), Ad)
    # Will replace these below
    Bu = np.kron(np.vstack([np.zeros((1, N)), np.eye(N)]), np.zeros((nx, nu)))

    # Simulate forward and linearize along that
    uprev = u0 # TODO: consider using the previous solution (scott suggestion)
    for j in range(N):
        # For lin provide different A, B for the blocks (not kron)
        Ad, Bd = getLin(uprev, tf0, y0)
        # Change y0
        y0 = Ad.dot(y0) + Bd.dot(np.array([uprev,tf0]))
        
        Ax[nx*(j+1):nx*(j+2), nx*j:nx*(j+1)] = Ad
        Bu[nx*(j+1):nx*(j+2), nu*j:nu*(j+1)] = Bd
        
        # print(Ad)
        # print(Ax[nx*(0+1):nx*(0+2), nx*0:nx*(0+1)])
        # print(Bd)
        # print(Bu[nx*(0+1):nx*(0+2), nu*0:nu*(0+1)])
        # sys.exit(0)
        uprev = -uprev # alternate u sign

    # stack
    Aeq = np.hstack([Ax, Bu])
    A = np.vstack([Aeq, Aineq])
    # sys.exit(0)
    # print(A.shape)

    # Update initial state
    l[:nx] = -ympc[ti-1, :]
    u[:nx] = -ympc[ti-1, :]

    return A, l, u


# MPC stuff -----------

# Constraints
# u0 = np.array([0,0])

umin = np.array([-10,0.1])
umax = np.array([10,100])
xmin = -np.inf * np.ones((nx))
xmax = np.inf * np.ones((nx))

# Objective function
Q = np.diag([0., 10., 10., 10., 0., 0., 0.])
QN = Q
R = 10.0 * np.eye(2)

# Initial and reference states
y0 = np.zeros(nx)
y0[0] = -24.5
yr = np.array([0.,0., 0.,0.,0.,0.,0.])
tvec = range(30)
# Arrays for logging
y = np.zeros((len(tvec), nx))
ympc = np.zeros_like(y)
y[0, :] = ympc[0,:] = y0

# Prediction horizon
N = 1 # FIXME: testing

# Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
# - quadratic objective
P = sparse.block_diag([np.kron(np.eye(N), Q), QN, np.kron(np.eye(N), R)])
# print(np.kron(np.eye(N), Q).shape)
# print(QN.shape)
# print(np.kron(np.eye(N), R).shape)
# print(Pp.shape)
# - linear objective
q = np.hstack([np.kron(np.ones(N), -Q.dot(yr)), -QN.dot(yr),
               np.zeros(N*nu)])

# - linear dynamics
# print(getLin(10, 5, np.zeros((7)))[0])
# print(getLin(-10, 5, np.zeros((7)))[0])
# print(getLin(10, 5, np.zeros((7)))[1])
# print(getLin(-10, 5, np.zeros((7)))[1])
# sys.exit(0)

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
Ad, Bd = getLin(0, 0, np.zeros((nx)))

# LTI version: identical Ad, Bd
Ax = np.kron(np.eye(N+1),-np.eye(nx)) + np.kron(np.eye(N+1, k=-1), Ad)
Bu = np.kron(np.vstack([np.zeros((1, N)), np.eye(N)]), Bd)
# print(Ax.shape)
# print(Bu.shape)
Aeq = np.hstack([Ax, Bu])
A = sparse.vstack([Aeq, Aineq])


prob.setup(P, q, A, l, u, warm_start=True, eps_abs=1.0e-02, eps_rel=1.0e-02)

# LTV version: need to generate new A,B as time goes on
u0 = 5
tf0 = 5
umpc = np.array([u0, tf0]) # keep track of previous input

for ti in tvec:
    # Alternate the nominal input
    u0 = -u0

    if ti == 0:
        continue

    # need u(t)
    # t = tvec[ti]
    # unom = np.array([15 * np.sin(2 * np.pi * 170 * t)])

    Ad0, Bd0 = getLin(u0, tf0, y[ti-1,:])
    uvecnom = np.array([u0, tf0])
    # Simulate OL system
    y[ti, :] = Ad0.dot(y[ti-1, :]) + Bd0.dot(uvecnom)

    # MPC ---
    # Need to update the dynamics
    A, l, u = getCondensed(u0, tf0, ympc[ti-1,:], l, u)
    prob.update(Ax=A, l=l, u=u)

    # Solve
    res = prob.solve()
    # Check solver status
    if res.info.status != 'solved':
        # FIXME: says primal infeasible. try to print out objective and constraints for the OL solution
        # Use exactly the A, l, u, and the nominal input
        Ad = A[nx*(0+1):nx*(0+2), nx*0:nx*(0+1)]
        Bd = A[nx*(0+1):nx*(0+2), nx*(N+1) + nu*0:nx*(N+1) + nu*(0+1)]
        ympcNextWithNomInput = Ad.dot(ympc[ti-1,:]) + Bd.dot(uvecnom)
        xtest = np.hstack((ympc[ti-1,:], ympcNextWithNomInput, uvecnom))
        # print(xnom)
        # print(A.shape)
        # print(xnom.shape)
        test1 = A.dot(xtest)
        print('Worst violations: ', np.min(test1 - l), np.argmin(test1 - l), np.min(u - test1), np.argmin(u - test1))
        # these are dx, dphi
        raise ValueError('OSQP did not solve the problem!')
    
    # Apply first control input to the plant
    umpc = res.x[-N*nu:-(N-1)*nu]
    print(res.x)

    # Simulate
    # Ad1, Bd1 = getLin(umpc[0], umpc[1], ympc[ti-1,:])
    # ympc[ti, :] = Ad1.dot(ympc[ti-1, :]) + Bd1.dot(umpc)


# PLOT ---

plt.figure()
np = 3

plt.subplot(np, 1, 1)
plt.plot(tvec, y[:, 0], label='ol')
plt.plot(tvec, ympc[:, 0], label='mpc')
plt.ylabel('sigma')
plt.legend()
# plt.title('timestep = '+str(dt))

plt.subplot(np, 1, 2)
plt.plot(y[:, 1], y[:, 2], label='ol')
plt.plot(ympc[:, 1], ympc[:, 2], label='mpc')
plt.ylabel('xz traj')
plt.legend()

plt.subplot(np, 1, 3)
plt.plot(tvec, y[:, 3], label='ol')
plt.plot(tvec, ympc[:, 3], label='ol')
plt.ylabel('phi')
plt.legend()

plt.show()

# sys.exit(0)
