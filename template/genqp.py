import osqp
import numpy as np
import scipy.sparse as sp
from scipy.spatial.transform import Rotation
np.set_printoptions(precision=4, suppress=True, linewidth=200)

def qpSetupDense(n, m):
    "Helper to set up a dense QP"
    model = osqp.OSQP()
    P = sp.csc_matrix(np.ones((n,n)))
    q = np.zeros(n)
    l = -np.inf * np.ones(m)
    u = np.inf * np.ones(m)
    A = sp.csc_matrix(np.ones((m,n)))
    model.setup(P=P, q=q, A=A, l=l, u=u, eps_rel=1e-4, eps_abs=1e-4, verbose=False)
    return model

def mpcDirtran(m, N, dt, snom, y0, Qfdiag, ydes, gms, umin, umax):
    """See https://github.com/avikde/robobee3d/pull/181.
    snom should be an N, shaped array"""
    nq = 6 #q = (p,s)
    nu = 3
    ny = 2*nq
    # For dirtran
    nx = N * (ny + nu)
    assert len(snom) == N
    
    Ad = np.eye(ny)
    Ad[:nq, nq:] = dt * np.eye(nq)
    # B(s0) function
    Bs = lambda s : np.block([
        [1/m * np.asarray(np.reshape(s, (3,1))), np.zeros((3,2))], 
        [np.zeros((2,1)), np.eye(2)], 
        [np.zeros((1,1)), np.zeros((1,2))]])
    Bds = lambda s : np.vstack((np.zeros((nq,nu)), Bs(s))) * dt
    cd = dt * np.reshape(np.array([0, 0, -gms, 0, 0, 0]), nq, 1)

    # Construct dynamics constraint
    A = np.zeros((N*ny + nu, nx))
    l = np.zeros(N*ny + nu)
    u = np.zeros(N*ny + nu)
    for k in range(N):
        # x(k+1) = Ad*xk + Bd(sk)*uk
        A[k*ny:(k+1)*ny, k*ny:(k+1)*ny] = -np.eye(ny) # for -x1...xN+1 on the LHS
        if k > 0:
            A[k*ny:(k+1)*ny, (k-1)*ny:(k)*ny] = Ad # for Ad*x(k-1)
        A[k*ny:(k+1)*ny, (N*ny + k*nu):(N*ny + (k+1)*nu)] = Bds(snom[k]) # Bd(sk)

        l[k*ny:(k+1)*ny] = u[k*ny:(k+1)*ny] = -cd
        # only in the first eqn
        if k == 0:
            l[k*ny:(k+1)*ny] += -Ad @ np.asarray(y0)
            u[k*ny:(k+1)*ny] += -Ad @ np.asarray(y0)
    # Input limits
    A[N*ny:,-nu:] = np.eye(nu)
    l[-nu:] = np.tile(umin, (N,1))
    A = sp.csc_matrix(A)
    # print(A, c)

    # construct objective. 
    P = np.zeros((nx, nx))
    q = np.zeros(nx)
    # only final cost
    for i in range(ny):
        P[nx - ny + i, nx - ny + i] = Qfdiag[i]
    P = sp.csc_matrix(P)
    q[-ny:] = -np.asarray(ydes)

    # osqp
    model = osqp.OSQP()
    model.setup(P=P, q=q, A=A, l=l, u=u, eps_rel=1e-4, eps_abs=1e-4, verbose=False)
    return model

if __name__ == "__main__":
    # # WLQP gen
    # prob = qpSetupDense(4,4)

    prob = mpcDirtran(100, 3, 2, [[0.1, 0.1, 1], [0.2, 0.1, 1], [0.3, 0.1, 1], ], [1, 0.2, 0.1, 0.1, 0.2, 0.9, 0, 0, 0, 0, 0, 0], [1, 1, 1, 0.1, 0.1, 0.1, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3], [2, 0.2, 0.1, 0, 0, 1, 0, 0, 0, 0, 0, 0], -9.81e-5, [-100, -100, -100], [100, 100, 100])
    
    # # codegen
    # try:
    #     prob.codegen('gen', project_type='', force_rewrite=True, parameters='matrices', FLOAT=True, LONG=False)
    # except:
    #     # No worries if python module failed to compile
    #     pass
