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

def mpcDirtran(m, N, dt, snom, y0, Qfdiag, ydes):
    """See https://github.com/avikde/robobee3d/pull/181.
    snom should be an N, shaped array"""
    nq = 2 #q = (p,s)
    nu = 2
    ny = 2*nq
    # For dirtran
    nx = N * (ny + nu)
    assert len(snom) == N
    
    Ad = np.eye(ny)
    Ad[:nq, nq:] = dt * np.eye(nq)
    # B(s0) function
    Bs = lambda s : np.array([[1/m * s, 0], [0, 1]])
    Bds = lambda s : np.vstack((np.zeros((nq,nu)), Bs(s))) * dt

    # Construct dynamics constraint
    A = np.zeros((N*ny, nx))
    c = np.zeros(N*ny)
    for k in range(N):
        # x(k+1) = Ad*xk + Bd(sk)*uk
        A[k*ny:(k+1)*ny, k*ny:(k+1)*ny] = np.eye(ny) # for x1...xN+1 on the LHS
        if k > 0:
            A[k*ny:(k+1)*ny, (k-1)*ny:(k)*ny] = -Ad # for -Ad*x(k-1)
        A[k*ny:(k+1)*ny, (N*ny + k*nu):(N*ny + (k+1)*nu)] = -Bds(snom[k]) # -Bd(sk)
        # only in the first eqn
        if k == 0:
            c[k*ny:(k+1)*ny] = Ad @ np.asarray(y0)
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
    model.setup(P=P, q=q, A=A, l=c, u=c, eps_rel=1e-4, eps_abs=1e-4, verbose=False)
    return model

if __name__ == "__main__":
    # # WLQP gen
    # prob = qpSetupDense(4,4)

    prob = mpcDirtran(100, 3, 2, [0.1, 0.2, 0.3], [1, 0.1, 0, 0], [1, 0.1, 1e-3, 1e-3], [2, 0.2, 0, 0])
    
    # # codegen
    # try:
    #     prob.codegen('gen', project_type='', force_rewrite=True, parameters='matrices', FLOAT=True, LONG=False)
    # except:
    #     # No worries if python module failed to compile
    #     pass
