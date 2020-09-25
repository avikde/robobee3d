import osqp
import numpy as np
import scipy.sparse as sp
from scipy.spatial.transform import Rotation

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

def mpcDirtran(m, N, dt, snom):
    """See https://github.com/avikde/robobee3d/pull/181.
    snom should be an N, shaped array"""
    nq = 2 #q = (p,s)
    nu = 2
    assert len(snom) == N
    
    Ad = np.eye(nq*2)
    Ad[:nq, nq:] = dt * np.eye(nq)
    Ad = sp.csc_matrix(Ad)
    # B(s0) function
    Bs = lambda s : np.array([[1/m * s, 0], [0, 1]])
    Bds = lambda s : sp.csc_matrix(np.vstack((np.zeros((nq,nu)), Bs(s))))

    print(Ad, Bds(0.1))

if __name__ == "__main__":
    # # WLQP gen
    # prob = qpSetupDense(4,4)
    # # codegen
    # try:
    #     prob.codegen('gen', project_type='', force_rewrite=True, parameters='matrices', FLOAT=True, LONG=False)
    # except:
    #     # No worries if python module failed to compile
    #     pass

    mpcDirtran(100, 3, 2, [0.1, 0.2, 0.3])
