import osqp
import numpy as np
import scipy.sparse as sp

def qpSetupDense(n, m):
    "Helper to set up a dense QP"
    model = osqp.OSQP()
    P = sp.csc_matrix(np.ones((n,n)))
    q = np.zeros(n)
    l = -np.inf * np.ones(m)
    u = np.inf * np.ones(m)
    A = sp.csc_matrix(np.ones((m,n)))
    model.setup(P=P, q=q, A=A, l=l, u=u, eps_rel=1e-2, eps_abs=1e-2, verbose=False)
    return model

class WrenchLinQP(object):
    def __init__(self, n, m):
        self.model = qpSetupDense(n, m)
        self.n = n
        self.m = m
    
    def update(self, p0, h0, B0, w0, dw_du0, dt, Qd, pdes):
        A = np.eye(self.n)
        # general form
        a0 = p0 - dt * h0 + dt * B0 @ w0
        A1 = dt * B0 @ dw_du0

        P = A1.T @ np.diag(Qd) @ A1
        q = A1.T @ np.diag(Qd) @ (a0 - pdes)
        u = 100 * np.ones(self.n)
        l = -100 * np.ones(self.n)
        # update OSQP
        Px = P[np.tril_indices(P.shape[0])] # need in col order
        self.model.update(Px=Px, q=q, l=l, u=u, Ax=np.ravel(A))
        res = self.model.solve()
        return res.x
        
