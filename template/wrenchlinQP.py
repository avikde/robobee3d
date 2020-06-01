import osqp
import autograd.numpy as np
from autograd import jacobian
import scipy.sparse as sp
from scipy.spatial.transform import Rotation

# TODO: this can be factorized
from ca6dynamics import dynamicsTerms, wrenchMap
dw_du = jacobian(wrenchMap)

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
    
    def update(self, p0, h0, B0, w0, u0, dt, Qd, pdes):
        A = np.eye(self.n)
        # general form
        a0 = p0 - dt * h0 + dt * B0 @ w0
        A1 = dt * B0 @ dw_du(u0)

        P = A1.T @ np.diag(Qd) @ A1
        q = A1.T @ np.diag(Qd) @ (a0 - pdes)
        u = 100 * np.ones(self.n)
        l = -100 * np.ones(self.n)
        # update OSQP
        Px = P[np.tril_indices(P.shape[0])] # need in col order
        self.model.update(Px=Px, q=q, l=l, u=u, Ax=np.ravel(A))
        res = self.model.solve()
        return res.x + u0

    def test(self):
        dq = np.zeros(6)
        M0, h0, B0 = dynamicsTerms(np.zeros(3), Rotation.from_euler('x', 0), dq)
        p0 = M0 @ dq
        u0 = np.array([1.0,0.0,0.0,1.0,0.0,0.0])
        dt = 2
        Qd = 0.1 * np.ones(6)
        pdes = np.array([0,0,10,0,0,0])
        w0 = np.zeros(6)
        print(self.update(p0, h0, B0, w0, u0, dt, Qd, pdes))
        
