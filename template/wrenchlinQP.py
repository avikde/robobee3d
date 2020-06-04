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
        self.uprev = np.zeros(n)
    
    def update(self, p0, h0, B0, w0, dt, Qd, pdes):
        A = np.eye(self.n)
        # general form
        a0 = p0 - dt * h0 + dt * B0 @ w0
        A1 = dt * B0 @ dw_du(self.uprev)

        P = A1.T @ np.diag(Qd) @ A1
        q = A1.T @ np.diag(Qd) @ (a0 - pdes)
        u = np.array([1e-3,1e-3,1e-3,1e-3,1e-3,1e-3])
        l = -u
        # update OSQP
        Px = P[np.tril_indices(P.shape[0])] # need in col order
        self.model.update(Px=Px, q=q, l=l, u=u, Ax=np.ravel(A))
        res = self.model.solve()
        self.uprev = res.x + self.uprev
        return self.uprev
    
    def updateFromState(self, Rb, dq, pdes):
        M0, h0, B0 = dynamicsTerms(Rb, dq)
        p0 = M0 @ dq
        dt = 2
        Qd = np.array([0.1,0.1,0.1,100,100,100])
        w0 = np.zeros(6)
        return self.update(p0, h0, B0, w0, dt, Qd, pdes)

    def test(self):
        Rb = Rotation.from_euler('x', 0)
        dq = np.zeros(6)
        self.uprev = np.array([1.0,0.0,0.0,1.0,0.0,0.0])
        pdes = np.array([0,0,10,0,0,0])
        print(self.updateFromState(Rb, dq, pdes))
        
