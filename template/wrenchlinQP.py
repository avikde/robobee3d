import osqp
import autograd.numpy as np
from autograd import jacobian
import scipy.sparse as sp
from scipy.spatial.transform import Rotation

# TODO: this can be factored
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
    model.setup(P=P, q=q, A=A, l=l, u=u, eps_rel=1e-4, eps_abs=1e-4, verbose=False)
    return model

class WrenchLinQP(object):
    def __init__(self, n, m):
        self.model = qpSetupDense(n, m)
        self.n = n
        self.m = m
        self.u0 = np.zeros(n)
        self.w0 = np.zeros(n)
    
    def update(self, p0, h0, B0, dt, Qd, pdes):
        if self.n == 6:
            # assume ca6 model
            u = 1e-2 * np.ones(self.n)
            curDwDu = dw_du(self.u0)
        else:
            # For testing other models (assume w=u if n != 6)
            u = np.ones(self.n) * 1e-1
            curDwDu = np.eye(self.n)
            pdes = pdes[2:3] # assume z component

        # general form
        A = np.eye(self.n)
        a0 = p0 - dt * h0 + dt * B0 @ self.w0
        # print(curDwDu)
        A1 = dt * B0 @ curDwDu
        P = A1.T @ np.diag(Qd) @ A1
        q = A1.T @ np.diag(Qd) @ (a0 - pdes)
        l = -u
        # update OSQP
        Px = P[np.tril_indices(P.shape[0])] # need in col order
        self.model.update(Px=Px, q=q, l=l, u=u, Ax=np.ravel(A))
        res = self.model.solve()
        self.u0 = res.x + self.u0

        if self.n == 6:
            self.w0 = wrenchMap(self.u0)
            return self.u0
        else:
            self.w0 = self.u0
            r = np.zeros(6)
            r[0] = r[3] = self.u0[0]
            return r
    
    def updateFromState(self, Rb, dq, pdes):
        M0, h0, B0 = dynamicsTerms(Rb, dq)
        p0 = M0 @ dq
        dt = 2
        Qd = 0.1 * ones(self.n)

        if self.n != 6:
            # For testing other models (assume w=u if n != 6)
            p0 = np.array([M0[2,2] * dq[2]]) # assume z
            Qd = Qd[2:3]
            B0 = np.eye(1)
            h0 = np.zeros(1)

        return self.update(p0, h0, B0, dt, Qd, pdes)

    def test(self):
        Rb = Rotation.from_euler('x', 0)
        dq = np.zeros(6)
        self.uprev = np.array([1.0,0.0,0.0,1.0,0.0,0.0])
        pdes = np.array([0,0,10,0,0,0])
        print(self.updateFromState(Rb, dq, pdes))
        
