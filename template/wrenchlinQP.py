import osqp
import autograd.numpy as np
from autograd import jacobian
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

class WrenchLinQP(object):
    def __init__(self, n, m, dynamicsTerms, wrenchMap, dwduMap=None, u0=None, dumax=None):
        self.model = qpSetupDense(n, m)
        self.n = n
        self.m = m
        self.u0 = np.zeros(n) if u0 is None else np.array(u0)
        self.dynamicsTerms = dynamicsTerms
        self.wrenchMap = wrenchMap
        self.w0 = self.wrenchMap(self.u0)
        if dwduMap is None:
            # Use autograd
            self.dwduMap = jacobian(self.wrenchMap)
        else:
            # Use the provided Jacobian function
            self.dwduMap = dwduMap
        # How much u can change in one step
        if dumax is None:
            self.U = 1e-2 * np.ones(self.n)
        else:
            self.U = np.array(dumax)
    
    def update(self, p0, h0, B0, dt, Qd, pdes):
        if self.n >= 4:
            # assume ca6/sdab model
            curDwDu = self.dwduMap(self.u0)
            # print(curDwDu)
            # FIXME: made up Jac for testing. out = 6
            curDwDu = np.array([
                [0,0,0,0],
                [0,0,0,0],
                [0.03,0,0,0],
                [0,0,0,0],
                [0,0,0,0],
                [0,0,0,0]
                ])
        else:
            # For testing other models (assume w=u if n != 6)
            curDwDu = np.eye(self.n)
            pdes = pdes[2:3] # assume z component

        # general form
        A = np.eye(self.n)
        a0 = p0 - dt * h0 + dt * B0 @ self.w0
        # print(curDwDu)
        A1 = dt * B0 @ curDwDu
        P = A1.T @ np.diag(Qd) @ A1
        q = A1.T @ np.diag(Qd) @ (a0 - pdes)
        # print("p0=", a0, "pdes=", pdes)
        L = -self.U
        # update OSQP
        Px = P[np.tril_indices(P.shape[0])] # need in col order
        self.model.update(Px=Px, q=q, l=L, u=self.U, Ax=np.ravel(A))
        res = self.model.solve()
        self.u0 = res.x + self.u0
        # FIXME: 
        self.u0[1:].fill(0)

        if self.n >= 4:
            self.w0 = self.wrenchMap(self.u0)
            print(res.x, self.w0[[2,4]], self.u0)
            return self.u0
        else:
            self.w0 = self.u0
            r = np.zeros(6)
            r[0] = r[3] = self.u0[0]
            return r
    
    def updateFromState(self, t, q, dq, pdes, dt=2):
        M0, h0, B0 = self.dynamicsTerms(q, dq) # B=I, since the input is B*w(u)
        p0 = M0 @ dq
        Qd = np.hstack((1.0*np.ones(3), 0.1*np.ones(3)))

        if self.n < 4:
            # For testing other models (assume w=u if n<4)
            p0 = np.array([M0[2,2] * dq[2]]) # assume z
            Qd = Qd[2:3]
            B0 = np.eye(1)
            h0 = np.zeros(1)

        return self.update(p0, h0, B0, dt, Qd, pdes)

    def test(self):
        q = [0.,0,0,0,0,0,1]
        dq = np.zeros(6)
        self.uprev = np.array([1.0,0.0,0.0,1.0,0.0,0.0])
        pdes = np.array([0,0,10,0,0,0])
        print(self.updateFromState(0., q, dq, pdes))
        
