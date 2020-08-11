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
    def __init__(self, n, m, dynamicsTerms, wrenchMap, dwduMap=None, u0=None, dumax=None, umin=None, umax=None):
        self.model = qpSetupDense(n, m)
        self.n = n
        self.m = m
        self.u0 = np.zeros(n) if u0 is None else np.array(u0)
        self.dynamicsTerms = dynamicsTerms
        self.wrenchMap = wrenchMap
        self.w0 = self.wrenchMap(self.u0)
        self.umin = umin
        self.umax = umax
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
    
    def update(self, h0, Qd, pdotdes):
        """See https://github.com/avikde/robobee3d/pull/166"""
        if self.n >= 4:
            # assume ca6/sdab model
            curDwDu = self.dwduMap(self.u0)
        else:
            # For testing other models (assume w=u if n != 6)
            curDwDu = np.eye(self.n)
            pdes = pdes[2:3] # assume z component

        # general form
        A = np.eye(self.n) # for input rate limit

        # Objective: momentum reference dynamics
        a0 = self.w0 - h0 - pdotdes
        # print(curDwDu)
        A1 = curDwDu
        P = A1.T @ np.diag(Qd) @ A1
        q = A1.T @ np.diag(Qd) @ a0
        L = -self.U
        U = np.copy(self.U)

        # Input limits (not just rate limits)
        if self.umin is not None:
            for i in range(len(self.umin)):
                if self.u0[i] < self.umin[i]:
                    L[i] = 0 # do not reduce further
                elif self.u0[i] > self.umax[i]:
                    U[i] = 0 # do not increase further

        # update OSQP
        Px = P[np.tril_indices(P.shape[0])] # need in col order
        self.model.update(Px=Px, q=q, l=L, u=U, Ax=np.ravel(A))
        res = self.model.solve()
        self.u0 = res.x + self.u0
        # self.u0[0] = 110 + 3 * (pdes[2] - p0[2])

        if self.n >= 4:
            self.w0 = self.wrenchMap(self.u0)
            return self.u0
        else:
            self.w0 = self.u0
            r = np.zeros(6)
            r[0] = r[3] = self.u0[0]
            return r
    
    def updateFromState(self, t, q, dq, pdes, kpmom=np.array([0,0,1,0.1,0.1,0.1])):
        M0, h0 = self.dynamicsTerms(q, dq) # B=I, since the input is B*w(u)
        p0 = M0 @ dq
        Qd = np.hstack((1.0*np.ones(3), 0.1*np.ones(3)))

        # Momentum reference dynamics https://github.com/avikde/robobee3d/pull/166 TODO: incorporate as MPC
        pdotdes = kpmom * (pdes - p0)

        if self.n < 4:
            # For testing other models (assume w=u if n<4)
            p0 = np.array([M0[2,2] * dq[2]]) # assume z
            Qd = Qd[2:3]
            h0 = np.zeros(1)

        return self.update(h0, Qd, pdotdes)

    def test(self):
        q = [0.,0,0,0,0,0,1]
        dq = np.zeros(6)
        self.uprev = np.array([1.0,0.0,0.0,1.0,0.0,0.0])
        pdes = np.array([0,0,10,0,0,0])
        print(self.updateFromState(0., q, dq, pdes))
        
if __name__ == "__main__":
    # Autogen code
    prob = qpSetupDense(4,4)
    # codegen
    try:
        prob.codegen('gen', project_type='', force_rewrite=True, parameters='matrices', FLOAT=True, LONG=False)
    except:
        # No worries if python module failed to compile
        pass


