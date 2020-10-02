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

class UprightMPC:
    def __init__(self, m, N, dt, snom, y0, Qfdiag, ydes, gms, umin, umax):
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
            [1/m * np.reshape(s, (3,1)), np.zeros((3,2))], 
            [np.zeros((2,1)), np.eye(2)], 
            [np.zeros((1,1)), np.zeros((1,2))]])
        Bds = lambda s : np.vstack((np.zeros((nq,nu)), Bs(s))) * dt
        cd = dt * np.hstack((np.zeros(6), np.array([0, 0, -gms, 0, 0, 0])))

        # Construct dynamics constraint
        self.A = np.zeros((N*ny + N*nu, nx))
        self.l = np.zeros(N*ny + N*nu)
        self.u = np.zeros(N*ny + N*nu)
        for k in range(N):
            # x(k+1) = Ad*xk + Bd(sk)*uk
            self.A[k*ny:(k+1)*ny, k*ny:(k+1)*ny] = -np.eye(ny) # for -x1...xN+1 on the LHS
            if k > 0:
                self.A[k*ny:(k+1)*ny, (k-1)*ny:(k)*ny] = Ad # for Ad*x(k-1)
            self.A[k*ny:(k+1)*ny, (N*ny + k*nu):(N*ny + (k+1)*nu)] = Bds(snom[k]) # Bd(sk)

            self.l[k*ny:(k+1)*ny] = self.u[k*ny:(k+1)*ny] = -cd
            # only in the first eqn
            if k == 0:
                self.l[k*ny:(k+1)*ny] += -Ad @ np.asarray(y0)
                self.u[k*ny:(k+1)*ny] += -Ad @ np.asarray(y0)
        # Input limits
        self.A[N*ny:,-N*nu:] = np.eye(N*nu)
        self.l[-N*nu:] = np.tile(umin, N)
        self.u[-N*nu:] = np.tile(umax, N)
        self.A = sp.csc_matrix(self.A)
        # print(A, c)

        # construct objective. 
        self.P = np.zeros((nx, nx))
        self.q = np.zeros(nx)
        # only final cost
        for i in range(ny):
            self.P[nx - ny + i, nx - ny + i] = Qfdiag[i]
        self.P = sp.csc_matrix(self.P)
        self.q[-ny:] = -np.asarray(ydes)
        
    def toOSQP(self):
        # osqp
        model = osqp.OSQP()
        model.setup(P=self.P, q=self.q, A=self.A, l=self.l, u=self.u, eps_rel=1e-4, eps_abs=1e-4, verbose=False)
        return model

if __name__ == "__main__":
    # # WLQP gen
    # prob = qpSetupDense(4,4)

    up = UprightMPC(100, 3, 2, [[0.1, 0.1, 1], [0.2, 0.1, 1], [0.3, 0.1, 1], ], [1, 0.2, 0.1, 0.1, 0.2, 0.9, 0, 0, 0, 0, 0, 0], [1, 1, 1, 0.1, 0.1, 0.1, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3], [2, 0.2, 0.1, 0, 0, 1, 0, 0, 0, 0, 0, 0], -9.81e-5, [-100, -100, -100], [100, 100, 100])
    
    # # codegen
    # try:
    #     prob.codegen('gen', project_type='', force_rewrite=True, parameters='matrices', FLOAT=True, LONG=False)
    # except:
    #     # No worries if python module failed to compile
    #     pass
