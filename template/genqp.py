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
    nq = 6 #q = (p,s)
    nu = 3
    ny = 2*nq
    
    cd = lambda self, g, m: dt * np.hstack((np.zeros(6), np.array([0, 0, -g/m, 0, 0, 0])))

    def __init__(self, N, dt, snom, y0, Qfdiag, ydes, g, m, ms, umin, umax):
        """See https://github.com/avikde/robobee3d/pull/181.
        snom should be an N, shaped array"""
        # For dirtran
        self.nx = N * (self.ny + self.nu)
        assert len(snom) == N

        # # dummy values: will be updated in the QP
        # dt = 2
        
        Ad = np.eye(self.ny)
        Ad[:self.nq, self.nq:] = dt * np.eye(self.nq)
        # B(s0) function
        Bs = lambda s : np.block([
            [1/m * np.reshape(s, (3,1)), np.zeros((3,2))], 
            [np.zeros((2,1)), 1/ms * np.eye(2)], 
            [np.zeros((1,1)), np.zeros((1,2))]])
        Bds = lambda s : np.vstack((np.zeros((self.nq,self.nu)), Bs(s))) * dt

        # Construct dynamics constraint
        self.A = np.zeros((self.nx, self.nx))
        self.l = np.zeros(self.nx)
        self.u = np.zeros(self.nx)
        for k in range(N):
            # x(k+1) = Ad*xk + Bd(sk)*uk
            self.A[k*self.ny:(k+1)*self.ny, k*self.ny:(k+1)*self.ny] = -np.eye(self.ny) # for -x1...xN+1 on the LHS
            if k > 0:
                self.A[k*self.ny:(k+1)*self.ny, (k-1)*self.ny:(k)*self.ny] = Ad # for Ad*x(k-1)
            self.A[k*self.ny:(k+1)*self.ny, (N*self.ny + k*self.nu):(N*self.ny + (k+1)*self.nu)] = Bds(snom[k]) # Bd(sk)

            self.l[k*self.ny:(k+1)*self.ny] = self.u[k*self.ny:(k+1)*self.ny] = -self.cd(g, m)
            # only in the first eqn
            if k == 0:
                self.l[k*self.ny:(k+1)*self.ny] += -Ad @ np.asarray(y0)
                self.u[k*self.ny:(k+1)*self.ny] += -Ad @ np.asarray(y0)
        # Input limits
        self.A[N*self.ny:,-N*self.nu:] = np.eye(N*self.nu)
        self.l[-N*self.nu:] = np.tile(umin, N)
        self.u[-N*self.nu:] = np.tile(umax, N)
        self.A = sp.csc_matrix(self.A)
        # print(A, c)

        # construct objective. 
        self.P = np.zeros((self.nx, self.nx))
        self.q = np.zeros(self.nx)
        # only final cost
        for i in range(self.ny):
            self.P[self.nx - self.ny + i, self.nx - self.ny + i] = Qfdiag[i]
        self.P = sp.csc_matrix(self.P)
        self.q[-self.ny:] = -np.asarray(ydes)
    
    def update(self, dt, snom, y0, Qfdiag, ydes, g, m, ms, umin, umax):
        # TODO: Axidx, Pxidx make in init
        # should not need the arguments in constructor (just sets sparsity)

        # update l, u
        self.l = np.hstack((np.tile(-self.cd(g, m), N), np.tile(umin, N)))
        self.u = np.hstack((np.tile(-self.cd(g, m), N), np.tile(umax, N)))
        y0pdt = np.asarray(y0)
        y0pdt[:self.nq] += dt * y0pdt[self.nq:]
        self.l[:self.ny] -= y0pdt
        self.u[:self.ny] -= y0pdt

        # To test dynamics constraint after need to update sparse csc_matrix
        print(self.A.nnz, len(self.A.data))
    
    def dynamics(self, y, u, dt, g, m, ms, s0):
        # Not needed for optimization, just to check
        q, dq = y[:6], y[6:12]
        # s = q[3:]
        uT = u[0] # thrust
        uM = u[1:] # moment
        ddq = np.hstack((-g/m * np.array([0,0,1]) + uT/m * np.asarray(s0), 1/ms * uM, 0))
        dq1 = dq + dt * ddq
        q1 = q + dt * dq
        return np.hstack((q1, dq1))

    def dynamicsTest(self, N, dt, g, m, ms, snom, y0):
        ys = np.zeros((N, self.ny))
        us = np.random.rand(N, self.nu)
        y0 = np.asarray(y0)

        for k in range(N):
            ys[k,:] = self.dynamics(y0, us[k,:], dt, g, m, ms, snom[k])
            y0 = ys[k,:]
        
        # reshape into an "x"
        xtest = np.hstack((np.ravel(ys), np.ravel(us)))
        # print(xtest.shape, up.nx)
        e1 = self.A @ xtest - self.l
        e2 = self.u - self.A @ xtest
        Z = np.zeros(N*self.ny)
        if np.allclose(e1[:N*self.ny], Z) and np.allclose(e2[:N*self.ny], Z):
            print('Dynamics test passed')
        else:
            print('Dynamics test FAILED', e1, e2)
        
    def toOSQP(self):
        # osqp
        model = osqp.OSQP()
        model.setup(P=self.P, q=self.q, A=self.A, l=self.l, u=self.u, eps_rel=1e-4, eps_abs=1e-4, verbose=False)
        return model

if __name__ == "__main__":
    # # WLQP gen
    # prob = qpSetupDense(4,4)
    N = 3
    dt = 2.0
    g = 9.81e-3
    m = 100.0
    ms = 1.0
    umax = np.array([100.0, 100.0, 100.0])
    umin = -umax
    snom = [[0.1, 0.1, 1], [0.2, 0.1, 1], [0.3, 0.1, 1]]
    y0 = [1, 0.2, 0.1, 0.1, 0.2, 0.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    Qfdiag = [1, 1, 1, 0.1, 0.1, 0.1, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3]
    ydes = [2, 0.2, 0.1, 0.4, 0.1, 1, 0.1, -0.1, -0.2, -0.3, -0.4, -0.5]

    up = UprightMPC(N, dt, snom, y0, Qfdiag, ydes, g, m, ms, umin, umax)

    up.update(dt, snom, y0, Qfdiag, ydes, g, m, ms, umin, umax)

    up.dynamicsTest(N, dt, g, m, ms, snom, y0)
    
    # # codegen
    # try:
    #     prob.codegen('gen', project_type='', force_rewrite=True, parameters='matrices', FLOAT=True, LONG=False)
    # except:
    #     # No worries if python module failed to compile
    #     pass
