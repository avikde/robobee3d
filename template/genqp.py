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

    def __init__(self, N):
        """See https://github.com/avikde/robobee3d/pull/181.
        snom should be an N, shaped array"""
        # For dirtran
        self.N = N
        self.nx = N * (self.nq + self.nu)

        # dummy values: will be updated in the QP
        dt = 2
        snom = [np.ones(3) for k in range(N)]
        smin = -np.ones(3)
        smax = np.ones(3)
        q0 = np.ones(self.nq)
        Qfdiag = np.ones(self.nq)
        Rdiag = np.ones(self.nu)
        qdes = -np.ones(self.nq)

        # B(s0) function
        Bs = lambda s : np.block([
            [np.reshape(s, (3,1)), np.zeros((3,2))], 
            [np.zeros((2,1)), np.eye(2)], 
            [np.zeros((1,1)), -s[:2]/s[2]]]) * dt

        # Construct dynamics constraint
        ncon = N * (self.nq + 3)
        self.A = np.zeros((ncon, self.nx))
        self.l = np.zeros(ncon)
        self.u = np.zeros(ncon)
        for k in range(N):
            # q(k+1) = q(k) + Bd(sk)*vk
            self.A[k*self.nq:(k+1)*self.nq, k*self.nq:(k+1)*self.nq] = -np.eye(self.nq) # for -q1...qN on the LHS
            if k > 0:
                self.A[k*self.nq:(k+1)*self.nq, (k-1)*self.nq:(k)*self.nq] = np.eye(self.nq) # for q
            self.A[k*self.nq:(k+1)*self.nq, (N*self.nq + k*self.nu):(N*self.nq + (k+1)*self.nu)] = Bs(snom[k]) # Bd(sk)

            # only in the first eqn
            if k == 0:
                self.l[k*self.nq:(k+1)*self.nq] += -np.asarray(q0)
                self.u[k*self.nq:(k+1)*self.nq] += -np.asarray(q0)

            # s limits
            self.A[(N*self.nq+3*k):(N*self.nq+(k+1)*3), (k)*self.nq:(k)*self.nq+3] = np.eye(3)
            self.l[(N*self.nq+3*k):(N*self.nq+(k+1)*3)] = smin
            self.u[(N*self.nq+3*k):(N*self.nq+(k+1)*3)] = smax
        self.A = sp.csc_matrix(self.A)
        # print(A, c)

        # construct objective. 
        self.q = np.zeros(self.nx)
        self.Pdiag = np.zeros(self.nx)
        # final cost
        self.Pdiag[(N-1)*self.nq : N*self.nq] = Qfdiag
        for k in range(N):
            self.Pdiag[N*self.nq + k*self.nu : N*self.nq + (k+1)*self.nu] = Rdiag
        self.q[(N-1)*self.nq:N*self.nq] = -(np.asarray(Qfdiag) * np.asarray(qdes))

        self.saveAxidx()
    
    def saveAxidx(self):
        A1nnz = (2*N-1) * self.nq + 3*N
        Bnnz = 7
        A2nnz = N * Bnnz

        # To test dynamics constraint after need to update sparse csc_matrix
        assert self.A.nnz == A1nnz + A2nnz

        # put together - store these. Should be filled with stacked Bi = dt*(sx,sy,sz,1,-sx/sz,1,-sy/sz)
        self.Axidx = range(A1nnz, A1nnz + A2nnz)
    
    def update(self, dt, snom, y0, Qfdiag, ydes, g, m, ms, umin, umax):
        # should not need the arguments in constructor (just sets sparsity)

        # update l, u
        self.l = np.hstack((np.tile(-self.cd(g, m), N), np.tile(umin, N)))
        self.u = np.hstack((np.tile(-self.cd(g, m), N), np.tile(umax, N)))
        y0pdt = np.asarray(np.copy(y0))
        y0pdt[:self.nq] += dt * y0pdt[self.nq:]
        self.l[:self.ny] -= y0pdt
        self.u[:self.ny] -= y0pdt

        # update q
        self.q[(N-1)*self.ny:N*self.ny] = -(np.asarray(Qfdiag) * np.asarray(ydes))

        # update P
        self.P.data = np.asarray(Qfdiag) # replace the whole thing

        # update A
        # print("hi",np.hstack(snom))
        self.A.data[self.Axidx] = np.hstack((np.full(self.AxidxNdt, dt), dt/m*np.hstack(snom), np.full(self.AxidxNms, dt/ms)))
        # print("B0",self.A.toarray()[:self.ny, N*self.ny:N*self.ny+self.nu])
    
    def dynamics(self, yi, u, dt, g, m, ms, s0):
        # Not needed for optimization, just to check
        y = np.copy(yi)
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
    
    def controlTest(self, dt, Qfdiag, m, ms, umin, umax, Nsim):
        # Hovering test
        ydes = np.zeros(self.ny)
        ydes[5] = 1 # sz
        y0 = np.copy(ydes)
        y0[0] = -1 # lower x
        
        # For lateral
        ydes[3] = 1 # sx

        model = self.toOSQP()

        ys = np.zeros((Nsim, self.ny))
        xs = np.zeros((Nsim, self.nx))
        yy = np.copy(y0)
        
        # nominal s
        s0 = [0.1,0,1]
        snom = [s0 for i in range(self.N)]

        for k in range(Nsim):
            # print(snom)
            # Update controller: copy out of update() for C version
            self.update(dt, snom, yy, Qfdiag, ydes, g, m, ms, umin, umax)
            # l,u update if needed
            model.update(Px=self.P.data, Ax_idx=np.asarray(self.Axidx), Ax=self.A.data[self.Axidx], q=self.q, l=self.l, u=self.u)
            res = model.solve()
            # print(res.info.status)

            xs[k,:] = res.x
            ys[k,:] = self.dynamics(yy, xs[k,self.N*self.ny : self.N*self.ny + self.nu], dt, g, m, ms, yy[3:6])
            # # normalize s
            # ys[k,3:6] /= np.linalg.norm(ys[k,3:6])
            yy = np.copy(ys[k,:])

            # use previous solution
            snom = [xs[k,i*self.ny+3:i*self.ny+6] for i in range(self.N)]
            ydes[3:6] = snom[-1] # FIXME: how to set this?

        # utest = np.zeros(3)
        # xtest = np.hstack((self.dynamics(y0, utest, dt, g, m, ms, s0), utest))
        # utest2 = -100*np.ones(3)
        # xtest2 = np.hstack((self.dynamics(y0, utest2, dt, g, m, ms, s0), utest2))
        # obj = lambda x : 0.5 * x @ self.P @ x + self.q @ x
        # print(obj(xtest2), obj(xtest), xtest2 - xtest)

        # print(y0)
        # print(xs)
        print(ys)
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2)
        ax[0].plot(ys[:,:3])
        ax[1].plot(ys[:,3:6])
        plt.show()

if __name__ == "__main__":
    # # WLQP gen
    # prob = qpSetupDense(4,4)
    N = 3
    dt = 2.0
    g = 9.81e-3
    m = 100.0
    ms = 1.0
    umax = np.array([1000.0, 1000.0, 1000.0])
    umin = -umax
    snom = [[0.1, 0.1, 1], [0.2, 0.1, 1], [0.3, 0.1, 1]]# Does not affect controlTest (only for initializing matrices)
    y0 = [1, 0.2, 0.1, 0.1, 0.2, 0.9, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    Qfdiag = [1, 1, 1, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2, 1e-2, 1e-3, 1e-3, 1e-3]
    ydes = [2, 0.2, 0.1, 0.4, 0.1, 1, 0.1, -0.1, -0.2, -0.3, -0.4, -0.5]

    up = UprightMPC(N)
    up.update(dt, snom, y0, Qfdiag, ydes, g, m, ms, umin, umax)
    up.dynamicsTest(N, dt, g, m, ms, snom, y0)

    up.controlTest(dt, Qfdiag, m, ms, umin, umax, 10)
    
    # # codegen
    # try:
    #     prob.codegen('gen', project_type='', force_rewrite=True, parameters='matrices', FLOAT=True, LONG=False)
    # except:
    #     # No worries if python module failed to compile
    #     pass
