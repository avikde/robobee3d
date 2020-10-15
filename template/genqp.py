import osqp
import numpy as np
import scipy.sparse as sp
from scipy.spatial.transform import Rotation
from scipy.linalg import expm
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

skew = lambda v : np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    
def quadrotorNLVF(p, Rb, dq, u):
    mb = 100
    ib = 1000
    omega = dq[3:6]
    
    dv = u[0] * Rb @ np.array([0,0,1]) / mb
    domega = (-np.cross(omega, ib * omega) + u[1:4]) / ib

    return np.hstack((dv, domega))

def quadrotorNLDyn(p, Rb, dq, u, dt):
    ddq = quadrotorNLVF(p, Rb, dq, u)
    # Euler integrate
    p = p + dt * dq[0:3]
    Rb = Rb @ expm(skew(dq[3:6]) * dt)
    dq = dq + dt * ddq
    return p, Rb, dq

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

        # (I + dt*vT0*N)
        A0 = np.eye(self.nq) + 1 * np.diag(np.ones(3), k=3) # using 1 for dt*vT0

        # Construct dynamics constraint
        ncon = N * (self.nq + 3)
        self.A = np.zeros((ncon, self.nx))
        self.l = np.zeros(ncon)
        self.u = np.zeros(ncon)
        for k in range(N):
            # q(k+1) = Ad(vT0)*q(k) + Bd(sk)*vk
            self.A[k*self.nq:(k+1)*self.nq, k*self.nq:(k+1)*self.nq] = -np.eye(self.nq) # for -q1...qN on the LHS
            if k > 0:
                self.A[k*self.nq:(k+1)*self.nq, (k-1)*self.nq:(k)*self.nq] = A0 # for q
            self.A[k*self.nq:(k+1)*self.nq, (N*self.nq + k*self.nu):(N*self.nq + (k+1)*self.nu)] = Bs(snom[k]) # Bd(sk)

            # only in the first eqn
            if k == 0:
                self.l[k*self.nq:(k+1)*self.nq] += -A0 @ q0
                self.u[k*self.nq:(k+1)*self.nq] += -A0 @ q0

            # s limits
            self.A[(N*self.nq+3*k):(N*self.nq+(k+1)*3), (k)*self.nq+3:(k)*self.nq+6] = np.eye(3)
            self.l[(N*self.nq+3*k):(N*self.nq+(k+1)*3)] = smin
            self.u[(N*self.nq+3*k):(N*self.nq+(k+1)*3)] = smax
        self.A = sp.csc_matrix(self.A)
        # print(A, c)

        # construct objective. 
        self.q = np.zeros(self.nx)
        self.P = np.zeros((self.nx, self.nx))
        # final cost
        np.fill_diagonal(self.P[(N-1)*self.nq : N*self.nq, (N-1)*self.nq : N*self.nq], Qfdiag)
        for k in range(N):
            np.fill_diagonal(self.P[N*self.nq + k*self.nu : N*self.nq + (k+1)*self.nu, N*self.nq + k*self.nu : N*self.nq + (k+1)*self.nu], Rdiag)
        self.P = sp.csc_matrix(self.P)
        self.q[(N-1)*self.nq:N*self.nq] = -(np.asarray(Qfdiag) * np.asarray(qdes))

        self.saveAxidx()
        self.toOSQP()
        self.resetNominal()

    def resetNominal(self):
        # Nominal traj management
        self.snom = [[0,0,1] for i in range(self.N)]
        self.vT0 = 0 # TODO: what to init at?
    
    def saveAxidx(self):
        A1nnz = 18*self.N-9
        Bnnz = 7
        A2nnz = self.N * Bnnz

        # To test dynamics constraint after need to update sparse csc_matrix
        assert self.A.nnz == A1nnz + A2nnz

        # put together - store these. 
        idxdtvT0 = [[18*k + i for i in [7,11,15]] for k in range(self.N-1)]
        idxdtvT0 = sum(idxdtvT0, []) # join the list of lists
        # Should be filled with stacked Bi = dt*(sx,sy,sz,1,-sx/sz,1,-sy/sz)
        idxBi = range(A1nnz, A1nnz + A2nnz)

        self.Axidx = list(idxdtvT0) + list(idxBi)
    
    def update(self, q0, qdes, Qfdiag, Rdiag, smin, smax, dt, snom, vT0):
        # should not need the arguments in constructor (just sets sparsity)
        A0 = np.eye(self.nq) + dt * vT0 * np.diag(np.ones(3), k=3)

        # update l, u
        self.l[:self.nq] = self.u[:self.nq] = -A0 @ np.asarray(q0)
        self.l[self.N*self.nq:] = np.tile(smin, self.N)
        self.u[self.N*self.nq:] = np.tile(smax, self.N)

        # update q
        self.q[(self.N-1)*self.nq:self.N*self.nq] = -(np.asarray(Qfdiag) * np.asarray(qdes))

        # update P
        self.P.data = np.hstack((Qfdiag, np.tile(Rdiag, self.N))) # replace the whole thing

        # update A
        # print("hi",np.hstack(snom))
        dtvT0data = np.full(3*(self.N-1), dt*vT0)
        Bidata = np.hstack([
            dt * np.array([snom[k][0],snom[k][1],snom[k][2],1,-snom[k][0]/snom[k][2],1,-snom[k][1]/snom[k][2]])
            for k in range(self.N)])
        self.A.data[self.Axidx] = np.hstack((dtvT0data, Bidata))

        # print(self.A[:N*self.nq,:N*self.nq].toarray())
        # print("B0",self.A.toarray()[:self.ny, N*self.ny:N*self.ny+self.nu])

        # OSQP solve ---
        self.model.update(Px=self.P.data, Ax_idx=np.asarray(self.Axidx), Ax=self.A.data[self.Axidx], q=self.q, l=self.l, u=self.u)
        res = self.model.solve()
        if 'solved' not in res.info.status:
            print(res.info.status)
        uu = res.x[self.N*self.nq : self.N*self.nq + self.nu]
        
        # Nominal traj management
        self.snom = [res.x[i*self.nq+3:i*self.nq+6] for i in range(self.N)]
        self.vT0 += uu[0]
        return res.x, uu
    
    def dynamics(self, qi, u, dt, s0, vT0):
        # Not needed for optimization, just to check
        q = np.copy(qi)
        s = q[3:]
        dvT = u[0] # delta thrust
        vM = u[1:] # moment
        dq = np.hstack((vT0 * s + dvT * np.asarray(s0), vM, np.dot(-np.asarray(s0[:2])/s0[2], vM)))
        return q + dt * dq

    def dynamicsNLVF(self, q, u):
        """For simulation nonlinear vector field. u is assumed to be full vT, vM"""
        s = q[3:]
        vT = u[0] # actual thrust
        vM = u[1:] # moment
        return np.hstack((vT * s, vM, np.dot(-np.asarray(s[:2])/s[2], vM)))

    def dynamicsTest(self, dt, snom, q0, vT0):
        qs = np.zeros((self.N, self.nq))
        us = np.random.rand(self.N, self.nu)
        qq = np.copy(q0)

        for k in range(N):
            qs[k,:] = self.dynamics(qq, us[k,:], dt, snom[k], vT0)
            qq = qs[k,:]
        
        # reshape into an "x"
        xtest = np.hstack((np.ravel(qs), np.ravel(us)))
        # print(xtest.shape, up.nx)
        e1 = self.A @ xtest - self.l
        e2 = self.u - self.A @ xtest
        Z = np.zeros(N*self.nq)
        if np.allclose(e1[:N*self.nq], Z) and np.allclose(e2[:N*self.nq], Z):
            print('Dynamics test passed')
        else:
            print('Dynamics test FAILED', e1, e2)
        
    def toOSQP(self):
        # osqp
        self.model = osqp.OSQP()
        self.model.setup(P=self.P, q=self.q, A=self.A, l=self.l, u=self.u, eps_rel=1e-4, eps_abs=1e-4, verbose=False)
    
    def controlTest(self, dt, Qfdiag, Rdiag, smin, smax, tend, dtsim=0.5, simmodel=0):
        """simmodel=1 means use nonlinear VF with small dtsim.
        simmodel=2 means second order quadrotor"""
        # Hovering test
        qdes = np.zeros(self.nq)
        qdes[5] = 1 # sz

        # Initial condition for anchor
        if simmodel == 2:
            y0 = np.zeros(9) # p, dq
            y0[:3] = np.array([-1, 0.5, -1])
            Rb0 = np.eye(3)
            RRb = np.copy(Rb0)
        else:
            y0 = np.copy(qdes)
            y0[:3] = np.array([-1, 0.5, -1])
        yy = np.copy(y0)
        
        # For lateral
        # qdes[3] = 1 # sx

        Nsim = int(tend//dtsim if simmodel > 0 else tend//dt)
        tt = np.linspace(0, tend, num=Nsim)
        ys = np.zeros((Nsim, len(y0)))
        xTs = np.zeros((Nsim, self.nx))
        uTs = np.zeros((Nsim, 3))

        self.resetNominal()

        for k in range(Nsim):
            # "Template projection"
            if simmodel == 2:
                qT = np.hstack((yy[:3], RRb @ np.array([0,0,1]))) # p,s
            else:
                qT = yy
            # print(snom)

            # Update controller: copy out of update() for C version
            xTs[k,:], uu = self.update(qT, qdes, Qfdiag, Rdiag, smin, smax, dt, self.snom, self.vT0)
            uTs[k,:] = np.hstack((self.vT0, uu[1:]))

            # Convert back to anchor
            if simmodel == 2:
                uA = np.hstack((uTs[k,:], 0)) # yaw moment
            else:
                uA = uTs[k,:]

            # Integrate dynamics
            if simmodel == 2:
                p, RRb, dq = quadrotorNLDyn(yy[:3], RRb, yy[3:], uA, dtsim)
                ys[k,:] = np.hstack((p, dq))
            elif simmodel == 1:
                dydt = self.dynamicsNLVF(yy, uA) # no local lin stuff
                ys[k,:] = yy + dtsim * dydt
            else:
                ys[k,:] = self.dynamics(yy, uA, dt, yy[3:6], self.vT0)

            if simmodel < 2:
                # normalize s
                ys[k,3:6] /= np.linalg.norm(ys[k,3:6])
            yy = np.copy(ys[k,:])

        # utest = np.zeros(3)
        # xtest = np.hstack((self.dynamics(y0, utest, dt, g, m, ms, s0), utest))
        # utest2 = -100*np.ones(3)
        # xtest2 = np.hstack((self.dynamics(y0, utest2, dt, g, m, ms, s0), utest2))
        # obj = lambda x : 0.5 * x @ self.P @ x + self.q @ x
        # print(obj(xtest2), obj(xtest), xtest2 - xtest)

        # print(y0)
        # print(qs)
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3)
        ax[0].plot(tt, ys[:,:3])
        ax[0].set_ylabel('q')
        ax[1].plot(tt, ys[:,3:6])
        ax[1].set_ylabel('s')
        ax[2].plot(tt, uTs)
        ax[2].set_ylabel('u')
        plt.show()

if __name__ == "__main__":
    # # WLQP gen
    # prob = qpSetupDense(4,4)
    N = 3
    dt = 3
    g = 9.81e-3
    m = 100.0
    ms = 1.0
    smax = np.array([0.5, 0.5, 1.5])
    smin = -smax
    snom = [np.ones(3) for i in range(N)]
    q0 = [1, 0.2, 0.1, 0.1, 0.2, 0.9]
    qdes = [2, 0.2, 0.1, 0.4, 0.1, 1]
    Qfdiag = [100, 100, 10, 100,100,100]
    Rdiag = [100, 100, 100]
    vT0 = 1

    up = UprightMPC(N)
    up.update(q0, qdes, Qfdiag, Rdiag, smin, smax, dt, snom, vT0)
    up.dynamicsTest(dt, snom, q0, vT0)

    up.controlTest(dt, Qfdiag, Rdiag, smin, smax, 200, simmodel=1)
    
    # # codegen
    # try:
    #     prob.codegen('gen', project_type='', force_rewrite=True, parameters='matrices', FLOAT=True, LONG=False)
    # except:
    #     # No worries if python module failed to compile
    #     pass
