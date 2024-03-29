
"""This file creates the python controller, and can return py ver, C ver, and a reactive"""

ny = 6
nu = 3

# Basic constituents of dynamics A0, B0 (only sparsity matters)
def getA0(T0):
    A0 = np.zeros((6, 6))
    A0[:3,3:] = T0*np.eye(3)
    return A0

def getB0(s0, Btau):
    return np.block([
        [np.reshape(s0, (3,1)),np.zeros((3,2))],
        [np.zeros((3,1)), Btau]
    ])
    
c0 = lambda g : np.array([0,0,-g,0,0,0])

e3h = skew([0,0,1])

def initConstraint(N, nx, nc):
    # these will be updated
    T0 = 1
    dt = 1
    s0 = np.ones(3)
    Btau = np.ones((3,2))

    A = np.zeros((nc, nx))
    P = np.eye(nx) # Q, R stacked
    # Dense blocks for WLQP
    P[N*ny:N*ny+ny,N*ny:N*ny+ny] = np.ones((ny,ny))
    P[N*ny:N*ny+ny,N*(2*ny+nu):N*(2*ny+nu)+4] = np.ones((ny,4))
    P[N*(2*ny+nu):N*(2*ny+nu)+4,N*ny:N*ny+ny] = np.ones((4,ny))
    P[N*(2*ny+nu):N*(2*ny+nu)+4,N*(2*ny+nu):N*(2*ny+nu)+4] = np.ones((4,4))

    # cols of A are broken up like this (partition nx)
    n1 = N*ny
    n2 = 2*N*ny
    n3 = 2*N*ny + nu*N
    # rows of A broken up:
    nc1 = N*ny # after this; yddot dynamics
    nc2 = 2*N*ny # after this; thrust lims
    nc3 = 2*N*ny + N # after this; Delta-u lims
    
    for k in range(N):
        # ykp1 = yk + dt*dyk equations
        A[k*ny:(k+1)*ny, k*ny:(k+1)*ny] = -np.eye(ny)
        A[k*ny:(k+1)*ny, n1 + k*ny:n1 + (k+1)*ny] = dt * np.eye(ny)
        if k>0:
            A[k*ny:(k+1)*ny, (k-1)*ny:(k)*ny] = np.eye(ny)
        
        # dykp1 equation
        A[nc1 + k*ny:nc1 + (k+1)*ny, n1 + k*ny:n1 + (k+1)*ny] = -np.eye(ny)
        A[nc1 + k*ny:nc1 + (k+1)*ny, n2 + k*nu:n2 + (k+1)*nu] = getB0(s0, Btau)
        if k>0:
            A[nc1 + k*ny:nc1 + (k+1)*ny, n1 + (k-1)*ny:n1 + (k)*ny] = np.eye(ny)
        if k>1:
            A[nc1 + k*ny:nc1 + (k+1)*ny, (k-2)*ny:(k-1)*ny] = getA0(T0*dt)
        
        # thrust lim
        A[nc2+k, n2+3*k] = 1
    
    # Delta-u WLQP rate limit constraint
    A[nc3:nc3+4, n3:n3+4] = np.eye(4)
    
    return sp.csc_matrix(A), sp.csc_matrix(P), P # for debugging

def updateConstraint(N, A, dt, T0, s0s, Btaus, y0, dy0, g, Tmax, delUL, delUU):
    nc = A.shape[0]

    # Update vector
    l = np.zeros(nc)
    y1 = y0 + dt * dy0
    l[:ny] = -y1
    for k in range(N):
        if k == 0:
            l[ny*N+k*ny : ny*N+(k+1)*ny] = -dy0 - dt*getA0(T0) @ y0 - dt*c0(g)
        elif k == 1:
            l[ny*N+k*ny : ny*N+(k+1)*ny] = -dt*getA0(T0) @ y1 - dt*c0(g)
        else:
            l[ny*N+k*ny : ny*N+(k+1)*ny] = -dt*c0(g)
    # copy for dynamics
    u = np.copy(l)
    # thrust lims
    for k in range(N):
        l[2*N*ny+k] = -T0
        u[2*N*ny+k] = Tmax-T0
    # Delta-u input (rate) limits
    l[2*N*ny+N:2*N*ny+N+4] = delUL
    u[2*N*ny+N:2*N*ny+N+4] = delUU

    # Left 1/4
    AxidxT0dt = []
    n2 = 2*ny + 3 # nnz in each block col on the left
    for k in range(N-2):
        AxidxT0dt += [n2*k + i for i in [8,11,14]]

    # Middle 2/4
    n1 = (2*N-1)*ny + (N-2)*3 # All the nnz in the left third
    n2 = 3*ny # nnz in each of the first N-1 block cols in the middle third
    Axidxdt = []
    for k in range(N):
        if k < N-1:
            Axidxdt += [n1 + n2*k + i for i in [0,3,6,9,12,15]]
        else:
            Axidxdt += [n1 + n2*k + i for i in [0,2,4,6,8,10]]
    
    # Right 3/4
    n1 += 3*ny*(N-1) + 2*ny # all nnz in the left and middle third
    n2 = 10 # nnz in each B0 + 1 for thrust lim
    Axidxs0 = []
    AxidxBtau = []
    for k in range(N):
        Axidxs0 += [n1 + n2*k + i for i in range(3)]
        AxidxBtau += [n1 + n2*k + 4 + i for i in range(6)]
    # No need to update rightmost

    # Last check
    assert A.nnz == n1 + n2*N + 4

    # Update
    A.data[AxidxT0dt] = dt*T0
    A.data[Axidxdt] = dt
    A.data[Axidxs0] = dt*np.hstack((s0s))
    A.data[AxidxBtau] = dt*np.hstack([np.ravel(Btau,order='F') for Btau in Btaus])

    Axidx = np.hstack((AxidxT0dt, Axidxdt, Axidxs0, AxidxBtau))
    # print("nAdata =",len(Axidx))

    # print(A[:,2*N*ny:2*N*ny+6].toarray())
    return A, l, u, Axidx

def getUpperTriang(P1):
    n = P1.shape[0]
    P1data = np.zeros(n*(n+1)//2)
    kk = 0
    for j in range(n):
        for i in range(j+1):
            P1data[kk] = P1[i,j]
            kk += 1
    return P1data

def updateObjective(N, Qyr, Qyf, Qdyr, Qdyf, R, ydes, dydes, Qw, dwdu, w0t, M0t, Pdense):
    # In the last column, need to stack the columns of the first matrix with the upper triang part of the second matrix
    mat1 = -M0t.T @ Qw @ dwdu # dy1,delu block
    mat2 = dwdu.T @ Qw @ dwdu # delu,delu block
    lastcol = np.zeros(6*4 + 4*(4+1)//2)
    offs = 0
    for j in range(mat1.shape[1]):
        lastcol[offs : offs+6] = mat1[:,j]
        offs += 6
        lastcol[offs : offs+j+1] = mat2[:(j+1),j]
        offs += j+1

    # Block diag components - see notes
    Pdata = np.hstack((
        np.hstack([Qyr for k in range(N-1)]),
        Qyf,
        getUpperTriang(np.diag(Qdyr) + M0t.T @ Qw @ M0t),# dy1,dy1 block upper triang
        np.hstack([Qdyr for k in range(N-2)]),
        Qdyf,
        np.hstack([R for k in range(N)]),
        lastcol
    ))
    # Dense P update for debugging
    ii = np.diag_indices(N*(2*ny+nu))
    Pdense[:N*(2*ny+nu), :N*(2*ny+nu)][ii] = np.hstack((
        np.hstack([Qyr for k in range(N-1)]),
        Qyf,
        np.hstack([Qdyr for k in range(N-1)]),
        Qdyf,
        np.hstack([R for k in range(N)])
    ))
    # Off diag parts
    Pdense[N*ny:(N+1)*ny, N*ny:(N+1)*ny] += M0t.T @ Qw @ M0t # keep the Qyr
    Pdense[N*(2*ny+nu):, N*(2*ny+nu):] = dwdu.T @ Qw @ dwdu
    Pdense[N*ny:(N+1)*ny, N*(2*ny+nu):] = -M0t.T @ Qw @ dwdu
    Pdense[N*(2*ny+nu):, N*ny:(N+1)*ny] = -dwdu.T @ Qw @ M0t

    q = np.hstack((
        np.hstack([-Qyr*ydes for k in range(N-1)]),
        -Qyf*ydes,
        -Qdyr*dydes - M0t.T @ Qw @ w0t, # dy1
        np.hstack([-Qdyr*dydes for k in range(N-2)]),
        -Qdyf*dydes,
        np.zeros(N*len(R)),
        dwdu.T @ Qw @ w0t
    ))
    return Pdata, q, Pdense

def openLoopX(N, dt, T0, s0s, Btaus, y0, dy0, g):
    ys = np.zeros((N,ny))
    dys = np.zeros((N,ny))
    us = np.random.rand(N,nu)

    yk = np.copy(y0)
    dyk = np.copy(dy0)
    for k in range(N):
        # at k=0, yy=y0, 
        dykp1 = dyk + dt*(getA0(T0) @ yk + getB0(s0s[k], Btaus[k]) @ us[k,:] + c0(g)) # dy[k+1]

        dys[k,:] = dykp1
        ykp1 = yk + dt * dyk # y[k+1]
        ys[k,:] = ykp1 + dt * dykp1 # y[k+2]

        # For next k
        yk = np.copy(ykp1)
        dyk = np.copy(dykp1)

    # stack
    x = np.hstack((np.ravel(ys), np.ravel(dys), np.ravel(us), np.zeros(4)))
    # print(ys, x)
    return x

class UprightMPC2():
    def __init__(self, N, dt, g, TtoWmax, ws, wds, wpr, wpf, wvr, wvf, wthrust, wmom, umin, umax, dumax, mb, Ib, Qw, controlRate):
        self.N = N

        nx = self.N * (2*ny + nu) + 4
        nc = 2*self.N*ny + self.N + 4

        self.A, self.P, self.Pdense = initConstraint(N, nx, nc)

        Qyr = np.hstack((np.full(3,wpr), np.full(3,ws)))
        Qyf = np.hstack((np.full(3,wpf), np.full(3,ws)))
        Qdyr = np.hstack((np.full(3,wvr), np.full(3,wds)))
        Qdyf = np.hstack((np.full(3,wvf), np.full(3,wds)))
        R = np.hstack((wthrust,np.full(2,wmom)))

        self.dt = dt
        self.Wts = [Qyr, Qyf, Qdyr, Qdyf, R]
        self.g = g
        self.Tmax = TtoWmax * g # use thrust-to-weight ratio to set max specific thrust

        # Create OSQP
        self.model = osqp.OSQP()
        self.model.setup(P=self.P, A=self.A, l=np.zeros(nc), eps_rel=1e-7, eps_abs=1e-7, verbose=False)

        # Manage linearization point
        self.T0 = 0 # mass-specific thrust
        self.Ibi = np.diag(1/Ib)
        self.M0 = np.hstack((mb,mb,mb,Ib))

        # WLQP state
        self.umin = umin
        self.umax = umax
        self.dumax = dumax / controlRate
        self.Qw = np.diag(Qw)
        # what "u" is depends on w(u). Here in python testing with w(u) = [0,0,u0,u1,u2,u3]
        self.u0 = np.zeros(4)
        
    def codegen(self, dirname='uprightmpc2/gen'):
        try:
            self.model.codegen(dirname, project_type='', force_rewrite=True, parameters='matrices', FLOAT=True, LONG=False)
        except:
            # No worries if python module failed to compile
            pass

    def testDyn(self, T0sp, s0s, Btaus, y0, dy0):
        # Test
        self.A, l, u, Axidx = updateConstraint(self.N, self.A, self.dt, T0sp, s0s, Btaus, y0, dy0, self.g, self.Tmax, np.zeros(4), np.zeros(4))
        xtest = openLoopX(self.N, self.dt, T0sp, s0s, Btaus, y0, dy0, self.g)
        print((self.A @ xtest - l)[:2*self.N*ny])
    
    def update1(self, T0sp, s0s, Btaus, y0, dy0, ydes, dydes, dwdu0, w0t, M0t):
        # WLQP Delta-u limit
        delUL = np.copy(-self.dumax)
        delUU = np.copy(self.dumax)
        # Input limits
        for i in range(4):
            if self.u0[i] < self.umin[i]:
                delUL[i] = 0
            elif self.u0[i] > self.umax[i]:
                delUU[i] = 0
        # Update
        self.A, self.l, self.u, self.Axidx = updateConstraint(self.N, self.A, self.dt, T0sp, s0s, Btaus, y0, dy0, self.g, self.Tmax, delUL, delUU)
        self.Pdata, self.q, self.Pdense = updateObjective(self.N, *self.Wts, ydes, dydes, self.Qw, dwdu0, w0t, M0t, self.Pdense)
        
        # OSQP solve ---
        self.model.update(Px=self.Pdata, Ax=self.A.data, q=self.q, l=self.l, u=self.u)
        res = self.model.solve()
        if 'solved' not in res.info.status:
            print(res.info.status)
        self.obj_val = res.info.obj_val
        # Functions for debugging
        self.obj = lambda x : 0.5 * x.T @ self.Pdense @ x + self.q.T @ x
        self.viol = lambda x : np.amin(np.hstack((self.A @ x - self.l, self.u - self.A @ x)))
        return res.x
    
    def update2(self, p0, R0, dq0, pdes, dpdes, sdes):
        # what "u" is depends on w(u). Here in python testing with w(u) = [0,0,u0,u1,u2,u3]
        w0 = np.hstack((0,0,0,0,0,0))#self.u0))
        dwdu0 = np.vstack((np.zeros((2,4)), np.eye(4)))
        # At current state
        s0 = np.copy(R0[:,2])
        s0s = [s0 for i in range(self.N)]
        Btau = (-R0 @ e3h @ self.Ibi)[:,:2] # no yaw torque
        Btaus = [Btau for i in range(self.N)]
        ds0 = -R0 @ e3h @ dq0[3:6] # omegaB

        y0 = np.hstack((p0, s0))
        dy0 = np.hstack((dq0[:3], ds0))
        ydes = np.hstack((pdes, sdes))
        dydes = np.hstack((dpdes, 0, 0, 0))

        h0 = np.hstack((R0.T @ np.array([0, 0, self.M0[0] * self.g]), np.zeros(3)))
        T0 = np.eye(6)
        T0[3:,3:] = e3h @ R0.T
        # M0t = M0*T0/dt
        M0t = np.diag(self.M0) @ T0 / self.dt
        # w0t = w0 - h0 + M0*dq0/dt
        w0t = w0 - h0 + (self.M0 * dq0) / self.dt

        self.prevsol = self.update1(self.T0, s0s, Btaus, y0, dy0, ydes, dydes, dwdu0, w0t, M0t)
        utilde = self.prevsol[2*ny*self.N : 2*ny*self.N+nu]
        self.T0 += utilde[0]

        # WLQP update u0
        delu = self.prevsol[(2*ny + nu)*self.N:]
        self.u0 += delu

        return np.hstack((self.T0, utilde[1:]))
    
    def getAccDes(self, R0, dq0):
        dy1des = self.prevsol[ny*self.N : ny*self.N+ny] # from horiz
        # # Coordinate change for the velocity
        # bTw = lambda dq : np.hstack((R0.T @ dq[:3], dq[3:6]))
        dq1des = np.hstack((dy1des[:3], e3h @ R0.T @ dy1des[3:6])) # NOTE omegaz is lost
        # return (bTw(dq1des) - bTw(dq0)) / self.dt
        return (dq1des - dq0) / self.dt # return in world frame
    
    def update(self, p0, R0, dq0, pdes, dpdes, sdes):
        # Version of above that computes the desired body frame acceleration
        u = self.update2(p0, R0, dq0, pdes, dpdes, sdes)
        return u, self.getAccDes(R0, dq0), self.u0