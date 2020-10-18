import osqp
import numpy as np
import scipy.sparse as sp
np.set_printoptions(precision=4, suppress=True, linewidth=200)
from genqp import quadrotorNLDyn

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

def initConstraint(N, nx, nc):
    # these will be updated
    T0 = 1
    dt = 1
    s0 = np.ones(3)
    Btau = np.ones((3,2))

    A = np.zeros((nc, nx))
    P = np.eye(nx) # Q, R stacked

    # cols of A are broken up like this (partition nx)
    n1 = N*ny
    n2 = 2*N*ny
    # rows of A broken up:
    nc1 = N*ny
    
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
    
    return sp.csc_matrix(A), sp.csc_matrix(P)

def updateConstraint(N, A, dt, T0, s0s, Btaus, y0, dy0, g):
    nc = A.shape[0]

    # Update vector
    lu = np.zeros(nc)
    y1 = y0 + dt * dy0
    lu[:ny] = -y1
    for k in range(N):
        if k == 0:
            lu[ny*N+k*ny : ny*N+(k+1)*ny] = -dy0 - dt*getA0(T0) @ y0 - dt*c0(g)
        elif k == 1:
            lu[ny*N+k*ny : ny*N+(k+1)*ny] = -dt*getA0(T0) @ y1 - dt*c0(g)
        else:
            lu[ny*N+k*ny : ny*N+(k+1)*ny] = -dt*c0(g)

    # Left third
    AxidxT0dt = []
    n2 = 2*ny + 3 # nnz in each block col on the left
    for k in range(N-2):
        AxidxT0dt += [n2*k + i for i in [8,11,14]]

    # Middle third
    n1 = (2*N-1)*ny + (N-2)*3 # All the nnz in the left third
    n2 = 3*ny # nnz in each of the first N-1 block cols in the middle third
    Axidxdt = []
    for k in range(N):
        if k < N-1:
            Axidxdt += [n1 + n2*k + i for i in [0,3,6,9,12,15]]
        else:
            Axidxdt += [n1 + n2*k + i for i in [0,2,4,6,8,10]]
    
    # Right third
    n1 += 3*ny*(N-1) + 2*ny # all nnz in the left and middle third
    n2 = 9 # nnz in each B0
    Axidxs0 = []
    AxidxBtau = []
    for k in range(N):
        Axidxs0 += [n1 + n2*k + i for i in range(3)]
        AxidxBtau += [n1 + n2*k + 3 + i for i in range(6)]

    # Last check
    assert A.nnz == n1 + n2*N

    # Update
    A.data[AxidxT0dt] = dt*T0
    A.data[Axidxdt] = dt
    A.data[Axidxs0] = dt*np.hstack((s0s))
    A.data[AxidxBtau] = dt*np.hstack([np.ravel(Btau,order='F') for Btau in Btaus])

    # print(A[:,36:42].toarray())
    return A, lu

def updateObjective(N, P, Qyr, Qyf, Qdyr, Qdyf, R, ydes, dydes):
    P.data = np.hstack((
        np.hstack([Qyr for k in range(N-1)]),
        Qyf,
        np.hstack([Qdyr for k in range(N-1)]),
        Qdyf,
        np.hstack([R for k in range(N)])
    ))
    q = np.hstack((
        np.hstack([-Qyr*ydes for k in range(N-1)]),
        -Qyf*ydes,
        np.hstack([-Qdyr*dydes for k in range(N-1)]),
        -Qdyf*dydes,
        np.zeros(N*len(R))
    ))
    return P, q

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
    x = np.hstack((np.ravel(ys), np.ravel(dys), np.ravel(us)))
    # print(ys, x)
    return x

class UprightMPC2():
    def __init__(self, N, dt, Qyr, Qyf, Qdyr, Qdyf, R, g):
        self.N = N

        nx = self.N * (2*ny + nu)
        nc = 2*self.N*ny

        self.A, self.P = initConstraint(N, nx, nc)

        self.dt = dt
        self.Wts = [Qyr, Qyf, Qdyr, Qdyf, R]
        self.g = g

        # Create OSQP
        self.model = osqp.OSQP()
        self.model.setup(P=self.P, A=self.A, l=np.zeros(nc), eps_rel=1e-4, eps_abs=1e-4, verbose=False)

    def testDyn(self, T0sp, s0s, Btaus, y0, dy0):
        # Test
        self.A, lu = updateConstraint(self.N, self.A, self.dt, T0sp, s0s, Btaus, y0, dy0, self.g)
        xtest = openLoopX(self.N, self.dt, T0, s0s, Btaus, y0, dy0, self.g)
        print(self.A @ xtest - lu)
    
    def update(self, T0sp, s0s, Btaus, y0, dy0, ydes, dydes):
        self.A, lu = updateConstraint(self.N, self.A, self.dt, T0sp, s0s, Btaus, y0, dy0, self.g)
        self.P, q = updateObjective(self.N, self.P, *self.Wts, ydes, dydes)
        
        # OSQP solve ---
        self.model.update(Px=self.P.data, Ax=self.A.data, q=q, l=lu, u=lu)
        res = self.model.solve()
        if 'solved' not in res.info.status:
            print(res.info.status)
        return res.x

def controlTest(mdl, tend, dtsim=0.5):
    # Initial conditions
    p = np.array([0, 0, -1])
    Rb = np.eye(3)
    dq = np.zeros(6)
    
    tt = np.arange(tend, step=dtsim)
    Nt = len(tt)

    # for logging
    ys = np.zeros((Nt, 12))
    us = np.zeros((Nt, 3))

    for ti in range(Nt):
        u = np.array([1,0.1,0])
        p, Rb, dq = quadrotorNLDyn(p, Rb, dq, u, dtsim)
        ys[ti,:] = np.hstack((p, Rb[:,2], dq))
        us[ti,:] = u
    
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(3,2)
    ax = ax.ravel()
        
    ax[0].plot(tt, ys[:,:3])
    ax[0].set_ylabel('p')
    ax[1].plot(tt, ys[:,3:6])
    ax[1].set_ylabel('s')
    ax[2].plot(tt, us)
    ax[2].set_ylabel('u')
    ax[3].plot(tt, ys[:,6:9])
    ax[3].set_ylabel('v')
    ax[4].plot(tt, ys[:,9:12])
    ax[4].set_ylabel('omega')
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    T0 = 0.5
    dt = 0.5
    N = 3
    s0s = [[0.1,0.1,0.9] for i in range(N)]
    Btaus = [np.full((3,2),1.123) for i in range(N)]
    y0 = np.random.rand(6)
    dy0 = np.random.rand(6)
    g = 9.81e-3
    Qyr = np.array([1,1,1,1,1,1])
    Qyf = np.array([1,1,1,1,1,1])
    Qdyr = np.array([1,1,1,1,1,1])
    Qdyf = np.array([1,1,1,1,1,1])
    R = np.array([1,1,1])
    ydes = np.zeros_like(y0)
    dydes = np.zeros_like(y0)

    up = UprightMPC2(N, dt, Qyr, Qyf, Qdyr, Qdyf, R, g)
    up.testDyn(T0, s0s, Btaus, y0, dy0)

    controlTest(up, 50)

