import osqp
import numpy as np
import scipy.sparse as sp
from scipy.spatial.transform import Rotation
np.set_printoptions(precision=4, suppress=True, linewidth=200)
from genqp import quadrotorNLDyn, skew, Ib

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

    # cols of A are broken up like this (partition nx)
    n1 = N*ny
    n2 = 2*N*ny
    # rows of A broken up:
    nc1 = N*ny # after this; yddot dynamics
    nc2 = 2*N*ny # after this; s lims
    nc3 = 2*N*ny + 3*N # after this; thrust lims
    
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
        
        # s lim
        A[nc2+3*k : nc2+3*(k+1), k*ny+3:(k+1)*ny] = np.eye(3)
        # thrust lim
        A[nc3+k, n2+3*k] = 1
    
    return sp.csc_matrix(A), sp.csc_matrix(P)

def updateConstraint(N, A, dt, T0, s0s, Btaus, y0, dy0, g, smin, smax, Tmax):
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
    # s lims
    for k in range(N):
        l[2*N*ny+3*k:2*N*ny+3*(k+1)] = smin
        u[2*N*ny+3*k:2*N*ny+3*(k+1)] = smax
    # thrust lims
    for k in range(N):
        l[2*N*ny+3*N+k] = -T0
        u[2*N*ny+3*N+k] = Tmax-T0

    # Left third
    AxidxT0dt = []
    n2 = 2*ny + 6 # nnz in each block col on the left
    for k in range(N-2):
        AxidxT0dt += [n2*k + i for i in [8,12,16]]

    # Middle third
    n1 = (2*N-1)*ny + (N-2)*3 + 3*N # All the nnz in the left third
    n2 = 3*ny # nnz in each of the first N-1 block cols in the middle third
    Axidxdt = []
    for k in range(N):
        if k < N-1:
            Axidxdt += [n1 + n2*k + i for i in [0,3,6,9,12,15]]
        else:
            Axidxdt += [n1 + n2*k + i for i in [0,2,4,6,8,10]]
    
    # Right third
    n1 += 3*ny*(N-1) + 2*ny # all nnz in the left and middle third
    n2 = 10 # nnz in each B0 + 1 for thrust lim
    Axidxs0 = []
    AxidxBtau = []
    for k in range(N):
        Axidxs0 += [n1 + n2*k + i for i in range(3)]
        AxidxBtau += [n1 + n2*k + 4 + i for i in range(6)]

    # Last check
    assert A.nnz == n1 + n2*N

    # Update
    A.data[AxidxT0dt] = dt*T0
    A.data[Axidxdt] = dt
    A.data[Axidxs0] = dt*np.hstack((s0s))
    A.data[AxidxBtau] = dt*np.hstack([np.ravel(Btau,order='F') for Btau in Btaus])

    # print(A[:,2*N*ny:2*N*ny+6].toarray())
    return A, l, u

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
    def __init__(self, N, dt, g, smin, smax, TtoWmax, ws, wds, wpr, wpf, wvr, wvf, wthrust, wmom):
        self.N = N

        nx = self.N * (2*ny + nu)
        nc = 2*self.N*ny + 4*self.N

        self.A, self.P = initConstraint(N, nx, nc)

        Qyr = np.hstack((np.full(3,wpr), np.full(3,ws)))
        Qyf = np.hstack((np.full(3,wpf), np.full(3,ws)))
        Qdyr = np.hstack((np.full(3,wvr), np.full(3,wds)))
        Qdyf = np.hstack((np.full(3,wvf), np.full(3,wds)))
        R = np.hstack((wthrust,np.full(2,wmom)))

        self.dt = dt
        self.Wts = [Qyr, Qyf, Qdyr, Qdyf, R]
        self.g = g
        self.smin = smin
        self.smax = smax
        self.Tmax = TtoWmax * g # use thrust-to-weight ratio to set max specific thrust

        # Create OSQP
        self.model = osqp.OSQP()
        self.model.setup(P=self.P, A=self.A, l=np.zeros(nc), eps_rel=1e-4, eps_abs=1e-4, verbose=False)

        # Manage linearization point
        self.T0 = 0 # mass-specific thrust
        self.Ibi = np.linalg.inv(Ib)

    def testDyn(self, T0sp, s0s, Btaus, y0, dy0):
        # Test
        self.A, l, u = updateConstraint(self.N, self.A, self.dt, T0sp, s0s, Btaus, y0, dy0, self.g, self.smin, self.smax, self.Tmax)
        xtest = openLoopX(self.N, self.dt, T0, s0s, Btaus, y0, dy0, self.g)
        print((self.A @ xtest - l)[:2*self.N*ny])
    
    def update(self, T0sp, s0s, Btaus, y0, dy0, ydes, dydes):
        self.A, l, u = updateConstraint(self.N, self.A, self.dt, T0sp, s0s, Btaus, y0, dy0, self.g, self.smin, self.smax, self.Tmax)
        self.P, q = updateObjective(self.N, self.P, *self.Wts, ydes, dydes)
        
        # OSQP solve ---
        self.model.update(Px=self.P.data, Ax=self.A.data, q=q, l=l, u=u)
        res = self.model.solve()
        if 'solved' not in res.info.status:
            print(res.info.status)

        return res.x
    
    def update2(self, p0, R0, dq0, pdes, dpdes):
        # At current state
        s0 = np.copy(R0[:,2])
        s0s = [s0 for i in range(self.N)]
        Btau = (-R0 @ e3h @ self.Ibi)[:,:2] # no yaw torque
        Btaus = [Btau for i in range(self.N)]
        ds0 = -R0 @ e3h @ dq0[3:6] # omegaB

        y0 = np.hstack((p0, s0))
        dy0 = np.hstack((dq0[:3], ds0))
        ydes = np.hstack((pdes, 0, 0, 1))
        dydes = np.hstack((dpdes, 0, 0, 0))

        self.prevsol = self.update(self.T0, s0s, Btaus, y0, dy0, ydes, dydes)
        utilde = self.prevsol[2*ny*self.N : 2*ny*self.N+nu]
        self.T0 += utilde[0]

        return np.hstack((self.T0, utilde[1:]))
    
    def getAccDes(self, R0, dq0):
        dy1des = self.prevsol[ny*self.N : ny*self.N+ny] # from horiz
        # # Coordinate change for the velocity
        # bTw = lambda dq : np.hstack((R0.T @ dq[:3], dq[3:6]))
        dq1des = np.hstack((dy1des[:3], e3h @ R0.T @ dy1des[3:6])) # NOTE omegaz is lost
        # return (bTw(dq1des) - bTw(dq0)) / self.dt
        return (dq1des - dq0) / self.dt # return in world frame
    
    def updateGetAccdes(self, p0, R0, dq0, pdes, dpdes):
        # Version of above that computes the desired body frame acceleration
        self.update2(p0, R0, dq0, pdes, dpdes)
        return self.getAccDes(R0, dq0)

def reactiveController(p, Rb, dq, pdes, kpos=[1e-3,5e-1], kz=[1e-1,1e0], ks=[1e0,1e2]):
    # Pakpong-style reactive controller
    sdes = np.clip(kpos[0] * (pdes - p) - kpos[1] * dq[:3], np.full(3, -0.5), np.full(3, 0.5))
    sdes[2] = 1
    # sdes = np.array([0,0,1])
    omega = dq[3:]
    dp = dq[:3]
    s = Rb[:,2]
    ds = -Rb @ e3h @ omega
    # Template controller <- LATEST
    fz = kz[0] * (pdes[2] - p[2]) - kz[1] * dq[2]
    fTorn = ks[0] * (s - sdes) + ks[1] * ds
    fTorn[2] = 0
    fAorn = -e3h @ Rb.T @ fTorn
    return np.hstack((fz, fAorn[:2]))

def controlTest(mdl, tend, dtsim=0.2, useMPC=True, trajFreq=0, trajAmp=0, ascentIC=False):
    """trajFreq in Hz, trajAmp in mm"""
    # Initial conditions
    dq = np.zeros(6)
    if ascentIC:
        p = np.array([0, 0, -50])
        Rb = np.eye(3)
    else:
        p = np.array([0, 0, -1])
        Rb = Rotation.from_euler('xyz', np.ones(3)).as_matrix()
        dq[0] = 0.1
    pdes = np.zeros(3)
    dpdes = np.zeros(3)
    
    tt = np.arange(tend, step=dtsim)
    Nt = len(tt)

    # for logging
    ys = np.zeros((Nt, 12))
    us = np.zeros((Nt, 3))
    pdess = np.zeros((Nt, 3))
    accdess = np.zeros((Nt,6))

    trajOmg = 2 * np.pi * trajFreq * 1e-3 # to KHz, then to rad/ms
    ddqdes = None # test integrate ddq sim below

    for ti in range(Nt):
        # Traj to follow
        pdes[0] = trajAmp * np.sin(trajOmg * tt[ti])
        dpdes[0] = trajAmp * trajOmg * np.cos(trajOmg * tt[ti])

        # Call controller
        if useMPC:
            u = mdl.update2(p, Rb, dq, pdes, dpdes)
            accdess[ti,:] = mdl.getAccDes(Rb, dq)
            # # Alternate simulation by integrating accDes
            # ddqdes = accdess[ti,:]
        else:
            u = reactiveController(p, Rb, dq, pdes)
        # u = np.array([1,0.1,0])

        p, Rb, dq = quadrotorNLDyn(p, Rb, dq, u, dtsim, ddq=ddqdes)
        ys[ti,:] = np.hstack((p, Rb[:,2], dq))
        us[ti,:] = u
        pdess[ti,:] = pdes
    
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(4,2)
    ax = ax.ravel()
        
    ax[0].plot(tt, ys[:,:3])
    ax[0].plot(tt, pdess[:,0], 'k--', alpha=0.3)
    ax[0].set_ylabel('p')
    ax[1].plot(tt, ys[:,3:6])
    ax[1].axhline(y=0, color='k', alpha=0.3)
    ax[1].set_ylabel('s')
    ax[2].plot(tt, us[:,0])
    ax[2].axhline(y=0, color='k', alpha=0.3)
    ax[2].set_ylabel('Sp. thrust')
    ax[3].plot(tt, us[:,1:])
    ax[3].axhline(y=0, color='k', alpha=0.3)
    ax[3].set_ylabel('Moments')
    ax[4].plot(tt, ys[:,6:9])
    ax[4].axhline(y=0, color='k', alpha=0.3)
    ax[4].set_ylabel('v')
    ax[5].plot(tt, ys[:,9:12])
    ax[5].axhline(y=0, color='k', alpha=0.3)
    ax[5].set_ylabel('omega')
    ax[6].plot(tt, accdess[:,:3])
    ax[6].axhline(y=0, color='k', alpha=0.3)
    ax[6].set_ylabel('accdes pos')
    ax[7].plot(tt, accdess[:,3:])
    ax[7].axhline(y=0, color='k', alpha=0.3)
    ax[7].set_ylabel('accdes ang')
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    T0 = 0.5
    dt = 5
    N = 3
    s0s = [[0.1,0.1,0.9] for i in range(N)]
    Btaus = [np.full((3,2),1.123) for i in range(N)]
    y0 = np.random.rand(6)
    dy0 = np.random.rand(6)
    g = 9.81e-3

    # weights
    ws = 1e1
    wds = 1e3
    wpr = 1
    wpf = 5
    wvr = 1e3
    wvf = 2e3
    wthrust = 1e-1
    wmom = 1e-2

    ydes = np.zeros_like(y0)
    dydes = np.zeros_like(y0)
    smin = np.array([-2,-2,0.5])
    smax = np.array([2,2,1.5])
    TtoWmax = 2 # thrust-to-weight

    up = UprightMPC2(N, dt, g, smin, smax, TtoWmax, ws, wds, wpr, wpf, wvr, wvf, wthrust, wmom)
    up.testDyn(T0, s0s, Btaus, y0, dy0)

    # # Hover
    # controlTest(up, 500, useMPC=True)
    # Ascent
    controlTest(up, 500, useMPC=True, ascentIC=True)
    # # Traj
    # controlTest(up, 2000, useMPC=True, trajAmp=50, trajFreq=1)

