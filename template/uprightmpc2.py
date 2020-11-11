import osqp
import numpy as np
import scipy.sparse as sp
from scipy.spatial.transform import Rotation
np.set_printoptions(precision=4, suppress=True, linewidth=200)
from genqp import quadrotorNLDyn, skew, Ib
from uprightmpc2py import UprightMPC2C # C version
from time import perf_counter
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import progressbar

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
        xtest = openLoopX(self.N, self.dt, T0, s0s, Btaus, y0, dy0, self.g)
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

def reactiveController(p, Rb, dq, pdes, kpos=[5e-3,5e-1], kz=[1e-1,1e0], ks=[10e0,1e2], **kwargs):
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

def viewControlTestLog(log, log2=None, callShow=True, goal0=False, desTraj=False, vscale=0.4):
    def traj3plot(_ax, t, p, v, cmap, narrow=20):
        cnorm = t/t[-1]
        _ax.scatter(p[:,0], p[:,1], p[:,2], c=cnorm, cmap=cmap, marker='.', label='_nolegend_')
        ii = np.linspace(0, len(t), narrow, dtype=int, endpoint=False)
        v *= vscale
        _ax.quiver(p[ii,0], p[ii,1], p[ii,2], v[ii,0], v[ii,1], v[ii,2], color='b' if "Blue" in cmap else 'r')

    def aspectEqual3(_ax, xyz):
        X, Y, Z = xyz[:,0], xyz[:,1], xyz[:,2]
        # Create cubic bounding box to simulate equal aspect ratio
        max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
        Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
        Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
        Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
        # Comment or uncomment following both lines to test the fake bounding box:
        for xb, yb, zb in zip(Xb, Yb, Zb):
            _ax.plot([xb], [yb], [zb], 'w', label='_nolegend_')

    def posParamPlot(_ax):
        traj3plot(_ax, log['t'], log['y'][:,:3], log['y'][:,3:6], "Blues_r")
        aspectEqual3(_ax, log['y'][:,:3])
        if log2 is not None:
            traj3plot(_ax, log2['t'], log2['y'][:,:3], log2['y'][:,3:6], "Reds_r")
        # _ax.plot(log['t'], log['pdes'][:,0], 'k--', alpha=0.3)
        if goal0:
            _ax.plot([0], [0], [0], 'g*', markersize=10, zorder=10)
            _ax.legend(('MPC', 'Reactive', 'Goal'))
        else:
            _ax.legend(('MPC', 'Reactive'))
        if desTraj:
            _ax.plot(log['pdes'][:,0], log['pdes'][:,1], log['pdes'][:,2], 'k--', alpha=0.5, zorder=9)
        _ax.set_xlabel('x [mm]')
        _ax.set_ylabel('y [mm]')
        _ax.set_zlabel('z [mm]')

    def posPlot(_ax):
        _ax.plot(log['t'], log['y'][:,:3])
        if log2 is not None:
            _ax.plot(log2['t'], log2['y'][:,:3], '--')
        _ax.plot(log['t'], log['pdes'][:,0], 'k--', alpha=0.3)
        _ax.set_ylabel('p')
    def splot(_ax):
        _ax.plot(log['t'], log['y'][:,3:6])
        _ax.axhline(y=0, color='k', alpha=0.3)
        _ax.set_ylabel('s')
    def inputsPlot(_ax1, _ax2):
        _ax1.plot(log['t'], log['u'][:,0])
        _ax1.axhline(y=0, color='k', alpha=0.3)
        _ax1.set_ylabel('Sp. thrust')
        _ax2.plot(log['t'], log['u'][:,1:])
        _ax2.axhline(y=0, color='k', alpha=0.3)
        _ax2.set_ylabel('Moments')
    def velsPlot(_ax1, _ax2):
        _ax1.plot(log['t'], log['y'][:,6:9])
        _ax1.axhline(y=0, color='k', alpha=0.3)
        _ax1.set_ylabel('v')
        if _ax2 is not None:
            _ax2.plot(log['t'], log['y'][:,9:12])
            _ax2.axhline(y=0, color='k', alpha=0.3)
            _ax2.set_ylabel('omega')
    def accdesPlots(_ax1, _ax2):
        _ax1.plot(log['t'], log['accdes'][:,:3])
        _ax1.axhline(y=0, color='k', alpha=0.3)
        _ax1.set_ylabel('accdes pos')
        _ax2.plot(log['t'], log['accdes'][:,3:])
        _ax2.axhline(y=0, color='k', alpha=0.3)
        _ax2.set_ylabel('accdes ang')
    def wlqpuPlots(_ax):
        _ax.plot(log['t'], log['wlqpu'][:,:2])
        _ax.plot(log['t'], log['wlqpu'][:,2:],'--')
        _ax.legend(('0','1','2','3'))
        _ax.set_ylabel('wlqpu')
    
    # fig = plt.figure()
    # ax = [fig.add_subplot(3,3,i+1) for i in range(1,12)]
    # posPlot(ax[0])
    # splot(ax[1])
    # inputsPlot(ax[2], ax[3])
    # velsPlot(ax[4], ax[5])
    # accdesPlots(ax[6], ax[7])
    # wlqpuPlots(ax[8])
    # fig.tight_layout()
    
    fig = plt.figure()
    ax3d = fig.add_subplot(1,1,1,projection='3d')
    posParamPlot(ax3d)

    if callShow:
        plt.show()

def controlTest(mdl, tend, dtsim=0.2, useMPC=True, trajFreq=0, trajAmp=0, ascentIC=False, showPlots=True, tpert=None, speedTest=False, perchTraj=False, flipTask=False, taulim=100, **kwargs):
    """trajFreq in Hz, trajAmp in mm"""
    speedTestvdes = 2 # m/s
    # Initial conditions
    dq = np.zeros(6)
    if ascentIC or speedTest or perchTraj:
        p = np.array([0, 0, -50]) if ascentIC else np.array([-500*speedTestvdes, 0, 0])
        if perchTraj:
            p = np.array([-100, 0, 0])
        Rb = np.eye(3)
    else:
        p = np.array([0, 0, 0])
        Rb = np.eye(3) if flipTask else Rotation.from_euler('xyz', [0.5,-0.5,0]).as_matrix()
        dq[0] = 0.1
    pdes = np.zeros(3)
    dpdes = np.zeros(3)
    sdes = np.array([0,0,1])
    
    tt = np.arange(tend, step=dtsim)
    Nt = len(tt)

    # for logging
    log = {'t': tt, 'y': np.zeros((Nt, 12)), 'u': np.zeros((Nt, 3)), 'pdes': np.zeros((Nt, 3)), 'accdes': np.zeros((Nt,6)), 'wlqpu': np.zeros((Nt,4))}

    trajOmg = 2 * np.pi * trajFreq * 1e-3 # to KHz, then to rad/ms
    ddqdes = None # test integrate ddq sim below

    avgTime = 0.0

    for ti in range(Nt):
        # Traj to follow
        if flipTask:
            # rotation phase 0 to 1
            ph = np.clip((tt[ti] - 100) / 200, 0, 1)
            sdes = np.array([-np.sin(ph*2*np.pi), 0, np.cos(ph*2*np.pi)])
        elif perchTraj:
            if tt[ti] < 500:
                pdes[0] = -100 + 0.2 * tt[ti]
                dpdes[0] = 0.2
                # rotation phase 0 to 1
                ph = np.clip((tt[ti] - 450) / 100, 0, 1)
                sdes = np.array([-np.sin(ph*np.pi), 0, np.cos(ph*np.pi)])
            else:
                pdes[0] = dpdes[0] = 0
                sdes = np.array([-1,0,0])
        elif speedTest:
            if tt[ti] < 500:
                dpdes[0] = speedTestvdes
                pdes[0] = -500*speedTestvdes + speedTestvdes*(tt[ti])
            else:
                pdes[0] = 0
                dpdes[0] = 0
        else:
            pdes[0] = trajAmp * np.sin(trajOmg * tt[ti])
            dpdes[0] = trajAmp * trajOmg * np.cos(trajOmg * tt[ti])
            if trajAmp > 1e-3:
                pdes[2] = 0.1 * tt[ti]
                dpdes[2] = 0.1
                # Add perturbation for this traj
                if tpert is not None and tt[ti] > tpert:
                    dq[1] += 2
                    tpert = None

        # Call controller
        if useMPC:
            t1 = perf_counter()
            # what "u" is depends on w(u). Here in python testing with w(u) = [0,0,u0,u1,u2,u3]
            w0 = np.hstack((0,0,mdl.u0))
            dwdu0 = np.vstack((np.zeros((2,4)), np.eye(4)))
            u, log['accdes'][ti,:], uwlqp = mdl.update(p, Rb, dq, pdes, dpdes, sdes)
            log['wlqpu'][ti,:] = uwlqp
            avgTime += 0.01 * (perf_counter() - t1 - avgTime)
            # # Alternate simulation by integrating accDes
            # ddqdes = accdess[ti,:]
        else:
            u = reactiveController(p, Rb, dq, pdes, **kwargs)
        # u = np.array([1,0.1,0])
        # Input limit
        for i in range(2):
            u[i+1] = np.clip(u[i+1], -taulim, taulim)

        p, Rb, dq = quadrotorNLDyn(p, Rb, dq, u, dtsim, ddq=ddqdes)
        log['y'][ti,:] = np.hstack((p, Rb[:,2], dq))
        log['u'][ti,:] = u
        log['pdes'][ti,:] = pdes
    if useMPC:
        print("Time (ms):", avgTime * 1e3)
    if showPlots:
        viewControlTestLog(log)
    return log

def logMetric(log):
    # A metric to plot about how good the tracking was
    Nt = len(log['t'])
    perr = log['y'][:,:3]
    tau = log['u'][:,1:3]
    serr = log['y'][:,3:6]
    serr[:,2] -= 1.0
    err = 0
    eff = 0
    for i in range(Nt):
        err += np.dot(perr[i,:], perr[i,:])# + 10 * np.dot(serr[i,:], serr[i,:])
        eff += np.dot(tau[i,:], tau[i,:])# + 10 * np.dot(serr[i,:], serr[i,:])
    err /= Nt
    eff /= Nt
    return err, eff

def createMPC(N, **kwargs):
    """Returns the mdl"""
    dt = 5
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
    # WLQP inputs
    mb = 100
    # what "u" is depends on w(u). Here in python testing with w(u) = [0,0,u0,u1,u2,u3].
    # Setting first 2 elements of Qw -> 0 => should not affect objective as longs as dumax does not constrain.
    Qw = np.hstack((np.zeros(2), np.zeros(4)))
    umin = np.array([0, -0.5, -0.2, -0.1])
    umax = np.array([10, 0.5, 0.2, 0.1])
    dumax = np.array([10, 10, 10, 10]) # /s
    controlRate = 1000
    pyver = UprightMPC2(N, dt, g, TtoWmax, ws, wds, wpr, wpf, wvr, wvf, wthrust, wmom, umin, umax, dumax, mb, Ib.diagonal(), Qw, controlRate)
    # C version can be tested too
    popts = np.zeros(90)
    cver = UprightMPC2C(dt, g, TtoWmax, ws, wds, wpr, wpf, wvr, wvf, wthrust, wmom, mb, Ib.diagonal(), umin, umax, dumax, Qw, controlRate, 50, popts)
    return pyver, cver

def papPlots():
    # # Flip traj ---------------------
    # l1 = controlTest(up, 1000, useMPC=True, showPlots=False, flipTask=True)
    # viewControlTestLog(l1, desTraj=True, vscale=10)
    # fig, ax = plt.subplots(1,3, figsize=(7.5,2.5))
    # for i in range(0,3,2):
    #     ax[i].plot(1e-3*l1['t'], l1['y'][:,i], 'b')
    #     ax[i].plot(1e-3*l1['t'], l1['pdes'][:,i], 'k--', alpha=0.3)
    #     ax[i].set_xlabel('t [s]')
        
    # ax[1].plot(1e-3*l1['t'], 180/np.pi*np.arctan2(l1['y'][:,3], l1['y'][:,5]), 'b')
    # ax[1].plot([0, 0.1, 0.2], [0, 0, -180], 'k--', alpha=0.3)
    # ax[1].plot([0.2, 0.3, 1], [180, 0, 0], 'k--', alpha=0.3)
    # ax[1].set_ylabel('Angle [deg]')
    # ax[0].set_ylabel('x [mm]')
    # ax[2].set_ylabel('z [mm]')
    # fig.tight_layout()
    # plt.show()

    # # Perch traj ---------------------
    # l1 = controlTest(up, 550, useMPC=True, showPlots=False, perchTraj=True)
    # # viewControlTestLog(l1, desTraj=True, vscale=10)
    # fig, ax = plt.subplots(1,3, figsize=(7.5,2.5))
    # for i in range(0,3,2):
    #     ax[i].plot(1e-3*l1['t'], l1['y'][:,i], 'b')
    #     ax[i].plot(1e-3*l1['t'], l1['pdes'][:,i], 'k--', alpha=0.3)
    #     ax[i].set_xlabel('t [s]')
        
    # ax[1].plot(1e-3*l1['t'], 180/np.pi*np.arctan2(l1['y'][:,3], l1['y'][:,5]), 'b')
    # ax[1].plot([0, 0.45, 0.55], [0, 0, -90], 'k--', alpha=0.3)
    # ax[1].set_ylabel('Angle [deg]')
    # ax[0].set_ylabel('x [mm]')
    # ax[2].set_ylabel('z [mm]')
    # fig.tight_layout()
    # plt.show()

    def hoverTask(show3d, reactiveArgs1, reactiveArgs2=None):
        l1 = controlTest(up, 300, useMPC=True, showPlots=False)
        l2 = controlTest(up, 1000, useMPC=False, showPlots=False, **reactiveArgs1)
        if reactiveArgs2 is not None:
            l3 = controlTest(up, 1000, useMPC=False, showPlots=False, **reactiveArgs2)
        if show3d:
            viewControlTestLog(l1, log2=l2, goal0=True)
        else:
            fig, ax = plt.subplots(1,2, figsize=(5,2.5))
            for i in range(2):
                ax[i].plot(1e-3*l1['t'], l1['y'][:,i], 'b')
                ax[i].plot(1e-3*l2['t'], l2['y'][:,i], 'r')
                if reactiveArgs2 is not None:
                    ax[i].plot(1e-3*l3['t'], l3['y'][:,i], 'r--')
                ax[i].plot(1e-3*l2['t'], l2['pdes'][:,i], 'k--', alpha=0.3)
                ax[i].set_xlabel('t [s]')
            ax[0].set_ylabel('x [mm]')
            ax[1].set_ylabel('y [mm]')
            fig.tight_layout()
            plt.show()

    def sTask(show3d, **reactiveArgs):
        l1 = controlTest(up, 2000, useMPC=True, showPlots=False, trajAmp=50, trajFreq=1, tpert=1000)
        l2 = controlTest(up, 2000, useMPC=False, showPlots=False, trajAmp=50, trajFreq=1, tpert=1000, **reactiveArgs)
        if show3d:
            viewControlTestLog(l1, log2=l2, desTraj=True, vscale=20)
        else:
            fig, ax = plt.subplots(1,2, figsize=(5,2.5))
            for i in range(2):
                ax[i].plot(1e-3*l1['t'], 1e-3*l1['y'][:,i], 'b')
                ax[i].plot(1e-3*l2['t'], 1e-3*l2['y'][:,i], 'r')
                ax[i].plot(1e-3*l2['t'], 1e-3*l2['pdes'][:,i], 'k--', alpha=0.3)
                ax[i].set_xlabel('t [s]')
            ax[0].set_ylabel('x [m]')
            ax[1].set_ylabel('y [m]')
            fig.tight_layout()
            plt.show()

    # # Straight line acceleration -------
    # l1 = controlTest(up, 1000, useMPC=True, showPlots=False, speedTest=True)
    # l2 = controlTest(up, 1000, useMPC=False, showPlots=False, speedTest=True)
    # viewControlTestLog(l1, log2=l2, desTraj=True, vscale=50)
    # fig, ax = plt.subplots(1,2, figsize=(5,2.5))
    # ax[0].plot(1e-3*l1['t'], 1e-3*l1['y'][:,0], 'b')
    # ax[0].plot(1e-3*l2['t'], 1e-3*l2['y'][:,0], 'r')
    # ax[0].plot(1e-3*l2['t'], 1e-3*l2['pdes'][:,0], 'k--', alpha=0.3)
    # ax[0].set_ylabel('x [m]')
    # ax[0].set_xlabel('t [s]')
    # ax[1].plot(1e-3*l1['t'], l1['y'][:,6], 'b')
    # ax[1].plot(1e-3*l2['t'], l2['y'][:,6], 'r')
    # ax[1].set_ylabel('xdot [m/s]')
    # ax[1].set_xlabel('t [s]')
    # fig.tight_layout()
    # plt.show()

    # Hover tuning ---------
    def gainTuningReactiveSims(kwgain, k1range, k2range, kwfixedn, kwfixedv, npts=10):
        k1s = np.linspace(*k1range,num=npts)
        k2s = np.linspace(*k2range,num=npts)
        xv, yv = np.meshgrid(k1s, k2s, indexing='ij') # treat xv[i,j], yv[i,j]
        costs = np.zeros_like(xv)
        efforts = np.zeros_like(xv)
        # create a progress bar
        widgets = [
            'Progress: ', progressbar.Percentage(),
            ' ', progressbar.Bar(),
            ' ', progressbar.ETA(),
        ]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=np.prod(costs.shape))
        nrun = 0
        for i in range(len(k1s)):
            for j in range(len(k2s)):
                nrun += 1
                bar.update(nrun)
                try:
                    kwargs = {kwgain: [xv[i,j],yv[i,j]], kwfixedn: kwfixedv}
                    l2 = controlTest(up, 1000, useMPC=False, showPlots=False, taulim=10, **kwargs)
                    costs[i,j], efforts[i,j] = logMetric(l2)
                except KeyboardInterrupt:
                    raise
                except:
                    costs[i,j] = efforts[i,j] = np.nan
        np.savez(kwgain+str('.npz'), xv=xv, yv=yv, costs=costs, efforts=efforts)

    def gainTuningReactivePlots(maxcost=10):
        lmpc = controlTest(up, 1000, useMPC=True, showPlots=False)
        empc, effmpc = logMetric(lmpc)
                        
        def plot1(ax, dat):
            costs = np.clip(dat['costs'] / empc, 0, maxcost)
            im = ax.pcolormesh(dat['xv'], dat['yv'], costs, cmap='gray_r', shading='auto')
            fig.colorbar(im, ax=ax)
            
        fig, ax = plt.subplots(1,2,figsize=(9,4))
        plot1(ax[0], np.load('ks.npz'))
        ax[0].plot([15], [100], 'r*', ms=20)
        plot1(ax[1], np.load('kpos.npz'))
        ax[1].plot([0.01, 0.04], [1.0, 1.25], 'r*', ms=20)
        plt.show()

    def trackingEffortPlot(ffs):
        # Baseline
        lmpc = controlTest(up, 1000, useMPC=True, showPlots=False)
        empc, effmpc = logMetric(lmpc)
        costs2 = []
        effs2 = []
        for ff in ffs:
            dat = np.load(ff)
            costs = dat['costs'].ravel() / empc
            effs = dat['efforts'].ravel() / effmpc
            ii = np.where(costs < 10)[0]
            costs2.append(costs[ii])
            effs2.append(effs[ii])
        fig, ax = plt.subplots(1, figsize=(4,4))
        ax.scatter(costs2, effs2, color='r')
        ax.axhline(1, color='k', linestyle='dashed', alpha=0.3)
        ax.axvline(1, color='k', linestyle='dashed', alpha=0.3)
        ax.set_xlim((0,10))
        ax.set_ylim((0,10))
        ax.set_aspect('equal')
        ax.set_xlabel('Relative tracking error [ ]')
        ax.set_ylabel('Relative actuator effort [ ]')
        plt.show()

    # # Run and save data
    # # defaults kpos=[5e-3,5e-1], kz=[1e-1,1e0], ks=[10e0,1e2]
    # gainTuningReactiveSims('ks', [5e0,2e1], [2e1,2e2], 'kpos', [5e-3,5e-1])
    # gainTuningReactiveSims('kpos', [1e-3,8e-2], [1e-1,2e0], 'ks', [15,100])

    # gainTuningReactivePlots()

    # hoverTask(False, {'ks':[15,100], 'kpos':[0.01,1]}, {'ks':[15,100], 'kpos':[0.04,1.25]})
    # sTask(False, ks=[15,100], kpos=[0.01,1])

    trackingEffortPlot(['kpos.npz'])

if __name__ == "__main__":
    N = 3
    T0 = 0.5
    s0s = [[0.1,0.1,0.9] for i in range(N)]
    Btaus = [np.full((3,2),1.123) for i in range(N)]
    y0 = np.random.rand(6)
    dy0 = np.random.rand(6)

    ydes = np.zeros_like(y0)
    dydes = np.zeros_like(y0)
    TtoWmax = 2 # thrust-to-weight

    up, upc = createMPC(N)
    up.testDyn(T0, s0s, Btaus, y0, dy0)

    # # FIXME: test
    # p = np.random.rand(3)
    # R = np.random.rand(3, 3)
    # dq = np.random.rand(6)
    # pdes = np.random.rand(3)
    # dpdes = np.random.rand(3)
    # retc = upc.update(p, R, dq, pdes, dpdes)
    # cl, cu, cq = upc.vectors()
    # cP, cAdata, cAidx = upc.matrices()
    # ret = up.update(p, R, dq, pdes, dpdes)
    # # print(cAdata - up.A.data[cAidx])
    # print(ret[0], ret[1], ret[0]-retc[0], ret[1]-retc[1])

    # # Hover
    # controlTest(up, 500, useMPC=True)
    # # Ascent
    # controlTest(up, 500, useMPC=True, ascentIC=True)
    # # Traj
    # controlTest(up, 2000, useMPC=True, trajAmp=50, trajFreq=1)

    papPlots()