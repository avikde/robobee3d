import autograd.numpy as np
from autograd import jacobian, hessian
import scipy.sparse as sparse
import matplotlib.pyplot as plt
from matplotlib import animation
import sys, time
sys.path.append('..')
from controlutils.py.model import Model
import controlutils.py.ltvsystem as ltvsystem

class Wing2DOF(Model):
    nx = 4
    nu = 1
    y0 = np.zeros(nx)
    
    rescale = 1.0

    def aero(self, y, u, _params=[]):
        cbar = _params[0]
        CLmax = 1.8
        CDmax = 3.4
        CD0 = 0.4
        rho = 1.225
        R = 15e-3
        
        sigma, psi, dsigma, dpsi = tuple(y)
        cpsi = np.cos(psi)
        spsi = np.sin(psi)
        alpha = psi

        # aero force
        wing1 = np.array([sigma, 0])
        paero = wing1 + np.array([[cpsi, -spsi], [spsi, cpsi]]) @ np.array([0, -cbar])
        Jaero = np.array([[1, cbar * cpsi], [0, cbar * spsi]])
        CL = CLmax * np.sin(2 * alpha)
        CD = (CDmax + CD0)/2 - (CDmax - CD0)/2 * np.cos(2 * alpha)
        vaero = np.array([dsigma, 0])
        # TODO: confirm and expose design params as argument
        Faero = 1/2 * rho * cbar * R * (vaero.T @ vaero) * np.array([CD, CL]) * np.sign(-dsigma)

        return paero, Jaero, Faero

    def dydt(self, yin, u, _params=[]):
        ''' 
        See mma file flapping wing traj
        '''
        T = _params[1]
        Kscale = np.diag([self.rescale, 1, self.rescale, 1])
        y = np.linalg.inv(Kscale) @ yin
        # NOTE: for optimizing transmission ratio
        # Thinking of y = (sigma_actuator, psi, dsigma_actuator, dpsi)
        # u = (tau_actuator)
        # sigma = sigma_actuator * T; tau = tau_actuator / T
        sigma = y[0] * T
        psi = y[1]
        dsigma = y[2] * T
        dpsi = y[3]
        cpsi = np.cos(psi)
        spsi = np.sin(psi)

        # params
        cbar = _params[0]
        mspar = 0
        ka = 0
        khinge = 1e-3
        mwing = 5e-6
        Iwing = 1e-9#mwing * cbar**2
        bpsi = 5e-7

        # inertial terms
        M = np.array([[mspar + mwing, cbar * mwing * cpsi], [cbar * mwing * cpsi, Iwing + cbar**2 * mwing]])
        corgrav = np.array([ka * sigma - cbar * mwing * spsi * dpsi**2, khinge * psi])
        # non-lagrangian terms
        taudamp = np.array([0, -bpsi * dpsi])
        _, Jaero, Faero = self.aero(y, u, params)
        tauaero = Jaero.T @ Faero
        # input
        tauinp = np.array([u[0] / T, 0])

        ddq = np.linalg.inv(M) @ (-corgrav + taudamp + tauaero + tauinp)

        return Kscale @ np.hstack((y[2:], ddq))

    @property
    def limits(self):
        # This is based on observing the OL trajectory
        umin = np.array([-0.15])
        umax = -umin
        xmin = np.array([-0.02 * self.rescale, -1.2, -np.inf, -np.inf])
        xmax = -xmin
        return umin, umax, xmin, xmax


def flapkin(yui, xyoff, _params):
    """Returns wing positions and stuff for visualization"""
    # wing extents
    wing1 = np.array([yui[0], 0]) + np.asarray(xyoff)
    c, s = np.cos(yui[1]), np.sin(yui[1])
    wing2 = wing1 + np.array([[c, -s], [s, c]]) @ np.array([0, -2*_params[0]])
    # aero arrow extents
    _, _, Faero = m.aero(yui[:m.nx], yui[m.nx:], _params)
    pcop = (wing1 + wing2)/2
    aeroEnd = pcop + 0.1 / m.rescale * Faero
    return wing1, wing2, pcop, aeroEnd

def flapVisUpdate(yui, xyoff, _params, plwing, plaero):
    wing1, wing2, pcop, aeroEnd = flapkin(yui, xyoff, _params)
    plwing.set_xdata([wing1[0], wing2[0]])
    plwing.set_ydata([wing1[1], wing2[1]])
    plaero.set_xdata([pcop[0], aeroEnd[0]])
    plaero.set_ydata([pcop[1], aeroEnd[1]])

# Global
m = Wing2DOF()
m.rescale = 30.0
# params: [cbar, T(ransmission ratio), ]
params = np.array([5e-3, 1])
# --------------------------------------

# Wing traj opt using QP -------------------------------------------------
def dirTranForm(xtraj, N, nx, nu):
    # convert from the (N-1,nx+nu) array to the dirtran form with x0,..,x[N-1],u0,...,u[N-2] i.e. shorter by 1
    return np.hstack((
            np.reshape(xtraj[:,:nx], (N+1)*nx, 'C'),
            np.reshape(xtraj[:-1,nx:], (N)*nu, 'C')
        ))
def xuMatForm(dirtranx, N, nx):
    Xmat = np.reshape(dirtranx[:(N)*nx], (N,nx), 'C')
    Umat = dirtranx[(N)*nx:][:,np.newaxis]
    # Need to add another u since dirtran form has one fewer u. Just copy the last one
    Umat = np.vstack((Umat, Umat[-1,:]))
    return np.hstack((Xmat, Umat))

# 

def Jobjinst(y, u, params):
    _, _, Faero = m.aero(y, u, params)
    return -Faero[1] # minimization

def Jcost_dirtran(dirtranx, N, params):
    '''this is over a traj, no penalty'''
    c = 0
    for k in range(N):
        c += Jobjinst(dirtranx[(k*m.nx):((k+1)*m.nx)], dirtranx[((N+1)*m.nx + k*m.nu):((N+1)*m.nx + (k+1)*m.nu)], params)
    return c


class QOFAvgLift:
    def __init__(self, N, wx, wu, kdampx, kdampu):
        self.N = N
        self.nx = len(wx)
        self.nu = len(wu)
        self.wx = wx
        self.wu = wu
        self.kdampx = kdampx
        self.kdampu = kdampu
        # autodiff to get gradient of avg lift cost
        Jx = lambda x : Jcost_dirtran(x, N, params)
        self.DJfun = jacobian(Jx)
        self.D2Jfun = hessian(Jx)

    def getPq(self, xtraj):
        dirtranx = dirTranForm(xtraj, self.N, self.nx, self.nu)
        nX = (self.N+1) * self.nx + self.N*self.nu
        # vector of weights for the whole dirtran x
        w = np.hstack((np.tile(self.wx, self.N+1), np.tile(self.wu, self.N)))
        kdamp = np.hstack((np.tile(self.kdampx, self.N+1), np.tile(self.kdampu, self.N)))

        # Only regularization
        self.P = sparse.diags(w + kdamp).tocsc()
        self.q = -kdamp * dirtranx

        # # Test: weight hinge angle
        # wk = np.array([0,1000,0,0])
        # w = np.hstack((np.tile(wk, self.N+1), np.zeros(self.nu * self.N)))
        # self.P = sparse.diags(w).tocsc()
        # self.q = np.zeros_like(dirtranx)
        # # self.q = - w * (np.pi / 4)

        # # Test: weight stroke pos
        # wk = np.array([1000000,0,0,0])
        # w = np.hstack((np.tile(wk, self.N+1), np.zeros(self.nu * self.N)))
        # self.P = sparse.diags(w).tocsc()
        # self.q = np.zeros_like(dirtranx)
        # self.q = - w * (np.pi / 4)

        # Second order approx of aero force
        DJ = self.DJfun(dirtranx)
        # # D2J = sparse.csc_matrix(np.outer(DJ, DJ))
        # D2J = sparse.csc_matrix(self.D2Jfun(dirtranx))
        # self.P = D2J
        # # self.P += sparse.diags(w + kdamp).tocsc()
        # self.q = DJ - D2J @ dirtranx
        # # self.q -= kdamp * dirtranx

        # First-order: make objective to *reduce" cost instead of find optimum. Like value function vs. optimal cost.
        # min (J(x) - J(x0)) ~= DJ(x0) (x - x0) ~~ DJ(x0) x
        self.q += DJ

        return self.P, self.q
    
    def cost(self, xtraj):
        self.getPq(xtraj)
        dirtranx = dirTranForm(xtraj, self.N, self.nx, self.nu)
        return 0.5 * dirtranx @ self.P @ dirtranx + self.q @ dirtranx

class WingQP:
    def __init__(self, model, N, wx, wu, kdampx, kdampu, **settings):
        self.ltvsys = ltvsystem.LTVSolver(model)
        # Dynamics and constraints
        self.ltvsys.initConstraints(model.nx, model.nu, N, periodic=True, 
        polyBlocks=None)
        # Add N+1th row to xtraj?
        self.ltvsys.xtraj = np.vstack((self.ltvsys.xtraj, self.ltvsys.xtraj[-1,:]))
        # Back to standard
        self.ltvsys.initObjective(QOFAvgLift(N, wx, wu, kdampx, kdampu))
        self.ltvsys.initSolver(**settings)

    def update(self, xtraj):
        N = self.ltvsys.N
        nx = self.ltvsys.nx
        # TODO: check which traj mode
        u0 = xtraj[:,4][:,np.newaxis]
        # NOTE: confirmed that updateTrajectory correctly updates the traj, and that updateDynamics is updating the A, B
        xtraj = self.ltvsys.updateTrajectory(xtraj[:,:4], u0, ltvsystem.GIVEN_POINT_OR_TRAJ, params)
        self.ltvsys.updateObjective()
        self.dirtranx, res = self.ltvsys.solve(throwOnError=False)
        if res.info.status not in ['solved', 'solved inaccurate', 'maximum iterations reached']:
            self.ltvsys.debugResult(res)
            # dirtranx = dirTranForm(xtraj, N, nx, self.ltvsys.nu)
            # print(self.ltvsys.u - self.ltvsys.A @ dirtranx, self.ltvsys.A @ dirtranx - self.ltvsys.l)

            raise ValueError(res.info.status)
        # debug
        # print(self.ltvsys.u - self.ltvsys.A @ dirtranx, self.ltvsys.A @ dirtranx - self.ltvsys.l)
        # reshape into (N,nx+nu)
        return xuMatForm(self.dirtranx, N+1, nx)
        
    # Debug the solution
    def _dt(self, traj):
        return dirTranForm(traj, self.ltvsys.N, self.ltvsys.nx, self.ltvsys.nu) if len(traj.shape) > 1 else traj

    def debugConstraintViol(self, *args):
        """args are trajs in dirtran form or square"""
        # row index where the dynamics constraints end
        ridyn = (self.ltvsys.N + 1) * self.ltvsys.nx
        rixlim = ridyn + (self.ltvsys.N + 1) * self.ltvsys.nx
        riulim = rixlim + self.ltvsys.N * self.ltvsys.nu
        _, ax = plt.subplots(2)
        for i in range(len(args)):
            ax[0].plot(self.ltvsys.A @ self._dt(args[i]) - self.ltvsys.l, '.', label=str(i))
        ax[0].legend()
        for i in range(len(args)):
            ax[1].plot(self.ltvsys.u - self.ltvsys.A @ self._dt(args[i]), '.', label=str(i))
        for i in range(2):
            ax[i].axhline(0, color='k', alpha=0.3) # mark 0
            for j in range(0, ridyn, self.ltvsys.nx):
                ax[i].axvline(j, color='k', alpha=0.1)
            ax[i].axvline(ridyn, color='b', alpha=0.3)
            for j in range(ridyn, rixlim, self.ltvsys.nx):
                ax[i].axvline(j, color='k', alpha=0.1)
            ax[i].axvline(rixlim, color='r', alpha=0.3)
            for j in range(rixlim, riulim, self.ltvsys.nu):
                ax[i].axvline(j, color='k', alpha=0.1)
            ax[i].axvline(riulim, color='g', alpha=0.3)
        ax[1].legend()

    def plotTrajs(self, *args):
        """Helper function to plot a bunch of trajectories superimposed"""
        umin, umax, xmin, xmax = self.ltvsys.m.limits
        _, ax = plt.subplots(3)
        for arg in args:
            ax[0].plot(self.trajt, arg[:,0], '.-')
        for yy in [xmin[0], xmax[0], 0]:
            ax[0].axhline(y=yy, color='k', alpha=0.3)
        ax[0].set_ylabel('stroke (m)')
        for arg in args:
            ax[1].plot(self.trajt, arg[:,1], '.-')
        for yy in [xmin[1], xmax[1], np.pi/4, -np.pi/4]:
            ax[1].axhline(y=yy, color='k', alpha=0.3)
        ax[1].set_ylabel('hinge angle (rad)')
        for arg in args:
            ax[2].plot(self.trajt, arg[:,4], '.-')
        ax[2].axhline(y=umin[0], color='k', alpha=0.3)
        ax[2].axhline(y=umax[0], color='k', alpha=0.3)
        ax[2].set_ylabel('stroke force (N)')
        for arg in args:
            print('cost = ', self.ltvsys.qof.cost(arg))

""" Penalty method 
Inspired by Coros et al
-----------------------------------------------------------------------
"""

def Ind(x, eps):
    """Indicator function to use inequality constraints in a penalty method.
    Following Geilinger et al (2018) skaterbots. Ref Bern et al (2017).
    
    A scalar inequality f(x) <= b should be modeled as Ind(f(x) - b) and added to the cost."""
    if x <= -eps:
        return 0
    elif x > -eps and x < eps:
        return x**3/(6*eps) + x**2/2 + eps*x/2 + eps**2/6
    else:
        return x**2 + eps**2/3

def Jcosttraj_penalty(dirtranx, N, params, opt={}):
    '''this is over a traj. yu = (nx+nu,Nt)-shaped'''
    # Get all the relevant options from the dict
    PENALTY_DYNAMICS = opt.get('dynamics', 1e-6)
    PENALTY_PERIODIC = opt.get('periodic', 0)
    PENALTY_ULIM = opt.get('input', 0)
    PENALTY_XLIM = opt.get('state', 0)
    PENALTY_EPS = opt.get('eps', 0.1)
    OBJ_LIFT = opt.get('olift', 1)
    OBJ_DRAG = opt.get('odrag', 0)
    OBJ_MOM = opt.get('omom', 0)

    c = 0
    ykfun = lambda k : dirtranx[(k*m.nx):((k+1)*m.nx)]
    ukfun = lambda k : dirtranx[((N+1)*m.nx + k*m.nu):((N+1)*m.nx + (k+1)*m.nu)]
    # Objective
    for i in range(N):
        paero, _, Faero = m.aero(ykfun(i), ukfun(i), params)
        c += -OBJ_LIFT * Faero[1]
        c += OBJ_DRAG * Faero[0]
        c += OBJ_MOM * (-paero[0] * Faero[1] + paero[1] * Faero[0]) # moment

    # Dynamics constraint
    for i in range(N-1):
        dynErr = ykfun(i+1) - (ykfun(i) + m.dydt(ykfun(i), ukfun(i), params) * m.dt)
        c += PENALTY_DYNAMICS * dynErr.T @ dynErr

    # Periodicity
    periodicErr = ykfun(N) - ykfun(0)
    c += PENALTY_PERIODIC * periodicErr.T @ periodicErr

    # Inequality constraint for input limit
    umin, umax, ymin, ymax = m.limits
    Indv = np.vectorize(lambda x : Ind(x, PENALTY_EPS))
    for i in range(N-1):
        yi = ykfun(i)
        ui = ukfun(i)
        c += PENALTY_ULIM * np.sum(Indv(ui - umax) + Indv(-ui + umin))
        c += PENALTY_XLIM * np.sum(Indv(yi - ymax) + Indv(-yi + ymin))

    return c

class WingPenaltyOptimizer:
    """Works with dirtran form of x only"""
    WRT_TRAJ = 0
    WRT_PARAMS = 1
    WRT_TRAJ_PARAMS = 2

    def __init__(self, N):
        self.N = N
        self._Nx = (self.N+1) * m.nx + self.N*m.nu #dirtran size
    
    def update(self, traj0, params0, mode=WRT_TRAJ, **kwargs):
        start = time.perf_counter()
        # Some error checking
        assert len(traj0) == self._Nx
        assert len(params0) == len(params)

        if mode == self.WRT_PARAMS:
            J = lambda p : Jcosttraj_penalty(traj0, self.N, p, **kwargs) # wrt params
            x0 = params0
        elif mode == self.WRT_TRAJ:
            J = lambda traj : Jcosttraj_penalty(traj, self.N, params0, **kwargs)
            x0 = traj0
        elif mode == self.WRT_TRAJ_PARAMS:
            J = lambda trajp : Jcosttraj_penalty(trajp[:self._Nx], self.N, trajp[self._Nx:], **kwargs)
            x0 = np.hstack((traj0, params0))
        else:
            raise ValueError('Invalid mode')
        DJ = jacobian(J)
        D2J = hessian(J)
        J0 = J(x0)
        DJ0 = DJ(x0)
        D2J0 = D2J(x0)

        # descent direction
        # v = -DJ0 # gradient descent
        # Newton's method followed by backtracking line search http://www.stat.cmu.edu/~ryantibs/convexopt-S15/lectures/14-newton.pdf
        try:
            # regularization
            D2J0 += 1e-3 * np.eye(D2J0.shape[0])
            v = -np.linalg.inv(D2J0) @ DJ0
        except np.linalg.LinAlgError:
            # last u (last diag elem) is 0 - makes sense
            print(np.linalg.eigs(D2J0))
            raise
        # search for s
        alpha = 0.4
        beta = 0.9
        s = 1
        while J(x0 + s * v) > J0 + alpha * s * DJ0.T @ v:
            s = beta * s
        # debugging
        print("{:.3f}s; cost {:.1f} -> {:.1f}".format(time.perf_counter() - start, J0, J(x0 + s * v)))
        # perform Newton update
        return x0 + s * v
        
    def plotTrajs(self, *args):
        """Helper function to plot a bunch of trajectories superimposed"""
        umin, umax, xmin, xmax = m.limits
        trajt = range(self.N)
        yend = (self.N) * m.nx # N to ignore the last one
        ustart = (self.N+1) * m.nx
        uend = ustart + self.N*m.nu
        _, ax = plt.subplots(3)
        for arg in args:
            ax[0].plot(trajt, arg[0:yend:m.nx], '.-')
        for yy in [xmin[0], xmax[0], 0]:
            ax[0].axhline(y=yy, color='k', alpha=0.3)
        ax[0].set_ylabel('stroke (m)')
        for arg in args:
            ax[1].plot(trajt, arg[1:yend:m.nx], '.-')
        for yy in [xmin[1], xmax[1], np.pi/4, -np.pi/4]:
            ax[1].axhline(y=yy, color='k', alpha=0.3)
        ax[1].set_ylabel('hinge angle (rad)')
        for arg in args:
            ax[2].plot(trajt, arg[ustart:uend:m.nu], '.-')
        ax[2].axhline(y=umin[0], color='k', alpha=0.3)
        ax[2].axhline(y=umax[0], color='k', alpha=0.3)
        ax[2].set_ylabel('stroke force (N)')

# Create "cts" trajectories from traj (control) knot points ----

def knotPointControl(t, y, traj, ttraj):
    # dumb TODO: interpolate
    u0 = [0] # which knot point input to use
    if len(traj.shape) == 1:
        # deduce N
        N = (len(traj) - m.nx) // (m.nx + m.nu)
        ustart = (N+1)*m.nx
    for i in range(len(ttraj)):
        if ttraj[0] + t > ttraj[i]:
            if len(traj.shape) > 1:
                u0 = traj[i,m.nx:] # xu mat form
            else:
                u0 = traj[(ustart+i*m.nu):(ustart+(i+1)*m.nu)] # dirtran form
    return m.dydt(y, u0, params)

def createCtsTraj(dt, ttrajs, trajs):
    from scipy.integrate import solve_ivp
    # Create shorter timestep sims
    ctstrajs = []
    tvec = np.arange(0, ttrajs[-1] - ttrajs[0], dt)
    tf = tvec[-1]
    if len(trajs[0].shape) > 1:
        y0 = trajs[0][0,:m.nx] # xu mat form
    else:
        y0 = trajs[0][:m.nx] # dirtran form
    for opttraj in trajs:
        # Sim of an openloop controller
        sol = solve_ivp(lambda t, y: knotPointControl(t, y, opttraj, ttrajs), [0, tf], y0, dense_output=True, t_eval=tvec)
        ctstrajs.append(sol.y.T)
    return tvec, ctstrajs


# Animation -------------------------

def _init():
    global _plwings
    global _plaeros
    return tuple(_plwings + _plaeros)

def _animate(i):
    global _plwings
    global _plaeros
    global _xyoffs
    global _trajs
    for k in range(len(_plwings)):
        flapVisUpdate(_trajs[k][i,:], _xyoffs[k], params, _plwings[k], _plaeros[k])
    return tuple(_plwings + _plaeros)

def trajAnim(tvec, ctstrajs, save=False, fps=30):
    global _plwings
    global _plaeros
    global _xyoffs
    global _trajs
    global anim #The object created by FuncAnimation must be assigned to a global variable

    fig, _ax = plt.subplots()

    _trajs = ctstrajs
    _xyoffs = [[0, 0.05 * i] for i in range(len(_trajs))]
    _plwings = [_ax.plot([], [], 'b.-', linewidth=4)[0] for i in range(len(_trajs))]
    _plaeros = [_ax.plot([], [], 'r', linewidth=2)[0] for i in range(len(_trajs))]

    _ax.grid(True)
    _ax.set_aspect(1)
    _ax.set_xlim([-1, 1])
    _ax.set_ylim([-0.05, 0.05 * len(_trajs)])

    anim = animation.FuncAnimation(fig, _animate, init_func=_init, frames=len(tvec), interval=1000/fps, blit=True)
    if save:
        # Set up formatting for the movie files
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)
        import time
        timestamp = time.strftime('%Y%m%d%H%M%S', time.localtime())
        anim.save('trajOptWing2DOF'+timestamp+'.mp4', writer=writer)

    plt.tight_layout()

