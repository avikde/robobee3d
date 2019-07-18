import autograd.numpy as np
from autograd import jacobian, hessian
import scipy.sparse as sparse
import matplotlib.pyplot as plt
from matplotlib import animation
import sys, time
sys.path.append('..')
from controlutils.py.model import Model
import controlutils.py.ltvsystem as ltvsystem

class DoubleIntegrator(Model):
    """Test simple model for making sure opt wrt params works"""
    def dydt(self, yin, u, _params=[]):
        T = _params[0]
        ddq = tauinp = np.array([u[0] / T])
        return np.hstack((y[1:] * T, ddq))

    @property
    def limits(self):
        # This is based on observing the OL trajectory
        umin = np.array([-0.15])
        umax = -umin
        xmin = np.array([-0.02, -np.inf])
        xmax = -xmin
        return umin, umax, xmin, xmax

class Wing2DOF(Model):
    nx = 4
    nu = 1
    y0 = np.zeros(nx)
    
    rescale = 1.0

    def aero(self, y, u, _params=[]):
        cbar = _params[0]
        T = _params[1]
        CLmax = 1.8
        CDmax = 3.4
        CD0 = 0.4
        rho = 1.225
        R = 15e-3
        
        sigma = y[0] * T
        psi = y[1]
        dsigma = y[2] * T
        dpsi = y[3]
        cpsi = np.cos(psi)
        spsi = np.sin(psi)
        alpha = np.pi / 2 - psi

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
        cbar = _params[0]
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

        return Kscale @ np.array([dsigma, dpsi, ddq[0], ddq[1]])

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
params = np.array([5e-3, 1.5])
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

def Jcosttraj_penalty(traj, N, params, opt={}):
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
    PENALTY_ALLOW = opt.get('pen', True)
    PENALTY_H = opt.get('timestep', (0, 0.0001, 0.005)) # tuple of (weight, min, max)

    c = 0
    ykfun = lambda k : traj[(k*m.nx):((k+1)*m.nx)]
    ukfun = lambda k : traj[((N+1)*m.nx + k*m.nu):((N+1)*m.nx + (k+1)*m.nu)]
    h = traj[-1]  #timestep
    # Objective
    Favg = np.zeros(2)
    for i in range(N):
        paero, _, Faero = m.aero(ykfun(i), ukfun(i), params)
        Favg += Faero
    # Use the avg aero force for the objective
    c += -OBJ_LIFT * np.sqrt(Favg[1]**2 + 1e-10) #* (1/30e-3) / h
    # Minimize avg drag over a flap TODO: desired force
    c += OBJ_DRAG * np.sqrt(Favg[0]**2 + 1e-10)
    # FIXME: avg lift not working
    # c *= (1/30e-3) / h
    c += 1e3 * h**2

    # c += OBJ_DRAG * Favg[0]
    # c += OBJ_MOM * (-paero[0] * Faero[1] + paero[1] * Faero[0]) # moment
    # For the objectives, want "average", i.e. divide by the total time of the traj = h * N. Leaving out the N (it is constant): the only difference it makes is to the penalty coefficients.
    # c *= m.dt / h # *initial dt in order to not return the penalties
    # c += 1e6 * h

    # Inequality constraint for input limit
    umin, umax, ymin, ymax = m.limits
    Indv = lambda v : [Ind(vj, PENALTY_EPS) for vj in v] # do this instead of vectorize which gives weird errors
    for i in range(N-1):
        c += PENALTY_ULIM * np.sum(Indv(ukfun(i) - umax) + Indv(-ukfun(i) + umin))
        c += PENALTY_XLIM * np.sum(Indv(ykfun(i) - ymax) + Indv(-ykfun(i) + ymin))
    # Limits on the timestep
    hmin, hmax = PENALTY_H[1], PENALTY_H[2]
    c += PENALTY_H[0] * (Ind(1e6 * (h - hmax), PENALTY_EPS) + Ind(1e6 * (-h + hmin), PENALTY_EPS)) # Use us units here
        
    # Quadratic terms handled separately since we can use a Gauss-Newton approx
    # Dynamics constraint
    rs = [np.sqrt(PENALTY_DYNAMICS) * (ykfun(i+1) - (ykfun(i) + h * m.dydt(ykfun(i), ukfun(i), params))) for i in range(N-1)]
    # Periodicity
    rs.append(np.sqrt(PENALTY_PERIODIC) * (ykfun(N) - ykfun(0)))
    # stack into a vector
    r = np.hstack(rs)

    return c, r

# For printing https://stackoverflow.com/questions/17489353/printing-a-mixed-type-dictionary-with-format-in-python
class my_dict(dict):                                              
    def __str__(self):
        return str({k:round(v,2) if isinstance(v,float) else v  for k,v in self.items()})

class WingPenaltyOptimizer:
    """Works with dirtran form of x only"""
    WRT_TRAJ = 0
    WRT_PARAMS = 1
    WRT_TRAJ_PARAMS = 2

    # Options for the solver
    GRADIENT_DESCENT = 0
    NEWTON_METHOD = 1
    GAUSS_NEWTON = 2

    def __init__(self, N):
        self.N = N
        self._Nx = (self.N+1) * m.nx + self.N*m.nu + 1 #dirtran size + timestep h
    
    def update(self, traj0, params0, mode=WRT_TRAJ, opt={}, Niter=1):
        # Some error checking
        assert len(traj0) == self._Nx
        assert len(params0) == len(params)

        HESS_REG = opt.get('hessreg', 1e-3)
        method = opt.get('method', self.GAUSS_NEWTON)
        optnp = dict(opt, **{'pen':False})

        if mode == self.WRT_PARAMS:
            Jtup = lambda p : Jcosttraj_penalty(traj0, self.N, p, opt)
            x0 = params0
        elif mode == self.WRT_TRAJ:
            Jtup = lambda traj : Jcosttraj_penalty(traj, self.N, params0, opt)
            x0 = traj0
        elif mode == self.WRT_TRAJ_PARAMS:
            Jtup = lambda trajp : Jcosttraj_penalty(trajp[:self._Nx], self.N, trajp[self._Nx:], opt)
            x0 = np.hstack((traj0, params0))
        else:
            raise ValueError('Invalid mode')
        
        # separately get the non-quadratic and quadratic terms
        Jnq = lambda x : Jtup(x)[0]
        r = lambda x : Jtup(x)[1]

        if method != self.GAUSS_NEWTON:
            def J(x):
                # handle all costs the same
                rr = r(x)
                return Jnq(x) + rr.T @ rr
        else:
            J = Jnq
            # Approximate the gradient, Hessian for these terms with Jr
            Jr = jacobian(r)
            # # Composing gradients of ri(xi), where xi = (yi,y{i+1},ui)
            # # NOTE: leaving this here, but it seems this is actually slower
            # N, nx, nu = self.N, m.nx, m.nu
            # def ri(xi):
            #     # r0(y0,y1,u0)
            #     x = np.hstack((xi[:2*nx], np.zeros((N-1)*nx), xi[-nu:], np.zeros((N-1)*nu)))
            #     return r(x)[:nx]
            # Jri = jacobian(ri)
            # # Now assemble Jr from these components
            # szr = (N)*nx
            # Jr0 = np.zeros((szr, self._Nx))
            # for i in range(N-1):
            #     xi = np.hstack((x0[nx*i:nx*(i+2)], x0[(N+1)*nx+nu*i : (N+1)*nx+nu*(i+1)]))
            #     Jr0i = Jri(xi)
            #     # Now put in the correct location
            #     Jr0[nx*i : nx*(i+1), nx*i:nx*(i+2)] = Jr0i[:,:2*nx]
            #     Jr0[nx*i : nx*(i+1), (N+1)*nx + nu*i : (N+1)*nx + nu*(i+1)] = Jr0i[:,-nu:]
                
        DJ = jacobian(J)
        D2J = hessian(J)

        prof = my_dict({'e': 0, 'eJ': 0, 'eH': 0, 'ls': 0})
        
        for minorIter in range(Niter):
            t0 = time.perf_counter()
            J0 = J(x0)
            if method == self.GAUSS_NEWTON:
                r0 = r(x0)
            prof['e'] = time.perf_counter() - t0
            t0 = time.perf_counter()
            DJ0 = DJ(x0)
            if method == self.GAUSS_NEWTON:
                Jr0 = Jr(x0)
                DJ0 += 2 * Jr0.T @ r0
            prof['eJ'] = time.perf_counter() - t0
            t0 = time.perf_counter()
            D2J0 = D2J(x0)
            if method == self.GAUSS_NEWTON:
                D2J0 += 2 * Jr0.T @ Jr0
            prof['eH'] = time.perf_counter() - t0

            # descent direction
            if method == self.GRADIENT_DESCENT:
                v = -DJ0 # gradient descent
            elif method in [self.NEWTON_METHOD, self.GAUSS_NEWTON]:
                # Newton's method http://www.stat.cmu.edu/~ryantibs/convexopt-S15/lectures/14-newton.pdf
                try:
                    # regularization
                    D2J0 += HESS_REG * np.eye(D2J0.shape[0])
                    v = -np.linalg.inv(D2J0) @ DJ0
                except np.linalg.LinAlgError:
                    # last u (last diag elem) is 0 - makes sense
                    print(np.linalg.eigvals(D2J0))
                    raise

            t0 = time.perf_counter()
                    
            # backtracking line search
            alpha = 0.4
            beta = 0.9
            s = 1
            x1 = x0 + s * v
            J1 = J(x1)
            if J1 > J0 and mode == self.WRT_PARAMS:
                # FIXME: why is the direction backwards sometimes in the WRT_PARAMS mode??
                v = -v
                x1 = x0 + s * v
                J1 = J(x1)
            # search for s
            while J1 > J0 + alpha * s * DJ0.T @ v and s > 1e-6:
                s = beta * s
                x1 = x0 + s * v
                J1 = J(x1)

            prof['ls'] = time.perf_counter() - t0

            # debugging
            print(mode, minorIter, prof, "cost {:.1f} -> {:.1f}".format(J0, J1), end = " ")
            if mode == self.WRT_TRAJ:
                print("h {:.2f}ms -> {:.2f}ms".format(1e3*x0[-1], 1e3*x1[-1]), end = " ")
                assert x1[-1] > 0, "Negative timestep"
            print()

            # for next iteration
            x0 = x1

        # perform Newton update
        return x1, J, J1
        
    def plotTrajs(self, *args):
        """Helper function to plot a bunch of trajectories superimposed"""
        umin, umax, xmin, xmax = m.limits
        trajt = range(self.N) # timestep is the last elem
        yend = (self.N) * m.nx # N to ignore the last one
        ustart = (self.N+1) * m.nx
        uend = ustart + self.N*m.nu
        _, ax = plt.subplots(3)
        for arg in args:
            ax[0].plot(trajt * arg[-1], arg[0:yend:m.nx], '.-')
        for yy in [xmin[0], xmax[0], 0]:
            ax[0].axhline(y=yy, color='k', alpha=0.3)
        ax[0].set_ylabel('act. disp (m)')
        for arg in args:
            ax[1].plot(trajt * arg[-1], arg[1:yend:m.nx], '.-')
        for yy in [xmin[1], xmax[1], np.pi/4, -np.pi/4]:
            ax[1].axhline(y=yy, color='k', alpha=0.3)
        ax[1].set_ylabel('hinge angle (rad)')
        for arg in args:
            ax[2].plot(trajt * arg[-1], arg[ustart:uend:m.nu], '.-')
        ax[2].axhline(y=umin[0], color='k', alpha=0.3)
        ax[2].axhline(y=umax[0], color='k', alpha=0.3)
        ax[2].set_ylabel('act. force (N)')

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

# Plot stuff wrt params ---

def plotTrajWrtParams(p0s, p1s, traj0, N, dpen=1e3, paramsPath=None):
    """Debug convexity wrt params"""
    P0S, P1S = np.meshgrid(p0s, p1s)
    JS = np.zeros_like(P0S)
    fig, ax = plt.subplots(2)

    def Jp(p, _dpen):
        c, r = Jcosttraj_penalty(traj0, N, p, opt={'dynamics':_dpen, 'periodic':0, 'input':1e4, 'state': 1e0})
        return c + r.T @ r
        
    for ii in range(JS.shape[0]):
        for jj in range(JS.shape[1]):
            JS[ii,jj] = Jp([P0S[ii,jj], P1S[ii,jj]], _dpen=dpen)

    ax[0].plot(p0s, [Jp([cc, np.mean(p1s)], dpen) for cc in p0s], '.-')
    # ax[0].plot(p0s, [Jp([cc, np.mean(p1s)], 1e4) for cc in p0s], '.-', label='1e4')
    ax[0].set_ylabel("ci")

    ax[1].plot(p1s, [Jp([np.mean(p0s), ti], dpen) for ti in p1s], '.-')
    # ax[1].plot(p1s, [Jp([np.mean(p0s), ti], 1e4) for ti in p1s], '.-', label='1e4')
    # ax[1].legend()
    ax[1].set_ylabel("Ji")
    ax[1].set_xlabel("T")

    # from mpl_toolkits.mplot3d import axes3d
    fig = plt.figure()
    ax = fig.add_subplot(111)#, projection='3d')
    # ax.plot_surface(P0S, P1S, JS, cmap=plt.get_cmap('gist_earth'))
    ax.contourf(P0S, P1S, JS, 50, cmap=plt.get_cmap('gist_earth'))
    if paramsPath is not None:
        ax.plot([pi[0] for pi in paramsPath], [pi[1] for pi in paramsPath], 'r*-')
    ax.set_xlabel('cbar')
    ax.set_ylabel('T')

