
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
from controlutils.py.model import Model

class Wing2DOF(Model):
    """
    NOTE: for optimizing transmission ratio
    Thinking of y = (sigma_actuator, psi, dsigma_actuator, dpsi)
    u = (tau_actuator)
    sigma = sigma_actuator * T; tau = tau_actuator / T
    """
    ny = 4
    nu = 1
    # CONST
    R = 20.0
    γ = 5.0 # wing shape fitting

    def aero(self, y, u, _params=[]):
        CLmax = 1.8
        CDmax = 3.4
        CD0 = 0.4
        ρ = 1.225e-3 # [mg/(mm^3)]
        R = 15e-3
        
        # unpack
        cbar, T = _params
        σ, Ψ, σ̇, Ψ̇ = np.array([T, 1, T, 1]) * y # [mm, rad, mm/ms, rad/ms]
        cΨ = np.cos(Ψ)
        sΨ = np.sin(Ψ)
        α = np.pi/2 - Ψ # AoA
        
        """    
        Use this in Mathematica to debug the signs.
        Manipulate[
        \[Alpha] = -\[Pi]/2 + \[Psi];
        Graphics@
        Arrow[{{0, 0},
            {((CDmax + CD0)/2 - (CDmax - CD0)/2*Cos[2 \[Alpha]]), CLmax Sin[2 \[Alpha]]} Sign[-d\[Sigma]]}],
        {d\[Sigma], -1, 1}, {\[Psi], -\[Pi]/2, \[Pi]/2}]
        """
        # aero force
        wing1 = np.array([σ, 0])
        rcopnondim = 0.25 + 0.25 / (1 + exp(5.0*(1.0 - 4*(π/2 - np.abs(Ψ))/π))) # [(6), Chen (IROS2016)]
        paero = wing1 + np.array([[cΨ, -sΨ], [sΨ, cΨ]]) @ np.array([0, -rcopnondim*cbar])
        Jaero = np.array([[1, rcopnondim*cbar * cΨ], [0, rcopnondim*cbar * sΨ]])
        Caero = np.array([(CDmax+CD0)/2 - (CDmax-CD0)/2 * np.cos(2*α), CLmax * np.sin(2*α)])
        Faero = 1/2 * ρ * cbar * self.R * σ̇**2 * Caero * np.sign(-σ̇) #[mN]

        return paero, Jaero, Faero

    def dydt(self, y, u, _params=[]):
        """Continuous dynamics second order model"""
        cbar, T = _params
        σ, Ψ, σ̇, Ψ̇ = np.array([T, 1, T, 1]) * y # [mm, rad, mm/ms, rad/ms]
        cΨ = np.cos(Ψ)
        sΨ = np.sin(Ψ)

        # params
        mspar = 0 # [mg]
        mwing = 0.51 # [mg]
        Iwing = mwing * cbar^2 # cbar is in mm
        kσ = 0 # [mN/mm]
        bσ = 0 # [mN/(mm/ms)]
        kΨ = 5 # [mN-mm/rad]
        bΨ = 3 # [mN-mm/(rad/ms)]

        # inertial terms
        M = np.array([[mspar+mwing, cbar*mwing*cΨ], [cbar*mwing*cΨ, Iwing + cbar**2*mwing]])
        corgrav = np.array([kσ*σ - cbar*mwing*sΨ*Ψ̇**2, kΨ*Ψ])
        # non-lagrangian terms
        τdamp = np.array([-bσ * σ̇, -bΨ * Ψ̇])
        _, Jaero, Faero = self.aero(y, u, _params)
        τaero = Jaero.T @ Faero
        # input
        τinp = np.array([u[0] / T, 0])

        ddq = np.linalg.inv(M) @ (-corgrav + τdamp + τaero + τinp)

        return np.array([y[2], y[3], ddq[0] / T, ddq[1]])

    @property
    def limits(self):
        # This is based on observing the OL trajectory
        umax = np.array([75.0]) # [mN]
        umin = -umax
        xmax = np.array([300e-3, 1.5, 100, 100]) # [mm, rad, mm/ms, rad/ms]
        xmin = -xmax
        return umin, umax, xmin, xmax

    def createInitialTraj(self, opt, N, freq, posGains, params):
        σmax = self.limits[-1][0]
        def strokePosController(y, t):
            σdes = 0.9 * σmax * np.sin(freq * 2 * np.pi * t)
            return posGains[0] * (σdes - y[0]) - posGains[1] * y[2]
        
        strokePosControlVF = lambda t, y: self.dydt(y, [strokePosController(y, t)], params)
        # OL traj1
        teval = np.arange(0, 100, 0.1) # [ms]
        # Sim of an openloop controller
        sol = solve_ivp(strokePosControlVF, [teval[0], teval[-1]], [0.,1.,0.,0.], dense_output=True, t_eval=teval)

        # # Animate whole traj
        # Nt = length(sol.t)
        # @gif for k = 1:Nt
        #     yk = sol.u[k]
        #     uk = [strokePosController(yk, sol.t[k])]
        #     drawFrame(m, yk, uk, params)
        # end
        # # Plot
        # σt = plot(sol, vars=3, ylabel="act vel [m/s]")
        # Ψt = plot(sol, vars=2, ylabel="hinge ang [r]")
        # plot(σt, Ψt, layout=(2,1))
        # gui()

        starti = 170
        traj_s = np.s_[starti:starti + 2*(N+1):2]
        trajt = sol.t[traj_s]
        trajym = sol.y[:, traj_s] # ny x (N+1) array
        δt = trajt[1] - trajt[0]
        # Transcribe to x = [y1,...,yNp1, u1,...,uNp1]
        traj = np.hstack((
            trajym.reshape((N+1)*self.ny, order='F'), 
            [strokePosController(trajym[:,i], trajt[i]) for i in range(N+1)]
            ))

        if opt['vart']:
            traj = np.append(traj, δt)
        else:
            print("Initial traj δt=", δt, ", opt.fixedδt=", opt['fixedδt'])

        return trajt - trajt[0], traj

    def plotTrajs(self, opt, params, *args):
        """Helper function to plot a bunch of trajectories superimposed"""
        umin, umax, xmin, xmax = self.limits
        N, δt, yk, uk = self.modelInfo(opt, args[0])
        
        cbar, T = params

        trajt = range(N+1) * δt
        yend = (N+1)*self.ny
        uend = (N+1)*(self.ny + self.nu)
        
        _, ax = plt.subplots(2,2)
        ax = np.ravel(ax)
        for arg in args:
            ax[0].plot(trajt, arg[0:yend:self.ny] * T / (self.R/2), '.-')
        for yy in [xmin[0], xmax[0], 0]:
            ax[0].axhline(y=yy, color='k', alpha=0.3)
        ax[0].set_ylabel('stroke ang [r]')

        for arg in args:
            ax[1].plot(trajt, arg[1:yend:self.ny], '.-')
        for yy in [xmin[1], xmax[1], np.pi/4, -np.pi/4]:
            ax[1].axhline(y=yy, color='k', alpha=0.3)
        ax[1].set_ylabel('hinge ang [r]')

        for arg in args:
            ax[2].plot(trajt, arg[yend:uend:self.nu], '.-')
        ax[2].axhline(y=umin[0], color='k', alpha=0.3)
        ax[2].axhline(y=umax[0], color='k', alpha=0.3)
        ax[2].set_ylabel('stroke force [mN]')

        def aeroPlotVec(_traj):
            _, _, _yk, _uk = self.modelInfo(opt, _traj)
            Faerok = lambda k : self.aero(_yk(k), _uk(k), params)[-1]
            return np.vstack([Faerok(k) for k in range(N+1)])
        
        for arg in args:
            ax[3].plot(trajt, aeroPlotVec(arg)[:,1], '.-')
        ax[3].axhline(y=0, color='k', alpha=0.3)
        ax[3].set_ylabel('lift [mN]')
        
        plt.tight_layout()

    def eulerIntegrate(self, opt, traj, params):
        """Test applying euler integration to the initial traj"""
        trajei = traj.copy()
        N, δt, yk, uk = self.modelInfo(opt, trajei)

        for k in range(N):
            trajei[(k+1)*self.ny : (k+2)*self.ny] = yk(k) + δt * self.dydt(yk(k), uk(k), params)
        return trajei

