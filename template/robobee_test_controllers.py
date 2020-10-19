"""Each controller has a set of parameters and ranges as well as a function to evaluate it"""
import pybullet as p
import autograd.numpy as np
from scipy.spatial.transform import Rotation
from ca6dynamics import dynamicsTerms
from wlqppy import WLController
import valuefunc
from uprightmpc2 import UprightMPC2

class WaveformGenerator(object):
    """A simple waveform generator that keeps track of phase"""
    def __init__(self):
        self.ph = 0 # phase
        self.tprev = 0
    
    def step(self, t, f):
        self.ph += 2 * np.pi * f * (t - self.tprev)
        self.tprev = t
    
    def get(self, h2=0., h3=0.1):
        return (1 + h3) * np.sin(self.ph) + h3 * np.sin(3 * self.ph) + h2 * np.sin(2 * self.ph)

class RobobeeController(object):
    def __init__(self, params):
        # Params stored as (min, max, default) tuples
        self.params = params
        # Params using pybullet GUI (sliders)
        self.dbgIDs = {k : p.addUserDebugParameter(k, *params[k]) for k in params.keys()}
        self.pdes = np.zeros(6) # controller should set
        self.u4 = np.zeros(4) # for logging

    def P(self, k):
        try:
            return p.readUserDebugParameter(self.dbgIDs[k])
        except: # if in p.DIRECT mode, just return the default
            return self.params[k][-1]
    
    def update(self, t, q, dq):
        raise NotImplementedError

class OpenLoop(RobobeeController):
    """Open-loop controller"""
    def __init__(self):
        super(OpenLoop, self).__init__({'freq': (0, 0.3, 0.16), 'umean': (0, 200, 150), 'udiff': (-0.5, 0.5, 0), 'uoffs': (-0.5, 0.5, 0), 'testFL': (-10,10,0), 'testFR': (-10,10,0)})
        self.wf = WaveformGenerator()

    def update(self, t, q, dq):
        # force control
        umean = self.P('umean')
        udiff = self.P('udiff')
        uoffs = self.P('uoffs')
        w = self.wf.update(t, self.P('freq'))
        return np.array([1 + udiff, 1 - udiff]) * umean * (w + uoffs)

def positionControllerPakpongLike(posdes, qb, dqb):
    """Unified interface for a controller that spits out desired momentum"""
    Rb = Rotation.from_quat(qb[3:7])
    omega = dqb[3:6]

    Rm = Rb.as_matrix()
    zdes = np.array([0.,0.,1.]) # desired z vector
    # upright controller
    zdes[0:2] = np.clip(0.01 * (posdes[0:2] - qb[0:2]) - 1 * dqb[0:2], -0.5 * np.ones(2), 0.5 * np.ones(2))
    # zdes /= np.linalg.norm(zdes)

    ornError = np.array([
        [Rm[0,1], Rm[1,1], Rm[2,1]], 
        [-Rm[0,0], -Rm[1,0], -Rm[2,0]], 
        [0,0,0]]) @ zdes # Pakpong (2013) (6)
    Iomegades = -1e2*ornError - 1e3*omega
    
    return np.hstack((0, 0, 0.01 * (posdes[2] - qb[2]), Iomegades))

class WaypointHover(RobobeeController):
    """Simplified version of Pakpong (2013). u = [Vmean, uoffs, udiff, h2]"""
    def __init__(self, wrenchMapPoptsFile, constPdes=None, useh2=True):
        super(WaypointHover, self).__init__({'freq': (0, 0.3, 0.16)})
        self.wf = WaveformGenerator()
        self.posdes = np.array([0.,0.,100.])
        self.constPdes = constPdes
        self.useh2 = useh2

        # NOTE: This is not actually used: need to put in wlcontroller.cpp or wlcontroller.c
        popts = np.load(wrenchMapPoptsFile)

        self.wl = WLController(np.ravel(popts), 1000)
        self.u4 = [140.0,0.,0.,0.]

        # self.momentumController = self.manualMapping
        self.momentumController = self.wrenchLinWrapper
        self.positionController = positionControllerPakpongLike
        # For momentum reference. only need to get once for now for hover task
        # In this R (for quadrotors) weigh the pitch and yaw torques high
        self.S = valuefunc.quadrotorS(9.81e-3, Qpos=[0,0,10,0.1,0.1,0.1], Qvel=[1,1,10,0.1,0.1,0.1], Rdiag=[1,1,1,1])
        self.printCtr = 0

        # upright MPC
        dt = 5
        N = 3
        g = 9.81e-3
        ws = 1e1
        wds = 1e3
        wpr = 1
        wpf = 5
        wvr = 1e3
        wvf = 2e3
        wthrust = 1e-1
        wmom = 1e-2
        smin = np.array([-2,-2,0.5])
        smax = np.array([2,2,1.5])
        self.up = UprightMPC2(N, dt, g, smin, smax, ws, wds, wpr, wpf, wvr, wvf, wthrust, wmom)

    def templateVF(self, t, p, dp, s, ds):
        # TEST
        trajomg = 5e-3
        self.posdes = np.array([50 * np.sin(trajomg * t),0,100])
        dposdes = np.array([50 * trajomg * np.cos(trajomg * t),0,0])

        self.printCtr = (self.printCtr + 1) % 100
        # If we want to go to a certain position have to set sdes
        sdes = np.zeros(3)
        pos2err = p[0:2] - self.posdes[0:2]
        dpos2err = dp[0:2] - dposdes[:2]
        sdes[0:2] = -1e-4 * pos2err - 1e0 * dpos2err
        # if self.printCtr == 0:
        #     print(self.posdes[:2], pos2err, self.pos2errI, sdes[:2],  -1e-1 * pos2err, - 1e-1 * dp[0:2],  - 1e0 * self.pos2errI)
        sdes[0:2] = np.clip(sdes[0:2], -0.5 * np.ones(2), 0.5 * np.ones(2))
        fTorn = 10 * (s - sdes) + 1e3 * ds
        fTorn[2] = 0 # z element

        # for position z
        fTpos = 1e-1 * (self.posdes - p) - 1e1 * dp
        fTpos[:2] = np.array([0,0])

        return fTpos, fTorn
    
    def momentumReference(self, t, q0, dq0, M0, pdes):
        """Used in the C version; returns pdotdes"""
        # # Simple quadratic VF on momentum kpmom * ||p0 - pdes||^2 OLD WORKS
        # kpmom = np.array([0,0,1,0.1,0.1,0.1])
        # return kpmom * (pdes - M0 @ dq0)

        # New try template VF
        e3h = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
        Rb = Rotation.from_quat(q0[3:]).as_matrix()
        omega = dq0[3:]
        p = q0[:3]
        dp = dq0[:3]
        s = Rb @ np.array([0,0,1])
        ds = -Rb @ e3h @ omega

        # Upright MPC
        self.posdes = np.array([0,0,150])
        dpdes = np.zeros(3)
        ddqdes = self.up.updateGetAccdes(p, Rb, dq0, self.posdes, dpdes)
        ddqdes[:3] = Rb.T @ ddqdes[:3] # Convert to body frame?
        pdotdes = M0 @ ddqdes
        # print(pdotdes)
        return pdotdes

        # # Template controller <- LATEST
        # fTpos, fTorn = self.templateVF(t, p, dp, s, ds)
        # fAorn = -e3h @ Rb.T @ fTorn
        # return np.hstack((fTpos, fAorn))

        # # Here the u is Thrust,torques (quadrotor template)
        # pT = q0[:3]
        # Rb = Rotation.from_quat(q0[3:])
        # phiT = Rb.as_euler('xyz')
        # # https://github.com/avikde/robobee3d/pull/178
        # xA = np.hstack((pT - np.array([0,0,100]), phiT, dq0))
        # Dpi = np.block([[Rb.as_dcm(), np.zeros((3,3))], [np.zeros((3,3)), Rb.as_dcm()]]) @ np.linalg.inv(M0)
        # pddes = -Dpi.T @ self.S[6:,:] @ xA
        # # # Rotate
        # # bRw = Rotation.from_euler('z', phi0[2])
        # # pddes = np.hstack((bRw.apply(pddes[:3]), bRw.apply(pddes[3:])))
        # print(pddes)
        # return pddes

    def wrenchLinWrapper(self, *args):
        t, qb, dqb, pdes = args

        M0, h0 = dynamicsTerms(qb, dqb)
        pdotdes = self.momentumReference(t, qb, dqb, M0, pdes)
        self.u4 = self.wl.update(self.u4, h0, pdotdes)

        Vmean, uoffs, udiff, h2 = self.u4
        # # test
        # Vmean = 110# + 1000 * (pdes[2] - dqb[2])
        # udiff = 0
        # uoffs = 0
        if not self.useh2:
            h2 = 0
        
        self.wf.step(t, self.P('freq'))
        wL = self.wf.get(h2=h2)
        wR = self.wf.get(h2=-h2)
        udiff = np.clip(udiff, -0.3, 0.3)
        uoffs = np.clip(uoffs, -0.3, 0.3)
        return Vmean * np.array([(1 + udiff) * (wL + uoffs), (1 - udiff) * (wR + uoffs)])

    def update(self, t, q, dq):
        # unpack
        qb = q[-7:]
        dqb = dq[-6:]
        # momentum-based control
        self.pdes = self.constPdes if self.constPdes is not None else self.positionController(self.posdes, qb, dqb)
        return self.momentumController(t, qb, dqb, self.pdes)
        
    def manualMapping(self, t, qb, dqb, pdes):
        """Low level mapping to torques. Can be replaced"""
        self.wf.step(t, self.P('freq'))
        w = self.wf.get(t, self.P('freq')) # no h2 usage in this
        umean = 120 + 1000 * (pdes[2] - dqb[2])
        udiff = np.clip(0.01*pdes[3], -0.5, 0.5)
        uoffs = np.clip(0.1*pdes[4], -0.5, 0.5)
        self.u4 = np.array([umean, uoffs, udiff, 0.]) # last is h2
        return np.array([1 + udiff, 1 - udiff]) * umean * (w + uoffs)
