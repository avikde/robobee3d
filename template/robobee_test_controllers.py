"""Each controller has a set of parameters and ranges as well as a function to evaluate it"""
import pybullet as p
import autograd.numpy as np
from scipy.spatial.transform import Rotation
from ca6dynamics import dynamicsTerms
from wlqppy import WLController
import valuefunc
from template_controllers import createMPC, reactiveController
from genqp import quadrotorNLVF, Ib
import flight_tasks

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

class WaypointHover(RobobeeController):
    """Simplified version of Pakpong (2013). u = [Vmean, uoffs, udiff, h2]"""
    def __init__(self, wrenchMapPoptsFile, initialPos, constdqdes=None, useh2=False, useMPC=True, useWLQP=True, task='helix'):
        super(WaypointHover, self).__init__({'freq': (0, 0.3, 0.16)})
        self.wf = WaveformGenerator()
        self.posdes = np.array([0.,0.,100.])
        self.constdqdes = constdqdes
        self.useh2 = useh2
        self.useMPC = useMPC
        self.initialPos = np.asarray(initialPos)
        self.task = task
        self.useWLQP = useWLQP

        # NOTE: This is not actually used: need to put in wlcontroller.cpp or wlcontroller.c
        popts = np.load(wrenchMapPoptsFile)

        self.wl = WLController(np.ravel(popts), 1000)
        self.u4 = [140.0,0.,0.,0.]

        # For acc reference. only need to get once for now for hover task
        # In this R (for quadrotors) weigh the pitch and yaw torques high
        self.S = valuefunc.quadrotorS(9.81e-3, Qpos=[0,0,10,0.1,0.1,0.1], Qvel=[1,1,10,0.1,0.1,0.1], Rdiag=[1,1,1,1])
        self.printCtr = 0

        # upright MPC
        mpcopts = {'ws':1.5, 'wds':1e3, 'wpr':2e-2, 'wvr':2e1, 'wpf':4e-2, 'wvf':2e1, 'TtoWmax':3}#line
        self.up, _ = createMPC(**mpcopts, popts=np.ravel(popts))

    def templateVF(self, t, p, dp, s, ds, posdes, dposdes, kpos=[0.5e-3,5e-1], kz=[1e-3,2e-1], ks=[4e-3,0.3e0]):
        # TEST
        self.printCtr = (self.printCtr + 1) % 100
        # If we want to go to a certain position have to set sdes
        sdes = np.zeros(3)
        pos2err = p[0:2] - posdes[0:2]
        dpos2err = dp[0:2] - dposdes[:2]
        sdes[0:2] = -kpos[0] * pos2err - kpos[1] * dpos2err
        # if self.printCtr == 0:
        #     print(self.posdes[:2], pos2err, self.pos2errI, sdes[:2],  -1e-1 * pos2err, - 1e-1 * dp[0:2],  - 1e0 * self.pos2errI)
        sdes[0:2] = np.clip(sdes[0:2], -0.5 * np.ones(2), 0.5 * np.ones(2))
        fTorn = ks[0] * (s - sdes) + ks[1] * ds
        fTorn[2] = 0 # z element

        # for position z
        fTpos = kz[0] * (posdes - p) - kz[1] * dp
        fTpos[:2] = np.array([0,0])

        return fTpos, fTorn
    
    def accReference(self, t, q0, dq0):
        """Used in the C version; returns accdes"""
        # # Simple quadratic VF on vel kp * ||dq0 - dqdes||^2 OLD WORKS
        # kdq = np.array([0,0,1,0.1,0.1,0.1])
        # return kdq * (dqdes - dq0)

        # New try template VF
        e3h = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
        Rb = Rotation.from_quat(q0[3:]).as_matrix()
        p = np.asarray(q0[:3])
        omega = dq0[3:]
        dp = dq0[:3]
        s = Rb @ np.array([0,0,1])
        ds = -Rb @ e3h @ omega

        if self.task == 'helix':
            self.posdes, dpdes, sdes = flight_tasks.helix(t, self.initialPos)
        elif self.task == 'line':
            self.posdes, dpdes, sdes = flight_tasks.straightAcc(t, self.initialPos, vdes=2, tduration=750)
        elif self.task == 'flip':
            self.posdes, dpdes, sdes = flight_tasks.flip(t, self.initialPos)

        if self.useMPC:
            # Upright MPC
            uquad, ddqdes, uwlqp = self.up.update(p, Rb, dq0, self.posdes, dpdes, sdes)
            # ddqdes[:3] = Rb.T @ ddqdes[:3] # Convert to body frame?
            # ddqdes[3:] = Rb.T @ ddqdes[3:] # Convert to body frame?
            return ddqdes
        else:
            # Template controller <- LATEST
            fTpos, fTorn = self.templateVF(t, p, dp, s, ds, self.posdes, dpdes)
            fAorn = -e3h @ Rb.T @ fTorn
            return np.hstack((fTpos, fAorn))

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

    def accController(self, *args):
        t, qb, dqb = args

        M0, h0 = dynamicsTerms(qb, dqb)
        self.accdes = self.accReference(t, qb, dqb)

        if self.useWLQP:
            self.u4, w0 = self.wl.update(self.u4, h0, M0 @ self.accdes)
            self.up.T0 += 1 * (w0[2]/M0[2,2] - self.up.T0)
        else:
            # Manual mapping
            self.u4 = self.manualMapping(self.accdes)
            # self.u4 = uwlqp

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
        # control
        return self.accController(t, qb, dqb)
        
    def manualMapping(self, accdes, kx=1e2, ky=1e3, kV1=1e3, kV0=105):
        """Low level mapping to torques. Can be replaced"""
        Vmean = kV0 + (accdes[2] - 9.81e-3) * kV1
        uoffs = np.clip(ky * accdes[4], -0.5, 0.5)
        udiff = np.clip(kx * accdes[3], -0.3, 0.3)
        # print(accdes)
        return Vmean, uoffs, udiff, 0
