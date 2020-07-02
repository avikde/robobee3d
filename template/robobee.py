
import numpy as np
import pybullet as p
import pybullet_data
from time import monotonic
from scipy.spatial.transform import Rotation # TODO: eliminate

# Good parameters
CD0 = 0.4
CDmax = 3.4
CLmax = 1.8
AERO_REGULARIZE_EPS = 1e-10 # stops undefined AoA when no wind
# True in Chen (2017) science robotics, but differently calculated in Osborne (1951)
RHO = 1.225e-3 # [mg/(mm^3)]

rcopnondim = 0.5

def aerodynamics(theta, dtheta, lrSign, params):
    """Return aerodynamic force and instantaneous CoP in the body frame (both in R^3). If flip=False it will work for the left wing (along +y axis), and if flip=True it will """
    theta[0] *= lrSign
    dtheta[0] *= lrSign

    # Faero = 1/2 * ρ * dφ^2 * Aw^2 * AR * m.r2h^2 * (Caero[1]*eD + Caero[2]*eL*sign(-dφ)) # [mN]
    Aw = 54.4
    r1h = 0.49
    r2h = 0.551
    cbar = 3.2 #params['cbar']
    rcnd = rcopnondim

    cbar2 = cbar**2
    AR = Aw / cbar2
    Lw = Aw/cbar
    ycp = params['ycp']#params['Rwing']*r1h #

    s, c = np.sin(theta[0]), np.cos(theta[0])
    sh, ch = np.sin(theta[1]), np.cos(theta[1])
    pcopB = np.array([0,params['Roffs'],params['d']]) + \
        np.array([-ycp*s - cbar*rcnd*c*sh, ycp*c - cbar*rcnd*s*sh, -cbar*rcnd*ch])
    aoa = 0.5 * np.pi - theta[1]
    eD = np.array([dtheta[0]*c, dtheta[0]*s, 0])/(np.abs(dtheta[0]) + 1e-4) # so it is well-defined
    eL = np.array([0, 0, 1])
    Caero = [((CDmax + CD0)/2 - (CDmax - CD0)/2 * np.cos(2*aoa)), CLmax * np.sin(2*aoa)]
    FaeroB = 0.5 * RHO * dtheta[0]**2 * Aw**2 * AR * r2h**2 * (Caero[0]*eD + Caero[1]*eL*np.sign(-dtheta[0]))
    
    if lrSign < 0:
        yflip = np.diag([1,-1,1])
        FaeroB = yflip @ FaeroB
        pcopB = yflip @ pcopB
    return FaeroB, pcopB

def wrench(F, p):
    return np.hstack((F, np.cross(p, F)))

def actuatorModel(V, qact, dqact):
    """Return actuator force for applied voltage input and current actuator state"""
    T = 2.6666
    return 1 / T * (50./180. * V) # proportional model FIXME: why low

class WingFlapFilter(object):
    """Utilities to filter based on wing flap frequency"""
    def __init__(self, freq, timestep, n):
        self.buf = np.zeros((int(1/(freq * timestep)), n)) # each row will be a sample
        self.bufi = 0
        self.firstTime = True
    
    def update(self, x):
        """Operators on a vector x"""
        if self.firstTime:
            self.firstTime = False
            for i in range(self.buf.shape[0]):
                self.buf[i,:] = x
        else:
            self.buf[self.bufi,:] = x
        self.bufi = (self.bufi + 1) % len(self.buf)
        return np.mean(self.buf, axis=0)

class RobobeeSim():
    """Robobee simulator using pybullet"""

    # Parameters
    FAERO_DRAW_SCALE = 2.0
    pcomLastDraw = np.zeros(3)

    def __init__(self, connMode, camLock=True, slowDown=True, timestep=1.0, gui=0, filtfreq=0.15):
        self._camLock = camLock
        self._slowDown = slowDown
        self.TIMESTEP = timestep
        # Init sim
        physicsClient = p.connect(connMode)
        # Set up the visualizer
        p.configureDebugVisualizer(p.COV_ENABLE_WIREFRAME, 0)
        p.configureDebugVisualizer(p.COV_ENABLE_SHADOWS, 0)
        p.configureDebugVisualizer(p.COV_ENABLE_GUI, gui)
        # p.configureDebugVisualizer(p.COV_ENABLE_MOUSE_PICKING, 0)
        p.setRealTimeSimulation(0)
        p.setTimeStep(self.TIMESTEP)
        p.setGravity(0,0,-9.81e-3)
        
        # load background
        p.setAdditionalSearchPath(pybullet_data.getDataPath()) #optionally
        self.planeId = p.loadURDF("plane.urdf", globalScaling=100.0)
        self.reset(filtfreq)
    
    def reset(self, filtfreq):
        """Reset states, to restart the sim, for instance"""
        self.simt = 0
        # vectors for storing states
        self.q = np.zeros(11)
        self.dq = np.zeros(10)
        self.tLastPrint = 0
        self.tWallLastDraw = monotonic()
        self.wf = WingFlapFilter(filtfreq, self.TIMESTEP, 6) # filter for body vel

    def wTb(self, FaeroB, pcopB):
        # Body to world frame --
        pcom = self.q[4:7]
        Rb = Rotation.from_quat(self.q[7:11]) # scalar-last format
        pcopW = pcom + Rb.apply(pcopB)
        FaeroW = Rb.apply(FaeroB)
        # for external torque about wing hinge, use r X F
        # hingeTorque = np.cross(0.5 * self.cbar * Rb.apply(chordB), FaeroW)
        return FaeroW, pcopW#, hingeTorque


    def load(self, filename, basePosition, *args, **kwargs):
        '''startPos = xyz
        startOrientation = quaternion xyzw
        '''
        self.bid = p.loadURDF(filename, basePosition, *args, **kwargs, flags=p.URDF_USE_INERTIA_FROM_FILE)
        # See https://github.com/bulletphysics/bullet3/issues/2152
        p.changeDynamics(self.bid, -1, maxJointVelocity=10000)
        p.changeDynamics(self.bid, -1, linearDamping=0,	angularDamping=0)
        
        self.Nj = p.getNumJoints(self.bid)
        self.pcomLastDraw = basePosition
        
        self.jointId, self.urdfParams = self.getInfoFromURDF()
        
        # Get the joints to reasonable positions
        self.resetJoints(self.bid, range(4), np.zeros(4), np.zeros(4))
        
        # Passive hinge dynamics implemented as position control rather than joint dynamics
        p.setJointMotorControlArray(self.bid, [1,3], p.PD_CONTROL, targetPositions=[0,0], positionGains=self.urdfParams['khinge']*np.ones(2), velocityGains=self.urdfParams['bhinge']*np.ones(2))

        return self.bid

    def getInfoFromURDF(self):
        '''Get info from the URDF
        '''
        urdfParams = {}

        #  Since each link is connected to a parent with a single joint,
        # the number of joints is equal to the number of links. Regular links have link indices in the range
        # [0..getNumJoints()] Since the base is not a regular 'link', we use the convention of -1 as its link
        # index
        jointId = {}
        for j in range(self.Nj):
            jinfo = p.getJointInfo(self.bid, j)
            # unpack
            jointIndex, jointName, parentFramePos = jinfo[0], jinfo[1], jinfo[14]
            # Make map of joint name -> id
            jointId[jointName] = jointIndex
            # Get the 'd' parameter
            if jointName == b'lwing_stroke':
                urdfParams['d'] = parentFramePos[2]
                urdfParams['Roffs'] = parentFramePos[1]

        # Geometry info
        for shape in p.getVisualShapeData(self.bid):
            # All the shapes are of type p.GEOM_BOX (shape[2])
            linkId, dimensions, pos, orn = shape[1], shape[3], shape[5], shape[6]
            if linkId == jointId[b'lwing_hinge']: # link 1 and 3 are the wing membranes
                urdfParams['ycp'] = 0.5 * dimensions[0]
                urdfParams['cbar'] = dimensions[2]

        # Stiffness etc. if needed
        for j in range(-1, self.Nj):
            dinfo = p.getDynamicsInfo(self.bid, j)
            stiffness, damping = dinfo[9], dinfo[8]
            if j == jointId[b'lwing_hinge']:
                urdfParams['khinge'] = stiffness
                urdfParams['bhinge'] = damping
                # print(dinfo)
            elif j == jointId[b'lwing_stroke']:
                urdfParams['kstroke'] = stiffness

        print(urdfParams)
        return jointId, urdfParams

    def sampleStates(self):
        # get actual state
        for j in range(self.Nj):
            self.q[j], self.dq[j] = p.getJointState(self.bid, j)[0:2]
        self.q[4:7], self.q[7:11] = p.getBasePositionAndOrientation(self.bid)[0:2]
        v, omega = p.getBaseVelocity(self.bid)[0:2]
        self.dq[4:10] = self.wf.update(np.hstack((v, omega)))
        return self.simt, self.q, self.dq
        
    def resetJoints(self, bid, jarr, q, dq):
        '''Forcibly move joints to certain positions (q) and velocities (dq)'''
        for j in jarr:
            p.resetJointState(bid, j, q[j], dq[j])
            # make it so it can be torque controlled
            p.setJointMotorControl2(bid, j, controlMode=p.VELOCITY_CONTROL, targetVelocity=0, force=0)
            p.setJointMotorControl2(bid, j, controlMode=p.TORQUE_CONTROL, force=0)
    
    def visAero(self, aeroW, col, lifeTime):
        p.addUserDebugLine(aeroW[1], aeroW[1] + self.FAERO_DRAW_SCALE * np.array(aeroW[0]), lineColorRGB=col, lifeTime=lifeTime, lineWidth=2)

    def update(self, u, testF=None, forceControl=True):
        if testF is not None:
            p.setJointMotorControlArray(self.bid, [0,2], p.POSITION_CONTROL, targetPositions=[0,0], positionGains=[0.01,0.01], velocityGains=[0.1,0.1], forces=np.full(2, 1000000))
        else:
            if forceControl:
                # stroke stiffness
                for i in range(2):
                    # force in mN from voltage; add on stroke stiffness
                    u[i] = actuatorModel(u[i], self.q[2*i], self.dq[2*i]) - self.urdfParams['kstroke'] * self.q[2*i]
                p.setJointMotorControlArray(self.bid, [0,2], p.TORQUE_CONTROL, forces=u)
            else:
                p.setJointMotorControlArray(self.bid, [0,2], p.POSITION_CONTROL, targetPositions=u, positionGains=[1,1], velocityGains=[1,1], forces=np.full(2, 1000000))

        # qp, dqp = self.q[[1,3]], self.dq[[1,3]]
        # taup = -self.urdfParams['khinge'] * qp - self.urdfParams['bhinge'] * dqp
        # p.setJointMotorControlArray(self.bid, [1,3], p.TORQUE_CONTROL, forces=taup)

        # print(self.q[0:4], self.dq[:4])
        aero1B = aerodynamics(self.q[0:2], self.dq[0:2], 1, self.urdfParams)
        aero2B = aerodynamics(-self.q[2:4], -self.dq[2:4], -1, self.urdfParams)

        if testF is not None:
            aero1B = ([testF[0],0,0], aero1B[1])
            aero2B = ([testF[1],0,0], aero2B[1])
        aero1W = self.wTb(*aero1B)
        aero2W = self.wTb(*aero2B)

        # linkID = jointID
        p.applyExternalForce(self.bid, self.jointId[b'lwing_hinge'], *aero1W, p.WORLD_FRAME)
        p.applyExternalForce(self.bid, self.jointId[b'rwing_hinge'], *aero2W, p.WORLD_FRAME)

        # Bullet update
        p.stepSimulation()
        if self.simt < 1e-10 or self._camLock:
            # Reset camera to be at the correct distance (only first time)
            p.resetDebugVisualizerCamera(60, 45, -30, self.q[4:7])

        self.simt += self.TIMESTEP

        # Drawing stuff
        twall = monotonic()
        drawInterval = 1/50.
        if twall - self.tWallLastDraw > drawInterval:
            # draw debug
            self.visAero(aero1W, [1,1,0], drawInterval)
            self.visAero(aero2W, [1,0,1], drawInterval)
            self.tWallLastDraw = monotonic()
        
        if self.simt - self.tLastPrint > 0.01:
            # draw trail
            p.addUserDebugLine(self.pcomLastDraw, self.q[4:7], lineColorRGB=[0,0,1], lifeTime=0)
            self.pcomLastDraw = self.q[4:7].copy()
            self.tLastPrint = self.simt
        
        # Keyboard control options
        keys = p.getKeyboardEvents()
        # if ord('z') in keys and keys[ord('z')] & p.KEY_WAS_TRIGGERED:
        #     if self._slowDown == self.SLOWDOWN_FAST:
        #         self._slowDown = self.SLOWDOWN_SLOW
        #     else:
        #         self._slowDown = self.SLOWDOWN_FAST
        if ord('c') in keys and keys[ord('c')] & p.KEY_WAS_TRIGGERED:
            self._camLock = not self._camLock

    def openLoopLeft(self, Vamp, uoffs, f, tendMS=100, h2=0, h3=0):
        """A function for testing open-loop control input results"""
        print('Now testing:', np.array([Vamp, uoffs, f]))
        qw = []
        omega = 2 * np.pi * f
        self.reset(f)

        while self.simt < tendMS:
            ss = self.sampleStates()
            ph = omega * self.simt
            V = Vamp * ((1+h3)*np.sin(ph) + h3*np.sin(3*ph) + h2*np.sin(2*ph) + uoffs)
            tau = [V,V] # same to both wings
            # call aerodynamics to calculate avg wrench
            waero = wrench(*aerodynamics(self.q[0:2], self.dq[0:2], 1, self.urdfParams))
            # data contains wing kinematics and wrench
            qw.append(np.copy(np.hstack((self.q[:2], waero))))
            self.update(tau)
        
        return np.array(qw)