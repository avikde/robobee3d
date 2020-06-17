
import numpy as np
import pybullet as p
import pybullet_data
from scipy.spatial.transform import Rotation # TODO: eliminate

# Good parameters
CD0 = 0.4
CDmax = 3.4
CLmax = 1.8
AERO_REGULARIZE_EPS = 1e-10 # stops undefined AoA when no wind
# True in Chen (2017) science robotics, but differently calculated in Osborne (1951)
RHO = 1.225e-3 # density of air kg/m^3

rcopnondim = 0.5

def aerodynamics(theta, dtheta, lrSign, params):
    """Return aerodynamic force and instantaneous CoP in the body frame (both in R^3). If flip=False it will work for the left wing (along +y axis), and if flip=True it will """
    # Helper from DIckinson model; in order lift, drag
    Caero = lambda a : np.array([CLmax * np.sin(2*a), (CDmax + CD0) / 2.0 - (CDmax - CD0) / 2.0 * np.cos(2*a)])

    # These are all in the body frame
    Rspar = Rotation.from_euler('z', lrSign * theta[0])
    Rhinge = Rotation.from_euler('y', theta[1])
    # vector along wing
    sparVecB = Rspar.apply(np.array([0, 1, 0]))
    # wing chord unit vector
    chordB = Rspar.apply(Rhinge.apply(np.array([0,0,-1])))
    # TODO: add external wind vector
    # NOTE: assuming spar velocity is dominant
    wB = -np.cross(lrSign * dtheta[0] * np.array([0,0,1]), params['ycp'] * sparVecB)
    # print(dtheta[0] * np.array([0,0,1]), self.ycp * sparVecB)

    # COP: half of cbar down
    pcopB = np.array([0,params['Roffs'],params['d']]) + params['ycp'] * sparVecB + rcopnondim * params['cbar'] * chordB

    # Various directions in the notation of Osborne (1951)
    # l = vector along wing
    # w = relative wind vector
    # c = chord
    lwB = np.cross(sparVecB, wB)
    lwpB = wB - wB.dot(sparVecB) * sparVecB
    wnorm = np.sqrt(wB.dot(wB)) + AERO_REGULARIZE_EPS
    lwnorm = np.sqrt(lwB.dot(lwB)) + AERO_REGULARIZE_EPS
    lwpnorm = np.sqrt(lwpB.dot(lwpB)) + AERO_REGULARIZE_EPS

    # Lift/drag directions
    eD = lwpB / lwpnorm
    eL = np.array([0,0,1])

    # Calculate aero force
    aoa = np.arccos(chordB.dot(wB) / wnorm)
    Cf = Caero(aoa)
    # Cf *= 0.5 * rho * beta
    FaeroB = 0.5 * RHO * params['cbar'] * params['ycp'] * (Cf[0] * eL + Cf[1] * eD) * lwnorm**2
    if lrSign < 0:
        yflip = np.diag([1,-1,1])
        FaeroB = yflip @ FaeroB
        pcopB = yflip @ pcopB
    return FaeroB, pcopB

def actuatorModel(V, qact, dqact):
    """Return actuator force for applied voltage input and current actuator state"""
    return V # TODO:

class RobobeeSim():
    """Robobee simulator using pybullet"""

    # Parameters
    FAERO_DRAW_SCALE = 10.0
    simt = 0
    tLastDraw = 0
    tLastDraw2 = 0
    pcomLastDraw = np.zeros(3)
    tLastPrint = 0
    # vectors for storing states
    q = np.zeros(11)
    dq = np.zeros(10)

    def __init__(self, camLock=True, slowDown=True, timestep=1.0):
        self._camLock = camLock
        self._slowDown = slowDown
        self.TIMESTEP = timestep
        # Init sim
        physicsClient = p.connect(p.GUI)#or p.DIRECT for non-graphical version
        # Set up the visualizer
        p.configureDebugVisualizer(p.COV_ENABLE_WIREFRAME, 0)
        p.configureDebugVisualizer(p.COV_ENABLE_SHADOWS, 0)
        # p.configureDebugVisualizer(p.COV_ENABLE_GUI, 0)
        # p.configureDebugVisualizer(p.COV_ENABLE_MOUSE_PICKING, 0)
        p.setRealTimeSimulation(0)
        p.setTimeStep(self.TIMESTEP)
        p.setGravity(0,0,0)#-9.81e-3)
        
        # load background
        p.setAdditionalSearchPath(pybullet_data.getDataPath()) #optionally
        self.planeId = p.loadURDF("plane.urdf", globalScaling=100.0)

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
        
        # Passive hinge dynamics implemented as position control rather than joint dynamics FIXME:
        # p.setJointMotorControlArray(self.bid, [1,3], p.POSITION_CONTROL, targetPositions=[0,0], positionGains=[0.02,0.02], velocityGains=[0.02,0.02])
        # p.setJointMotorControlArray(self.bid, [1,3], p.PD_CONTROL, targetPositions=[0,0], positionGains=self.urdfParams['khinge']*np.ones(2), velocityGains=self.urdfParams['bhinge']*np.ones(2))

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
                print(dinfo)
            elif j == jointId[b'lwing_stroke']:
                urdfParams['kstroke'] = stiffness

        print(urdfParams)
        return jointId, urdfParams

    def sampleStates(self):
        # get actual state
        for j in range(self.Nj):
            self.q[j], self.dq[j] = p.getJointState(self.bid, j)[0:2]
        self.q[4:7], self.q[7:11] = p.getBasePositionAndOrientation(self.bid)[0:2]
        self.dq[4:7], self.dq[7:10] = p.getBaseVelocity(self.bid)[0:2]
        
    def resetJoints(self, bid, jarr, q, dq):
        '''Forcibly move joints to certain positions (q) and velocities (dq)'''
        for j in jarr:
            p.resetJointState(bid, j, q[j], dq[j])
            # make it so it can be torque controlled
            p.setJointMotorControl2(bid, j, controlMode=p.VELOCITY_CONTROL, targetVelocity=0, force=0)
            p.setJointMotorControl2(bid, j, controlMode=p.TORQUE_CONTROL, force=0)
    
    def visAero(self, aeroB, col):
        aeroW = self.wTb(*aeroB)
        p.addUserDebugLine(aeroW[1], aeroW[1] + self.FAERO_DRAW_SCALE * np.array(aeroW[0]), lineColorRGB=col, lifeTime=8 * self._slowDown * self.TIMESTEP)

    def update(self, u, testF=0, forceControl=False):
        if forceControl:
            # stroke stiffness
            u[0] -= self.urdfParams['kstroke'] * self.q[0]
            u[1] -= self.urdfParams['kstroke'] * self.q[2]
            # tau = [0,0]
            p.setJointMotorControlArray(self.bid, [0,2], p.TORQUE_CONTROL, forces=u)
        else:
            p.setJointMotorControlArray(self.bid, [0,2], p.POSITION_CONTROL, targetPositions=u, positionGains=[1,1], velocityGains=[1,1], forces=np.full(2, 1000000))
            # p.setJointMotorControlArray(self.bid, [0,2], p.POSITION_CONTROL, targetPositions=u, positionGains=[0.01,0.01], velocityGains=[0.1,0.1], forces=np.full(2, 1000000))
        # qp, dqp = self.q[[1,3]], self.dq[[1,3]]
        # taup = -self.urdfParams['khinge'] * qp - self.urdfParams['bhinge'] * dqp
        # p.setJointMotorControlArray(self.bid, [1,3], p.TORQUE_CONTROL, forces=taup)

        # print(self.q[0:4], self.dq[:4])
        aero1B = aerodynamics(self.q[0:2], self.dq[0:2], 1, self.urdfParams)
        aero2B = aerodynamics(-self.q[2:4], -self.dq[2:4], -1, self.urdfParams)

        if True:
            # linkID = jointID
            # FIXME: transformations https://github.com/avikde/robobee3d/pull/158
            p.applyExternalForce(self.bid, self.jointId[b'lwing_hinge'], np.diag([-1,-1,1])@aero1B[0], np.diag([-1,-1,1])@aero1B[1], p.LINK_FRAME)
            p.applyExternalForce(self.bid, self.jointId[b'rwing_hinge'], np.diag([-1,-1,1])@aero2B[0], np.diag([-1,-1,1])@aero2B[1], p.LINK_FRAME)
        else:
            # Need to convert to body frame (link -1)
            # p.applyExternalForce(bid, -1, Faeros[i], [0,0,0], p.LINK_FRAME)
            raise 'Not implemented'

        # Bullet update
        p.stepSimulation()
        if self.simt < 1e-10 or self._camLock:
            # Reset camera to be at the correct distance (only first time)
            p.resetDebugVisualizerCamera(60, 45, -30, self.q[4:7])

        self.simt += self.TIMESTEP

        # Drawing stuff
        if self.simt - self.tLastDraw > 2 * self.TIMESTEP:
            # draw debug
            self.visAero(aero1B, [1,1,0])
            self.visAero(aero2B, [1,0,1])
            self.tLastDraw = self.simt
        
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
