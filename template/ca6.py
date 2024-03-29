import time
import gzip, pickle
import numpy as np
import sys
from scipy.spatial.transform import Rotation
import pybullet as p
import pybullet_data
import viewlog
from wrenchlinQP import WrenchLinQP
from ca6dynamics import ycp, dynamicsTerms, wrenchMap
np.set_printoptions(precision=2, suppress=True, linewidth=200)

# Usage params
TIMESTEP = 1
SLOWDOWN = 0.01
FAERO_DRAW_SCALE = 10.0

# Init sim
physicsClient = p.connect(p.GUI)#or p.DIRECT for non-graphical version
# Set up the visualizer
p.configureDebugVisualizer(p.COV_ENABLE_WIREFRAME, 0)
p.configureDebugVisualizer(p.COV_ENABLE_SHADOWS, 0)
p.configureDebugVisualizer(p.COV_ENABLE_GUI, 0)
# p.configureDebugVisualizer(p.COV_ENABLE_MOUSE_PICKING, 0)
p.setRealTimeSimulation(0)
p.setTimeStep(TIMESTEP)
p.setGravity(0,0,-9.81e-3)

# load background
p.setAdditionalSearchPath(pybullet_data.getDataPath()) #optionally
planeId = p.loadURDF("plane.urdf", globalScaling=100.0)

# load robot
startPos = [0,0,100]
startOrientation = Rotation.from_euler('xy', [0.8,0.8])
bid = p.loadURDF("../urdf/sdabNW.urdf", startPos, startOrientation.as_quat(), useFixedBase=False)

simt = 0
tLastDraw1 = 0
data = viewlog.initLog()

# controller
wlqp = WrenchLinQP(6, 6, dynamicsTerms, wrenchMap)
# wlqp = WrenchLinQP(1,1)
# wlqp.test()

def ca6ApplyInput(u):
    """Input vector https://github.com/avikde/robobee3d/pull/154"""
    u1L, u2L, u3L, u1R, u2R, u3R = u # unpack
    FL, pL = [u3L,0,u1L], [u2L,ycp,0]
    FR, pR = [u3R,0,u1R], [u2R,-ycp,0]
    p.applyExternalForce(bid, -1, FL, pL, p.LINK_FRAME)
    p.applyExternalForce(bid, -1, FR, pR, p.LINK_FRAME)
    return FL, pL, FR, pR

def ca6DrawInput(FL, pL, FR, pR):
    lt = TIMESTEP * 0.04
    lw = 2
    p.addUserDebugLine(pL, np.array(pL) + FAERO_DRAW_SCALE * np.array(FL), lineColorRGB=[1,0,1], lifeTime=lt, lineWidth=lw)
    p.addUserDebugLine(pR, np.array(pR) + FAERO_DRAW_SCALE * np.array(FR), lineColorRGB=[1,0,0], lifeTime=lt, lineWidth=lw)

def getState():
    # quat is in xyzw order (same as scipy.Rotation)
    pcom, qorn = p.getBasePositionAndOrientation(bid)[0:2]
    Rb = Rotation.from_quat(qorn)
    # Get velocity
    vW, omegaW = p.getBaseVelocity(bid)[0:2]
    # convert to body frame
    vB = Rb.inv().apply(vW)
    omegaB = Rb.inv().apply(omegaW)
    return np.hstack((pcom, qorn)), np.hstack((vB, omegaB))

def wTb(q, Fb, pb):
    "Transform wrench to world frame"
    pw, Rb = q[:3], Rotation.from_quat(q[3:7])
    return Rb.apply(Fb), pw + Rb.apply(pb)

def transformWrenches(q, FL, pL, FR, pR):
    # join tuples with +
    return wTb(q, FL, pL) + wTb(q, FR, pR)

def testControl(q, dq, pdes):
    # u = [1,0,0,1,0,0]

    # mm = 1.0
    # ezb = Rb.apply([0,0,1])
    # dd = ezb[1] * (1.0)
    # pitchCtrl = ezb[0]
    # u = [mm + dd, pitchCtrl,0.0,mm-dd,pitchCtrl,-0.0]

    u = wlqp.updateFromState(0., q, dq, pdes)

    return u

def ornVF(Rb, omega, posErr):
    """return last elements of pdes"""
    # hat = lambda M : np.array([M[2,1], M[0,2], M[1,0]])
    Rdes = np.eye(3)#Rotation.from_euler('x', 0)
    Rm = Rb.as_matrix()
    zdes = np.array([0,0,1]) # desired z vector
    # ornError = hat(Rdes.T @ Rm - Rm.T @ Rdes)
    # ornError = -np.cross([0,0,1], Rb.inv().apply(zdes))
    ornError = np.array([[Rm[0,1], Rm[1,1], Rm[2,1]], [-Rm[0,0], -Rm[1,0], -Rm[2,0]], [0,0,0]]) @ zdes # Pakpong (2013) (6)
    # ornError = ornError / np.linalg.norm(ornError)
    # print(ornError, test)
    Iomegades = -20.0*ornError - 1.0*omega
    return Iomegades

while True:
    try:
        q, dq = getState()
        posErr = np.array([-20,20,0]) - q[:3]
        pdes = np.hstack(([0,0,10], ornVF(Rotation.from_quat(q[3:7]), dq[3:], posErr)))
        u = testControl(q, dq, pdes)
        data = viewlog.appendLog(data, simt, q, dq, u, pdes)
        wrenchesB = ca6ApplyInput(u)
        # Bullet update
        p.stepSimulation()
        simt += TIMESTEP
        
        # Other updates
        wrenchesW = transformWrenches(q, *wrenchesB)

        if True:#simt < 1 or _camLock:
            # Reset camera to be at the correct distance (only first time)
            p.resetDebugVisualizerCamera(100, 45, -30, q[:3])

        # Drawing stuff
        if simt - tLastDraw1 > 2 * TIMESTEP:
            ca6DrawInput(*wrenchesW)
            tLastDraw = simt
        
        # if self.simt - self.tLastPrint > 0.01:
        #     # draw trail
        #     p.addUserDebugLine(self.pcomLastDraw, self.q[4:7], lineColorRGB=[0,0,1], lifeTime=0)
        #     self.pcomLastDraw = self.q[4:7].copy()
        #     self.tLastPrint = self.simt
            
        time.sleep(SLOWDOWN * TIMESTEP)
    
        # if simt - tLastPrint > 10:
        #     # print(sim.simt, posErr)
        #     tLastPrint = simt
    
    except:
        viewlog.saveLog('../logs/ca6', data)
        raise
    
# p.disconnect()
