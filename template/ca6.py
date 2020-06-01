import time
import gzip, pickle
import numpy as np
import sys
from scipy.spatial.transform import Rotation
import pybullet as p
import pybullet_data
import viewlog
from wrenchlinQP import WrenchLinQP

# Usage params
TIMESTEP = 1
SLOWDOWN = 0.01
FAERO_DRAW_SCALE = 10.0
# model params
ycp = 10 # mm

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
startPos = [0,0,20]
startOrientation = Rotation.from_euler('xy', [0.5,0.5])
bid = p.loadURDF("../urdf/sdabNW.urdf", startPos, startOrientation.as_quat(), useFixedBase=False)

simt = 0
tLastDraw1 = 0
data = viewlog.initLog()

# controller
wlqp = WrenchLinQP(6,6)
wlqp.test()

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
    # self.dq[4:7], self.dq[7:10] = p.getBaseVelocity(bid)[0:2]
    Rb = Rotation.from_quat(qorn)
    return pcom, Rb

def wTb(pw, Rb, Fb, pb):
    "Transform wrench to world frame"
    return Rb.apply(Fb), pw + Rb.apply(pb)

def transformWrenches(pw, Rb, FL, pL, FR, pR):
    # join tuples with +
    return wTb(pw, Rb, FL, pL) + wTb(pw, Rb, FR, pR)

def testControl(pw, Rb):
    # u = [1,0,0,1,0,0]
    mm = 1.0
    ezb = Rb.apply([0,0,1])
    dd = ezb[1] * (1.0)
    pitchCtrl = ezb[0]
    u = [mm + dd, pitchCtrl,0.0,mm-dd,pitchCtrl,-0.0]
    return u

while True:
    try:
        ss = getState()
        u = testControl(*ss)
        data = viewlog.appendLog(data, simt, *ss, u)
        wrenchesB = ca6ApplyInput(u)
        # Bullet update
        p.stepSimulation()
        simt += TIMESTEP
        
        # Other updates
        wrenchesW = transformWrenches(*ss, *wrenchesB)

        if True:#simt < 1 or _camLock:
            # Reset camera to be at the correct distance (only first time)
            p.resetDebugVisualizerCamera(100, 45, -30, ss[0])

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
        break
    
viewlog.saveLog('../logs/', data)
# p.disconnect()
