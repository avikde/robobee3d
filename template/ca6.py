import time
import numpy as np
import sys
from scipy.spatial.transform import Rotation
import pybullet as p
import pybullet_data
# import SimInterface
# import FlappingModels3D

# Usage params
TIMESTEP = 1
SLOWDOWN = 0.01
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
startOrientation = Rotation.from_euler('x', 0.5)
bid = p.loadURDF("../urdf/sdabNW.urdf", startPos, startOrientation.as_quat(), useFixedBase=False)
p.resetDebugVisualizerCamera(100, 45, -30, [0,0,0])#q[4:7])

simt = 0
tLastPrint = 0

def ca6ApplyInput(u):
    """Input vector https://github.com/avikde/robobee3d/pull/154"""
    u1L, u2L, u3L, u1R, u2R, u3R = u # unpack
    p.applyExternalForce(bid, -1, [u3L,0,u1L], [u2L,ycp,0], p.LINK_FRAME)
    p.applyExternalForce(bid, -1, [u3R,0,u1R], [u2R,-ycp,0], p.LINK_FRAME)

while True:
    try:
        ca6ApplyInput([10,0,0,10,0,0])
        # Bullet update
        p.stepSimulation()
        # if simt < 1 or _camLock:
        #     # Reset camera to be at the correct distance (only first time)
        #     

        simt += TIMESTEP

        # # Drawing stuff
        # if simt - tLastDraw > 2 * self.TIMESTEP:
        #     # draw debug
        #     cols = [[1,1,0], [1,0,1]]
        #     for i in range(2):
        #         p.addUserDebugLine(pcops[i], pcops[i] + self.FAERO_DRAW_SCALE * Faeros[i], lineColorRGB=cols[i], lifeTime=3 * self._slowDown * self.TIMESTEP)
        #     self.tLastDraw = self.simt
        
        # if self.simt - self.tLastPrint > 0.01:
        #     # draw trail
        #     p.addUserDebugLine(self.pcomLastDraw, self.q[4:7], lineColorRGB=[0,0,1], lifeTime=0)
        #     self.pcomLastDraw = self.q[4:7].copy()
        #     self.tLastPrint = self.simt
            
        time.sleep(SLOWDOWN * TIMESTEP)
    
        if simt - tLastPrint > 10:
            # print(sim.simt, posErr)
            tLastPrint = simt
    
    except KeyboardInterrupt:
        break
    
p.disconnect()
