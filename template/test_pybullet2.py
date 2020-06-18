import pybullet as p
import pybullet_data
import time
import numpy as np

SIM_TIMESTEP = 1

# Init sim
physicsClient = p.connect(p.GUI)#or p.DIRECT for non-graphical version
# Set up the visualizer
p.configureDebugVisualizer(p.COV_ENABLE_WIREFRAME, 0)
p.configureDebugVisualizer(p.COV_ENABLE_SHADOWS, 0)
p.configureDebugVisualizer(p.COV_ENABLE_GUI, 1)
p.setRealTimeSimulation(0)
p.setTimeStep(SIM_TIMESTEP)

# load robot
startPos = [0,0,0]
startOrientation = p.getQuaternionFromEuler([0,0,0])
bid2 = p.loadURDF("test2link.urdf", startPos, startOrientation, useFixedBase=True)

simt = 0

idff = p.addUserDebugParameter("Test force", -1, 1, 0)

while True:
    try:
        p.stepSimulation()
        time.sleep(SIM_TIMESTEP * 1e-3)
    except:
        raise
