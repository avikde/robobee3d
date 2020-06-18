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
bid = p.loadURDF("test2link.urdf", startPos, startOrientation, useFixedBase=True)

# springy control for both joints
# p.resetJointState(bid, 0, 0.5, 0)
p.resetJointState(bid, 1, -0.5, 0)
p.setJointMotorControlArray(bid, [0,1], p.POSITION_CONTROL, targetPositions=[0,0], positionGains=[0.01,1], velocityGains=[0.1,0.99])

idff = p.addUserDebugParameter("Test force", -0.01, 0.01, 0)

while True:
    try:
        # Apply force to last link
        pend = [0,1.1,0.1]
        Fend = [p.readUserDebugParameter(idff),0,0]
        p.applyExternalForce(bid, 1, Fend, pend, p.WORLD_FRAME)
        p.addUserDebugLine(pend, np.array(pend) + 100*np.array(Fend), lineColorRGB=[1,0,0],lifeTime=0.1,lineWidth=2)
        p.stepSimulation()
        time.sleep(SIM_TIMESTEP * 1e-3)
    except:
        raise
