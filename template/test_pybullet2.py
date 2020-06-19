import pybullet as p
import pybullet_data
import time
import numpy as np

SIM_TIMESTEP = 1./240.

# Init sim
physicsClient = p.connect(p.GUI)#or p.DIRECT for non-graphical version
# Set up the visualizer
p.configureDebugVisualizer(p.COV_ENABLE_WIREFRAME, 0)
p.configureDebugVisualizer(p.COV_ENABLE_SHADOWS, 0)
p.configureDebugVisualizer(p.COV_ENABLE_GUI, 1)
p.setRealTimeSimulation(0)
p.setTimeStep(SIM_TIMESTEP)
#p.setGravity(0,0,-10)
# load robot
startPos = [0,0,0]
startOrientation = p.getQuaternionFromEuler([0,0,0])
bid = p.loadURDF("test2link.urdf", startPos, startOrientation, useFixedBase=True, flags=p.URDF_USE_INERTIA_FROM_FILE)
#cid = p.createConstraint(bid, -1, -1, -1, p.JOINT_FIXED, [0, 0, 0], [0, 0, 0], [0, 0, 1])

# springy control for both joints
# p.resetJointState(bid, 0, 0.5, 0)
# p.resetJointState(bid, 1, -0.5, 0)
p.setJointMotorControlArray(bid, [0,1], p.POSITION_CONTROL, targetPositions=[0,0],forces=[0,0], positionGains=[0.5,0.5], velocityGains=[0.1,0.1])

idff = p.addUserDebugParameter("Test force", -2, 2, 0)

while True:
    try:
        # Apply force to last link
        pend = [0,1.1,0.1]
        # world_com_pos = p.getLinkState(bid,1, computeForwardKinematics=True)[0]
        #print("world_com_pos=",world_com_pos)
        # pend = world_com_pos
        Fend = [p.readUserDebugParameter(idff),0,0]
        p.applyExternalForce(bid, 1, Fend, pend, p.WORLD_FRAME)
        p.addUserDebugLine(pend, np.array(pend) + 1*np.array(Fend), lineColorRGB=[1,0,0],lifeTime=0.1,lineWidth=2)
        p.stepSimulation()
        # time.sleep(SIM_TIMESTEP)
    except:
        raise
