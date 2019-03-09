import pybullet as p
import pybullet_data
import time
import numpy as np

SIM_TIMESTEP = 0.0001

# Init sim
physicsClient = p.connect(p.GUI)#or p.DIRECT for non-graphical version
# Set up the visualizer
p.configureDebugVisualizer(p.COV_ENABLE_WIREFRAME, 0)
p.configureDebugVisualizer(p.COV_ENABLE_SHADOWS, 0)
p.configureDebugVisualizer(p.COV_ENABLE_GUI, 0)
# p.configureDebugVisualizer(p.COV_ENABLE_MOUSE_PICKING, 0)
p.setRealTimeSimulation(0)
p.setTimeStep(SIM_TIMESTEP)


p.setAdditionalSearchPath(pybullet_data.getDataPath()) #optionally
planeId = p.loadURDF("plane.urdf")
# load robot
startPos = [0,0,1]
startOrientation = p.getQuaternionFromEuler([0,0,0])
bid = p.loadURDF("urdf/TwoJointRobot.urdf", startPos, startOrientation, useFixedBase=True)

simt = 0
ph = 0
tLastPrint = 0

while simt < 0.5:
	if simt < 0.2:
		freq = 50 # Hz
	elif simt < 0.4:
		freq = 100
	elif simt < 0.6:
		freq = 200
	elif simt < 0.8:
		freq = 400
	else:
		freq = 1000
	# No dynamics: reset positions
	ph = ph + 2 * np.pi * freq * SIM_TIMESTEP
	qdes = 0.5 * np.sin(ph)
	dqdes = 0.5 * 2 * np.pi * freq * np.cos(ph)
	# joint 1 - second link, joint 0 is fixed
	p.setJointMotorControlArray(bid, [1], p.POSITION_CONTROL, targetPositions=[qdes], targetVelocities=[dqdes], positionGains=[10000], velocityGains=[100], forces=[1000000000])

	p.stepSimulation()
	time.sleep(SIM_TIMESTEP)
	simt += SIM_TIMESTEP
	
	if simt - tLastPrint > 0.0005:
		qact = p.getJointState(bid, 1)[0]
		tLastPrint = simt
		print(simt, qdes, qact)

p.disconnect()

