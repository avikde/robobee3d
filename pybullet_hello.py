import pybullet as p
import time
import pybullet_data
import numpy as np

physicsClient = p.connect(p.GUI)#or p.DIRECT for non-graphical version
p.setAdditionalSearchPath(pybullet_data.getDataPath()) #optionally
# p.setGravity(0,0,-10)
planeId = p.loadURDF("plane.urdf")
startPos = [0,0,1]
startOrientation = p.getQuaternionFromEuler([0,0,0])

bid = p.loadURDF("urdf/robobee.urdf", startPos, startOrientation)

# Get info about urdf
Nj = p.getNumJoints(bid)
jointId = {}
for j in range(Nj):
	jinfo = p.getJointInfo(bid, j)
	jointId[jinfo[1]] = jinfo[0]
print(jointId)

#  Since each link is connected to a parent with a single joint,
# the number of joints is equal to the number of links. Regular links have link indices in the range
# [0..getNumJoints()] Since the base is not a regular 'link', we use the convention of -1 as its link
# index

# Simulation
t = 0.0
dt = 0.001

def aerodynamics(bRight):
	pwingB = 0
	# pwing = pcom + 

	pcopW = np.array([0, 0, 0])
	FaeroW = np.array([0, 0, 0])

	return pcopW, FaeroW

def applyAero(bRight):
	pcopW, Faero = aerodynamics(bRight)
	# p.appyExternalForce(bid, jointId[b'lwing_hinge'], [0, 0, 0], [0, 0, 0], p.WORLD_FRAME)

	# draw debug
	red = [1, 0, 0]
	p.addUserDebugLine(pcopW, pcopW + 0.1 * FaeroW, lineColorRGB=red, lifeTime=0.3)



for i in range(10000):
	# No dynamics: reset positions
	ph = 2 * np.pi * t
	p.resetJointState(bid, jointId[b'lwing_stroke'], 0.5 * np.sin(ph))
	p.resetJointState(bid, jointId[b'rwing_stroke'], 0.5 * np.sin(ph))
	p.resetJointState(bid, jointId[b'lwing_hinge'], 1 * np.cos(ph))
	p.resetJointState(bid, jointId[b'rwing_hinge'], -1 * np.cos(ph))

	# actual sim
	# applyAero(0)
	# applyAero(1)

	p.stepSimulation()
	time.sleep(dt)
	t += dt

	# cubePos, cubeOrn = p.getBasePositionAndOrientation(bid)
	# print(cubePos,cubeOrn)
p.disconnect()
