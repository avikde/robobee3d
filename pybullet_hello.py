import pybullet as p
import time
import pybullet_data
import numpy as np
import FlappingModels3D

physicsClient = p.connect(p.GUI)#or p.DIRECT for non-graphical version
p.setAdditionalSearchPath(pybullet_data.getDataPath()) #optionally
# p.setGravity(0,0,-10)
planeId = p.loadURDF("plane.urdf")
startPos = [0,0,1]
startOrientation = p.getQuaternionFromEuler([0,0,0])

bid = p.loadURDF("urdf/sdab.xacro.urdf", startPos, startOrientation, useFixedBase=True)

# Get info about urdf
Nj = p.getNumJoints(bid)
jointId = {}
for j in range(Nj):
	jinfo = p.getJointInfo(bid, j)
	jointId[jinfo[1]] = jinfo[0]
print(jointId)

# vectors for storing states
q = np.zeros(11)
dq = np.zeros(10)

bee = FlappingModels3D.QuasiSteadySDAB()

#  Since each link is connected to a parent with a single joint,
# the number of joints is equal to the number of links. Regular links have link indices in the range
# [0..getNumJoints()] Since the base is not a regular 'link', we use the convention of -1 as its link
# index

# Simulation
simt = 0.0
dt = 0.001

tLastDraw = 0

# Helpers ---

def sampleStates():
	# get actual state
	for j in range(Nj):
		q[j], dq[j] = p.getJointState(bid, j)[0:2]
	q[4:7], q[7:11] = p.getBasePositionAndOrientation(bid)[0:2]
	dq[4:7], dq[7:10] = p.getBaseVelocity(bid)[0:2]

def applyAero(t, q, dq, bRight):
	pcopW, FaeroW = bee.aerodynamics(q, dq, bRight)
	# p.appyExternalForce(bid, jointId[b'lwing_hinge'], [0, 0, 0], [0, 0, 0], p.WORLD_FRAME)
	return pcopW, FaeroW

# ---

for i in range(10000):
	# No dynamics: reset positions
	ph = 2 * np.pi * simt
	th0 = 0.5 * np.sin(ph)
	dth0 = np.pi * np.cos(ph)
	th1 = np.cos(ph)
	dth1 = -2 * np.pi * np.sin(ph)
	p.resetJointState(bid, jointId[b'lwing_stroke'], th0, dth0)
	p.resetJointState(bid, jointId[b'rwing_stroke'], th0, dth0)
	p.resetJointState(bid, jointId[b'lwing_hinge'], th1, dth1)
	p.resetJointState(bid, jointId[b'rwing_hinge'], -th1, -dth1)

	# actual sim
	sampleStates()

	pcop1, Faero1 = applyAero(simt, q, dq, -1)
	pcop2, Faero2 = applyAero(simt, q, dq, 1)
	# applyAero(simt, 1)

	p.stepSimulation()
	time.sleep(dt)
	simt += dt
	
	if simt - tLastDraw > 0.05:
		# draw debug
		red = [1, 1, 0]
		p.addUserDebugLine(pcop1, pcop1 + 0.5 * Faero1, lineColorRGB=red, lifeTime=0.2)
		p.addUserDebugLine(pcop2, pcop2 + 0.5 * Faero2, lineColorRGB=[1,0,1], lifeTime=0.2)
		tLastDraw = simt
		# print(simt, dth0)

p.disconnect()
