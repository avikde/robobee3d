import pybullet as p
import time
import pybullet_data
import numpy as np
import FlappingModels3D

# sim parameters
FAERO_DRAW_SCALE = 20
SIM_SLOWDOWN = 500

# Init sim
physicsClient = p.connect(p.GUI)#or p.DIRECT for non-graphical version
# Set up the visualizer
p.configureDebugVisualizer(p.COV_ENABLE_WIREFRAME, 0)
p.configureDebugVisualizer(p.COV_ENABLE_SHADOWS, 0)
p.configureDebugVisualizer(p.COV_ENABLE_GUI, 0)
p.configureDebugVisualizer(p.COV_ENABLE_MOUSE_PICKING, 0)

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


# TODO: get the params from the URDF
bee = FlappingModels3D.QuasiSteadySDAB(0.006, 0.006, 0.005)

#  Since each link is connected to a parent with a single joint,
# the number of joints is equal to the number of links. Regular links have link indices in the range
# [0..getNumJoints()] Since the base is not a regular 'link', we use the convention of -1 as its link
# index

# Simulation
simt = 0.0
dt = 0.0001
tLastDraw = 0
# vectors for storing states
q = np.zeros(11)
dq = np.zeros(10)

# Helpers ---

def sampleStates():
	# get actual state
	for j in range(Nj):
		q[j], dq[j] = p.getJointState(bid, j)[0:2]
	q[4:7], q[7:11] = p.getBasePositionAndOrientation(bid)[0:2]
	dq[4:7], dq[7:10] = p.getBaseVelocity(bid)[0:2]
	
def simulatorUpdate():
	# Bullet update
	p.stepSimulation()
	# Reset camera to track
	p.resetDebugVisualizerCamera(0.15, 45, -30, q[4:7])

def applyAero(t, q, dq, bRight):
	pcopW, FaeroW = bee.aerodynamics(q, dq, bRight)
	# p.appyExternalForce(bid, jointId[b'lwing_hinge'], [0, 0, 0], [0, 0, 0], p.WORLD_FRAME)
	return pcopW, FaeroW

# ---

for i in range(10000):
	# No dynamics: reset positions
	freq = 170.0
	ph = 2 * np.pi * freq * simt
	th0 = 0.5 * np.sin(ph)
	dth0 = np.pi * freq * np.cos(ph)
	th1 = np.cos(ph)
	dth1 = -2 * np.pi * freq * np.sin(ph)
	p.resetJointState(bid, jointId[b'lwing_stroke'], th0, dth0)
	p.resetJointState(bid, jointId[b'rwing_stroke'], th0, dth0)
	p.resetJointState(bid, jointId[b'lwing_hinge'], th1, dth1)
	p.resetJointState(bid, jointId[b'rwing_hinge'], th1, dth1)

	# actual sim
	sampleStates()

	pcop1, Faero1 = applyAero(simt, q, dq, -1)
	pcop2, Faero2 = applyAero(simt, q, dq, 1)
	# applyAero(simt, 1)

	simulatorUpdate()
	time.sleep(SIM_SLOWDOWN * dt)
	simt += dt
	
	if simt - tLastDraw > 2 * dt:
		# draw debug
		red = [1, 1, 0]
		p.addUserDebugLine(pcop1, pcop1 + FAERO_DRAW_SCALE * Faero1, lineColorRGB=red, lifeTime=3 * SIM_SLOWDOWN * dt)
		p.addUserDebugLine(pcop2, pcop2 + FAERO_DRAW_SCALE * Faero2, lineColorRGB=[1,0,1], lifeTime=3 * SIM_SLOWDOWN * dt)
		tLastDraw = simt
		print("total lift =", (Faero1[2] + Faero2[2]) * 1e6)

p.disconnect()
