import pybullet as p
import time
import pybullet_data
import numpy as np
import FlappingModels3D

# sim parameters
FAERO_DRAW_SCALE = 100000
SIM_SLOWDOWN = 500

# Init sim
physicsClient = p.connect(p.GUI)#or p.DIRECT for non-graphical version
# Set up the visualizer
p.configureDebugVisualizer(p.COV_ENABLE_WIREFRAME, 0)
p.configureDebugVisualizer(p.COV_ENABLE_SHADOWS, 0)
p.configureDebugVisualizer(p.COV_ENABLE_GUI, 0)
# p.configureDebugVisualizer(p.COV_ENABLE_MOUSE_PICKING, 0)

p.setAdditionalSearchPath(pybullet_data.getDataPath()) #optionally
# p.setGravity(0,0,-10)
planeId = p.loadURDF("plane.urdf")
startPos = [0,0,1]
startOrientation = p.getQuaternionFromEuler([0,0,0])

bid = p.loadURDF("urdf/sdab.xacro.urdf", startPos, startOrientation, useFixedBase=True, flags=p.URDF_USE_INERTIA_FROM_FILE)

# Get info from the URDF
urdfParams = {}

#  Since each link is connected to a parent with a single joint,
# the number of joints is equal to the number of links. Regular links have link indices in the range
# [0..getNumJoints()] Since the base is not a regular 'link', we use the convention of -1 as its link
# index
Nj = p.getNumJoints(bid)
jointId = {}
for j in range(Nj):
	jinfo = p.getJointInfo(bid, j)
	# unpack
	jointIndex, jointName, parentFramePos = jinfo[0], jinfo[1], jinfo[14]
	# Make map of joint name -> id
	jointId[jointName] = jointIndex
	# Get the 'd' parameter
	if jointName == b'lwing_stroke':
		urdfParams['d'] = parentFramePos[2]

# Geometry info
for shape in p.getVisualShapeData(bid):
	# All the shapes are of type p.GEOM_BOX (shape[2])
	linkId, dimensions, pos, orn = shape[1], shape[3], shape[5], shape[6]
	if linkId == jointId[b'lwing_hinge']: # link 1 and 3 are the wing membranes
		urdfParams['rcp'] = 0.5 * dimensions[1]
		urdfParams['cbar'] = dimensions[2]

# Stiffness etc. if needed
# for j in range(-1, Nj):
# 	print("Link", j, "info:", p.getDynamicsInfo(bid, j))

print(jointId)
print(urdfParams)

bee = FlappingModels3D.QuasiSteadySDAB(urdfParams)


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

def resetAllJoints(q, dq):
	for j in range(4):
		p.resetJointState(bid, j, q[j], dq[j])
		# make it so it can be torque controlled
		# NOTE: need to comment this out if trying to reset states in the loop
		p.setJointMotorControl2(bid, j, controlMode=p.VELOCITY_CONTROL, force=0)

# ---

resetAllJoints(np.zeros(4), np.zeros(4))

for i in range(10000):
	# No dynamics: reset positions
	freq = 170.0
	ph = 2 * np.pi * freq * simt
	th0 = 0.5 * np.sin(ph)
	dth0 = np.pi * freq * np.cos(ph)
	th1 = np.cos(ph)
	dth1 = -2 * np.pi * freq * np.sin(ph)
	# resetAllJoints([th0,th1,th0,th1], [dth0,dth1,dth0,dth1])

	p.setJointMotorControlArray(bid, [0,2], p.POSITION_CONTROL, targetPositions=[th0,th0], positionGains=[1,1], velocityGains=[1,1])
	# TODO: this should be passive
	p.setJointMotorControlArray(bid, [1,3], p.POSITION_CONTROL, targetPositions=[th1,th1], positionGains=[1,1], velocityGains=[1,1])

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
		print("total lift =", (Faero1[2] + Faero2[2]) * 1e6,'dth =',dq[0:2])

p.disconnect()
