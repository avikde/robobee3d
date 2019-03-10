import pybullet as p
import time
import pybullet_data
import numpy as np
import FlappingModels3D
import sys

# sim parameters
FAERO_DRAW_SCALE = 20
SIM_SLOWDOWN_SLOW = 200
SIM_SLOWDOWN_FAST = 5
SIM_SLOWDOWN = SIM_SLOWDOWN_FAST
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
p.setGravity(0,0,-10)

# load background
p.setAdditionalSearchPath(pybullet_data.getDataPath()) #optionally
planeId = p.loadURDF("plane.urdf")
# load robot
startPos = [0,0,1]
startOrientation = p.getQuaternionFromEuler([0,0,0])
bid = p.loadURDF("urdf/sdab.xacro.urdf", startPos, startOrientation, useFixedBase=False, flags=p.URDF_USE_INERTIA_FROM_FILE)
# See https://github.com/bulletphysics/bullet3/issues/2152
p.changeDynamics(bid,	-1,	maxJointVelocity=10000)
p.changeDynamics(bid, -1, linearDamping=0,	angularDamping=0)

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
for j in range(-1, Nj):
	dinfo = p.getDynamicsInfo(bid, j)
	stiffness, damping = dinfo[9], dinfo[8]
	if j == jointId[b'lwing_hinge']:
		urdfParams['stiffnessHinge'] = stiffness
		urdfParams['dampingHinge'] = damping

print(jointId)
print(urdfParams)

bee = FlappingModels3D.QuasiSteadySDAB(urdfParams)


# Simulation
simt = 0.0
tLastDraw = 0
tLastPrint = 0
# vectors for storing states
q = np.zeros(11)
dq = np.zeros(10)
bCamLock = True

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
	if simt < 1e-10 or bCamLock:
		# Reset camera to be at the correct distance (only first time)
		p.resetDebugVisualizerCamera(0.12, 45, -30, q[4:7])

def applyAero(t, q, dq, lrSign):
	pcopW, FaeroW = bee.aerodynamics(q, dq, lrSign)
	jid = jointId[b'lwing_hinge']
	if lrSign > 0:
		jid = jointId[b'rwing_hinge']
		
	p.applyExternalForce(bid, jid, FaeroW, [0, 0, 0], p.WORLD_FRAME)
	return pcopW, FaeroW

def resetAllJoints(q, dq):
	for j in range(4):
		p.resetJointState(bid, j, q[j], dq[j])
		# make it so it can be torque controlled
		# NOTE: need to comment this out if trying to reset states in the loop
		p.setJointMotorControl2(bid, j, controlMode=p.VELOCITY_CONTROL, targetVelocity=0, force=0)
		p.setJointMotorControl2(bid, j, controlMode=p.TORQUE_CONTROL, force=0)

# ---

resetAllJoints(np.zeros(4), np.zeros(4))
# Passive hinge dynamics implemented as position control rather than joint dynamics
p.setJointMotorControlArray(bid, [1,3], p.POSITION_CONTROL, targetPositions=[0,0], positionGains=urdfParams['stiffnessHinge']*np.ones(2), velocityGains=urdfParams['dampingHinge']*np.ones(2))

for i in range(10000):
	# No dynamics: reset positions
	omega = 2 * np.pi * 170.0
	ph = omega * simt
	th0 = 1.0 * np.sin(ph)
	dth0 = omega * np.cos(ph)

	# POSITION_CONTROL uses Kp, Kd in [0,1]
	p.setJointMotorControlArray(bid, [0,2], p.POSITION_CONTROL, targetPositions=[th0,th0], positionGains=[1,1], velocityGains=[0.1,0.1], forces=[1000000,1000000])

	# actual sim
	sampleStates()

	pcop1, Faero1 = applyAero(simt, q, dq, -1)
	pcop2, Faero2 = applyAero(simt, q, dq, 1)
	# applyAero(simt, 1)

	simulatorUpdate()
	time.sleep(SIM_SLOWDOWN * SIM_TIMESTEP)
	simt += SIM_TIMESTEP
	
	if simt - tLastDraw > 2 * SIM_TIMESTEP:
		# draw debug
		p.addUserDebugLine(pcop1, pcop1 + FAERO_DRAW_SCALE * Faero1, lineColorRGB=[1,1,0], lifeTime=3 * SIM_SLOWDOWN * SIM_TIMESTEP)
		p.addUserDebugLine(pcop2, pcop2 + FAERO_DRAW_SCALE * Faero2, lineColorRGB=[1,0,1], lifeTime=3 * SIM_SLOWDOWN * SIM_TIMESTEP)
		tLastDraw = simt
	
	if simt - tLastPrint > 0.01:
		tLastPrint = simt
		# print(simt, (Faero1[2] + Faero2[2]) * 1e6, q[0:2], dth0)
		print(simt, (Faero1[2] + Faero2[2]) * 1e6, dq[6])
	
	# Keyboard control options
	keys = p.getKeyboardEvents()
	if ord('z') in keys and keys[ord('z')] & p.KEY_WAS_TRIGGERED:
		if SIM_SLOWDOWN == SIM_SLOWDOWN_FAST:
			SIM_SLOWDOWN = SIM_SLOWDOWN_SLOW
		else:
			SIM_SLOWDOWN = SIM_SLOWDOWN_FAST
	if ord('c') in keys and keys[ord('c')] & p.KEY_WAS_TRIGGERED:
		bCamLock = not bCamLock

p.disconnect()
