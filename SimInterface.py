'''Interface to sim engines such as pybullet, pydrake

TODO:
- check if import works
- pydrake interface

This is the only file that will have sim engine-specific imports
'''
import numpy as np
import pybullet as p
import pybullet_data

class PyBullet():
	# Parameters
	bCamLock = True
	SLOWDOWN_SLOW = 200
	SLOWDOWN_FAST = 5
	SLOWDOWN = SLOWDOWN_FAST
	TIMESTEP = 0.0001
	FAERO_DRAW_SCALE = 20
	simt = 0
	tLastDraw = 0
	tLastDraw2 = 0
	pcomLastDraw = np.zeros(3)
	tLastPrint = 0
	# vectors for storing states
	q = np.zeros(11)
	dq = np.zeros(10)
	# 
	POSITION_CONTROL = p.POSITION_CONTROL
	TORQUE_CONTROL = p.TORQUE_CONTROL

	def __init__(self):
		# Init sim
		physicsClient = p.connect(p.GUI)#or p.DIRECT for non-graphical version
		# Set up the visualizer
		p.configureDebugVisualizer(p.COV_ENABLE_WIREFRAME, 0)
		p.configureDebugVisualizer(p.COV_ENABLE_SHADOWS, 0)
		p.configureDebugVisualizer(p.COV_ENABLE_GUI, 0)
		# p.configureDebugVisualizer(p.COV_ENABLE_MOUSE_PICKING, 0)
		p.setRealTimeSimulation(0)
		p.setTimeStep(self.TIMESTEP)
		p.setGravity(0,0,-10)
		
		# load background
		p.setAdditionalSearchPath(pybullet_data.getDataPath()) #optionally
		self.planeId = p.loadURDF("plane.urdf")

	def loadURDF(self, filename, startPos, startOrientation, useFixedBase):
		'''startPos = xyz
		startOrientation = quaternion xyzw
		'''
		bid = p.loadURDF(filename, startPos, startOrientation, useFixedBase=useFixedBase, flags=p.URDF_USE_INERTIA_FROM_FILE)
		# See https://github.com/bulletphysics/bullet3/issues/2152
		p.changeDynamics(bid,	-1,	maxJointVelocity=10000)
		p.changeDynamics(bid, -1, linearDamping=0,	angularDamping=0)
		
		self.Nj = p.getNumJoints(bid)
		self.pcomLastDraw = startPos

		return bid

	def getInfoFromURDF(self, bid):
		'''Get info from the URDF
		'''
		urdfParams = {}

		#  Since each link is connected to a parent with a single joint,
		# the number of joints is equal to the number of links. Regular links have link indices in the range
		# [0..getNumJoints()] Since the base is not a regular 'link', we use the convention of -1 as its link
		# index
		jointId = {}
		for j in range(self.Nj):
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
		for j in range(-1, self.Nj):
			dinfo = p.getDynamicsInfo(bid, j)
			stiffness, damping = dinfo[9], dinfo[8]
			if j == jointId[b'lwing_hinge']:
				urdfParams['stiffnessHinge'] = stiffness
				urdfParams['dampingHinge'] = damping
			elif j == jointId[b'lwing_stroke']:
				urdfParams['stiffnessStroke'] = stiffness
				urdfParams['dampingStroke'] = damping

		return jointId, urdfParams

	def sampleStates(self, bid):
		# get actual state
		for j in range(self.Nj):
			self.q[j], self.dq[j] = p.getJointState(bid, j)[0:2]
		self.q[4:7], self.q[7:11] = p.getBasePositionAndOrientation(bid)[0:2]
		self.dq[4:7], self.dq[7:10] = p.getBaseVelocity(bid)[0:2]
		
	def resetJoints(self, bid, jarr, q, dq):
		'''Forcibly move joints to certain positions (q) and velocities (dq)'''
		for j in jarr:
			p.resetJointState(bid, j, q[j], dq[j])
			# make it so it can be torque controlled
			# NOTE: need to comment this out if trying to reset states in the loop
			p.setJointMotorControl2(bid, j, controlMode=p.VELOCITY_CONTROL, targetVelocity=0, force=0)
			p.setJointMotorControl2(bid, j, controlMode=p.TORQUE_CONTROL, force=0)

	def applyExternalForce(self, bid, jid, force, pos):
		p.applyExternalForce(bid, jid, force, pos, p.WORLD_FRAME)

	def disconnect(self):
		p.disconnect()

	def setJointArray(self, *args, **kwargs):
		p.setJointMotorControlArray(*args, **kwargs)

	def addUserDebugLine(self, *args, **kwargs):
		p.addUserDebugLine(*args, **kwargs)

	def update(self, bid, jointIndices, pcops, Faeros):
		for i in range(2):
		# FIXME: 0 and not pcop?
			self.applyExternalForce(bid, jointIndices[i], Faeros[i], [0,0,0])

		# Bullet update
		p.stepSimulation()
		self.simt += self.TIMESTEP
		if self.simt < 1e-10 or self.bCamLock:
			# Reset camera to be at the correct distance (only first time)
			p.resetDebugVisualizerCamera(0.12, 45, -30, self.q[4:7])

		# Drawing stuff
		if self.simt - self.tLastDraw > 2 * self.TIMESTEP:
			# draw debug
			cols = [[1,1,0], [1,0,1]]
			for i in range(2):
				p.addUserDebugLine(pcops[i], pcops[i] + self.FAERO_DRAW_SCALE * Faeros[i], lineColorRGB=cols[i], lifeTime=3 * self.SLOWDOWN * self.TIMESTEP)
			self.tLastDraw = self.simt
		
		if self.simt - self.tLastPrint > 0.01:
			# draw trail
			p.addUserDebugLine(self.pcomLastDraw, self.q[4:7], lineColorRGB=[0,0,1], lifeTime=0)
			self.pcomLastDraw = self.q[4:7].copy()
			print(self.simt, (Faeros[0][2] + Faeros[1][2]) * 1e6)
			self.tLastPrint = self.simt
		
		# Keyboard control options
		keys = p.getKeyboardEvents()
		if ord('z') in keys and keys[ord('z')] & p.KEY_WAS_TRIGGERED:
			if self.SLOWDOWN == self.SLOWDOWN_FAST:
				self.SLOWDOWN = self.SLOWDOWN_SLOW
			else:
				self.SLOWDOWN = self.SLOWDOWN_FAST
		if ord('c') in keys and keys[ord('c')] & p.KEY_WAS_TRIGGERED:
			self.bCamLock = not self.bCamLock
