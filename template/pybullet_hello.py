import time
import numpy as np
import sys
from scipy.spatial.transform import Rotation
import SimInterface
import FlappingModels3D

# Usage params
STROKE_FORCE_CONTROL = False # if false, use position control on the stroke
T_END = 1
AERO_WORLD_FRAME = True

sim = SimInterface.PyBullet(slowDown=False, camLock=True)
# load robot
startPos = [0,0,1]
startOrientation = Rotation.from_euler('x', 0)
bid = sim.loadURDF("../urdf/sdab.xacro.urdf", startPos, startOrientation.as_quat(), useFixedBase=False)
jointId, urdfParams = sim.getInfoFromURDF(bid)
print(jointId, urdfParams)
# Our model for control
bee = FlappingModels3D.QuasiSteadySDAB(urdfParams)

# Helper function: traj to track
def traj(t):
	return startPos + np.array([0.03 * np.sin(2*np.pi*t), 0, 0.5 * t])

# draw traj
tdraw = np.linspace(0, T_END, 20)
for ti in range(1, len(tdraw)):
	sim.addUserDebugLine(traj(tdraw[ti-1]), traj(tdraw[ti]), lineColorRGB=[0,0,0], lifeTime=0)
	
# ---

# Get the joints to reasonable positions
sim.resetJoints(bid, range(4), np.zeros(4), np.zeros(4))
# Passive hinge dynamics implemented as position control rather than joint dynamics
sim.setJointArray(bid, [1,3], sim.POSITION_CONTROL, targetPositions=[0,0], positionGains=urdfParams['stiffnessHinge']*np.ones(2), velocityGains=urdfParams['dampingHinge']*np.ones(2))

# conventional controller params
ctrl = {'thrust': 4e-5, 'strokedev': 0, 'ampl': 1.0, 'freq': 170}
tLastPrint = 0

while sim.simt < T_END:
	# Stroke kinematics (not force controlled yet)
	omega = 2 * np.pi * ctrl['freq']
	ph = omega * sim.simt
	th0 = ctrl['ampl'] * (np.sin(ph) + ctrl['strokedev'])
	dth0 = omega * ctrl['ampl'] * np.cos(ph)

	# POSITION_CONTROL uses Kp, Kd in [0,1]
	if STROKE_FORCE_CONTROL:
		tau = np.full(2, ctrl['thrust'] * (np.sin(ph) + ctrl['strokedev']))
		tau[0] += -urdfParams['stiffnessStroke'] * q[0] - urdfParams['dampingStroke'] * dq[0]
		tau[1] += -urdfParams['stiffnessStroke'] * q[2] - urdfParams['dampingStroke'] * dq[2]
		sim.setJointArray(bid, [0,2], sim.TORQUE_CONTROL, forces=tau)
	else:
		sim.setJointArray(bid, [0,2], sim.POSITION_CONTROL, targetPositions=[th0,th0], positionGains=[1,1], velocityGains=[0.1,0.1], forces=np.full(2, 1000000))

	# actual sim
	sim.sampleStates(bid)
	
	# Conventional controller
	# posErr = sim.q[4:7] - traj(sim.simt)
	# FIXME: simple path
	xDes = 0.03 * np.sin(2*np.pi*(sim.q[6] / 0.5))
	zdDes = 0.05
	ctrl['ampl'] = 0.9 - 0.1 * (sim.dq[6] - zdDes)
	desPitch = np.clip(-(200 * (sim.q[4] - xDes) + 0.1 * sim.dq[4]), -np.pi/4.0, np.pi/4.0)
	curPitch = Rotation.from_quat(sim.q[7:11]).as_euler('xyz')[1]
	pitchCtrl = 2 * (curPitch - desPitch) + 0.1 * sim.dq[8]
	ctrl['strokedev'] = np.clip(pitchCtrl, -0.4, 0.4)

	# Calculate aerodynamics
	pcop1, Faero1, Taero1 = bee.aerodynamics(sim.q, sim.dq, -1, worldFrame=AERO_WORLD_FRAME)
	pcop2, Faero2, Taero2 = bee.aerodynamics(sim.q, sim.dq, 1, worldFrame=AERO_WORLD_FRAME)
	sim.update(bid, [jointId[b'lwing_hinge'], jointId[b'rwing_hinge']], [pcop1, pcop2], [Faero1, Faero2], [Taero1, Taero2], worldFrame=AERO_WORLD_FRAME)
	time.sleep(sim._slowDown * sim.TIMESTEP)
	
	if sim.simt - tLastPrint > 0.01:
		# print(sim.simt, posErr)
		tLastPrint = sim.simt
	
sim.disconnect()
