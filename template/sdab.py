import time
import numpy as np
import pybullet as p
from robobee import RobobeeSim
np.set_printoptions(precision=2, suppress=True, linewidth=200)

# Usage params
STROKE_FORCE_CONTROL = False # if false, use position control on the stroke
T_END = 1000
AERO_WORLD_FRAME = True

bee = RobobeeSim(slowDown=0.01, camLock=True, timestep=0.1)
# load robot
startPos = [0,0,10]
startOrientation = p.getQuaternionFromEuler(np.zeros(3))
bid = bee.load("../urdf/sdab.xacro.urdf", startPos, startOrientation, useFixedBase=False)
# print(jointId, urdfParams)

# Helper function: traj to track
def traj(t):
	return startPos + np.array([30 * np.sin(0.002*np.pi*t), 0, 0.5 * t])

# draw traj
tdraw = np.linspace(0, T_END, 20)
for ti in range(1, len(tdraw)):
	p.addUserDebugLine(traj(tdraw[ti-1]), traj(tdraw[ti]), lineColorRGB=[0,0,0], lifeTime=0)
	
# ---


# conventional controller params
ctrl = {'thrust': 4e-5, 'strokedev': 0, 'ampl': np.pi/4, 'freq': 0.17}
tLastPrint = 0

while True:
	try:
		# Stroke kinematics (not force controlled yet)
		omega = 2 * np.pi * ctrl['freq']
		ph = omega * bee.simt
		th0 = ctrl['ampl'] * (np.sin(ph) + ctrl['strokedev'])
		dth0 = omega * ctrl['ampl'] * np.cos(ph)

		# POSITION_CONTROL uses Kp, Kd in [0,1]
		if STROKE_FORCE_CONTROL:
			# tau = np.full(2, ctrl['thrust'] * (np.sin(ph) + ctrl['strokedev']))
			# tau[0] += -urdfParams['stiffnessStroke'] * q[0]
			# tau[1] += -urdfParams['stiffnessStroke'] * q[2]
			tau = [0,0]
			p.setJointMotorControlArray(bid, [0,2], p.TORQUE_CONTROL, forces=tau)
		else:
			p.setJointMotorControlArray(bid, [0,2], p.POSITION_CONTROL, targetPositions=[th0,th0], positionGains=[1,1], velocityGains=[1,1], forces=np.full(2, 1000000))

		# actual sim
		bee.sampleStates(bid)
		
		# # Conventional controller
		# # posErr = sim.q[4:7] - traj(sim.simt)
		# # FIXME: simple path
		# xDes = 0.03 * np.sin(2*np.pi*(sim.q[6] / 0.5))
		# zdDes = 0.05
		# ctrl['ampl'] = 0.9 - 0.1 * (sim.dq[6] - zdDes)
		# desPitch = np.clip(-(200 * (sim.q[4] - xDes) + 0.1 * sim.dq[4]), -np.pi/4.0, np.pi/4.0)
		# curPitch = Rotation.from_quat(sim.q[7:11]).as_euler('xyz')[1]
		# pitchCtrl = 2 * (curPitch - desPitch) + 0.1 * sim.dq[8]
		# ctrl['strokedev'] = np.clip(pitchCtrl, -0.4, 0.4)

		# Calculate aerodynamics
		pcop1, Faero1, Taero1 = bee.aerodynamics(bee.q, bee.dq, -1, worldFrame=AERO_WORLD_FRAME)
		pcop2, Faero2, Taero2 = bee.aerodynamics(bee.q, bee.dq, 1, worldFrame=AERO_WORLD_FRAME)
		bee.update(bid, [pcop1, pcop2], [Faero1, Faero2], [Taero1, Taero2], worldFrame=AERO_WORLD_FRAME)
		time.sleep(bee._slowDown * bee.TIMESTEP)
	except:
		raise

