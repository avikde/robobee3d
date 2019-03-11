import time
import numpy as np
import sys
import SimInterface
import FlappingModels3D
from scipy.spatial.transform import Rotation

# Usage params
STROKE_FORCE_CONTROL = False # if false, use position control on the stroke

sim = SimInterface.PyBullet()
# load robot
startPos = [0,0,1]
startOrientation = Rotation.from_euler('x', 0)
bid = sim.loadURDF("urdf/sdab.xacro.urdf", startPos, startOrientation.as_quat(), useFixedBase=False)
jointId, urdfParams = sim.getInfoFromURDF(bid)
print(jointId, urdfParams)
# Our model for control
bee = FlappingModels3D.QuasiSteadySDAB(urdfParams)

# Helper function: traj to track
def traj(t):
	return startPos + np.array([0.1 * np.sin(1 * t), 0, 0.5 * t])
# # draw traj
# p.addUserDebugLine(traj(tLastDraw), traj(simt), lineColorRGB=[0,0,1], lifeTime=0)
	
# ---

# Get the joints to reasonable positions
sim.resetJoints(bid, range(4), np.zeros(4), np.zeros(4))
# Passive hinge dynamics implemented as position control rather than joint dynamics
sim.setJointArray(bid, [1,3], sim.POSITION_CONTROL, targetPositions=[0,0], positionGains=urdfParams['stiffnessHinge']*np.ones(2), velocityGains=urdfParams['dampingHinge']*np.ones(2))

while sim.simt < 1:
	# No dynamics: reset positions
	omega = 2 * np.pi * 170.0
	ph = omega * sim.simt
	th0 = 1.0 * np.sin(ph)
	dth0 = omega * np.cos(ph)

	# POSITION_CONTROL uses Kp, Kd in [0,1]
	if STROKE_FORCE_CONTROL:
		tau = np.full(2, 4e-5 * np.sin(ph))
		tau[0] += -urdfParams['stiffnessStroke'] * q[0] - urdfParams['dampingStroke'] * dq[0]
		tau[1] += -urdfParams['stiffnessStroke'] * q[2] - urdfParams['dampingStroke'] * dq[2]
		sim.setJointArray(bid, [0,2], sim.TORQUE_CONTROL, forces=tau)
	else:
	  sim.setJointArray(bid, [0,2], sim.POSITION_CONTROL, targetPositions=[th0,th0], positionGains=[1,1], velocityGains=[0.1,0.1], forces=np.full(2, 1000000))

	# actual sim
	sim.sampleStates(bid)

	pcop1, Faero1 = bee.aerodynamics(sim.q, sim.dq, -1)
	pcop2, Faero2 = bee.aerodynamics(sim.q, sim.dq, 1)
	sim.update(bid, [jointId[b'lwing_hinge'], jointId[b'rwing_hinge']], [np.zeros(3), np.zeros(3)], [Faero1, Faero2])
	time.sleep(sim.SLOWDOWN * sim.TIMESTEP)
	
sim.disconnect()
