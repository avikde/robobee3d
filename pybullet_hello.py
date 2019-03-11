import time
import numpy as np
import sys
from scipy.spatial.transform import Rotation
import SimInterface
import FlappingModels3D

# Usage params
STROKE_FORCE_CONTROL = False # if false, use position control on the stroke
T_END = 1

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
	return startPos + np.array([0.05 * np.sin(2*np.pi*t), 0, 0.2 * t])

# draw traj
tdraw = np.linspace(0, T_END, 20)
for ti in range(1, len(tdraw)):
	sim.addUserDebugLine(traj(tdraw[ti-1]), traj(tdraw[ti]), lineColorRGB=[0,0,0], lifeTime=0)
	
# ---

# Get the joints to reasonable positions
sim.resetJoints(bid, range(4), np.zeros(4), np.zeros(4))
# Passive hinge dynamics implemented as position control rather than joint dynamics
sim.setJointArray(bid, [1,3], sim.POSITION_CONTROL, targetPositions=[0,0], positionGains=urdfParams['stiffnessHinge']*np.ones(2), velocityGains=urdfParams['dampingHinge']*np.ones(2))

while sim.simt < T_END:
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
	sim.update(bid, [jointId[b'lwing_hinge'], jointId[b'rwing_hinge']], [pcop1, pcop2], [Faero1, Faero2])
	time.sleep(sim.SLOWDOWN * sim.TIMESTEP)
	
sim.disconnect()
