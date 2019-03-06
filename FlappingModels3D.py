import numpy as np
# trying to not rely on pybullet
from pyquaternion import Quaternion

class QuasiSteadySDAB:
	# Model with force control of the wing spar
	# Could make it modular so that lumped actuator models can be introduced as well

	def getCOP(self, q, bRight):
		# pass in current positions 11DOF: joints (4) + com (3 linear + 4 angular)
		pcop = q[4:7]
		quat = Quaternion(q[10], q[7], q[8], q[9])

		pcopB = np.array([0,0,0.1])

		return pcop + quat.rotate(pcopB)

	def aerodynamics(self, q, dq, bRight):
		pwingB = 0
		# pwing = pcom + 

		FaeroW = np.array([0, 0, 0])

		return self.getCOP(q, bRight), FaeroW
