import numpy as np
# trying to not rely on pybullet
from scipy.spatial.transform import Rotation as R

class QuasiSteadySDAB:
	# Model with force control of the wing spar
	# Could make it modular so that lumped actuator models can be introduced as well

	# parameters
	d = 0.1

	def getCOP(self, q, bRight):
		# pass in current positions 11DOF: joints (4) + com (3 linear + 4 angular)
		pcop = q[4:7]
		quat = R.from_quat(q[7:11]) # scalar-last format

		pcopB = np.array([0,0,self.d])
		# Body to world frame
		return pcop + quat.apply(pcopB)

	def aerodynamics(self, q, dq, bRight):
		pwingB = 0
		# pwing = pcom + 

		FaeroW = np.array([0, 0, 0])

		return self.getCOP(q, bRight), FaeroW
