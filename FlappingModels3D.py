import numpy as np
# trying to not rely on pybullet
from scipy.spatial.transform import Rotation as R

class QuasiSteadySDAB:
	# Model with force control of the wing spar
	# Could make it modular so that lumped actuator models can be introduced as well

	# parameters
	d = 0.1
	Rwing = 0.07

	def getCOP(self, q, lrSign):
		# pass in current positions 11DOF: joints (4) + com (3 linear + 4 angular)
		pcop = q[4:7]
		quat = R.from_quat(q[7:11]) # scalar-last format

		# vector from center to distance along spar
		# joint angles
		theta = q[0:2]
		if lrSign > 0:
			theta = -q[2:4]

		pcopB = np.array([0,0,self.d]) + R.from_euler('z', theta[0]).apply(np.array([0, lrSign * self.Rwing, 0]))
		# Body to world frame
		return pcop + quat.apply(pcopB)

	def aerodynamics(self, q, dq, lrSign):
		pwingB = 0
		# pwing = pcom + 

		FaeroW = np.array([0, 0, 0])

		return self.getCOP(q, lrSign), FaeroW
