import numpy as np
# trying to not rely on pybullet
from scipy.spatial.transform import Rotation

class QuasiSteadySDAB:
	# Model with force control of the wing spar
	# Could make it modular so that lumped actuator models can be introduced as well

	# Good parameters
	CD0 = 0.4
	CDmax = 3.4
	CLmax = 1.8
	# FIXME: parameters need to be updated
	d = 0.1
	ycp = 0.07
	
	def CF(self, a):
		# in order lift,drag
		return np.array([self.CLmax * np.sin(2*a), (self.CDmax + self.CD0) / 2.0 - (self.CDmax - self.CD0) / 2.0 * np.cos(2*a)])

	def aerodynamics(self, q, dq, lrSign):
		# pass in current positions 11DOF: joints (4) + com (3 linear + 4 angular)

		# vector from center to distance along spar
		# joint angles
		if lrSign > 0:
			theta = -q[2:4]
			dtheta = -dq[2:4]
		else:
			theta = q[0:2]
			dtheta = dq[0:2]

		# These are all in the body frame
		Rspar = Rotation.from_euler('z', theta[0])
		Rhinge = Rotation.from_euler('y', theta[1])
		# vector along wing
		sparVecB = Rspar.apply(np.array([0, lrSign, 0]))
		# wing chord unit vector
		chordB = Rspar.apply(Rhinge.apply(np.array([0,0,-1])))
		# TODO: add external wind vector
		# NOTE: assuming spar velocity is dominant
		wB = -np.cross(dtheta[0] * np.array([0,0,1]), self.ycp * sparVecB)
		# print(dtheta[0] * np.array([0,0,1]), self.ycp * sparVecB)

		# COP
		pcopB = np.array([0,0,self.d]) + self.ycp * sparVecB

		# Various directions in the notation of Osborne (1951)
		# l = vector along wing
		# w = relative wind vector
		# c = chord
		lwB = np.cross(sparVecB, wB)
		lwpB = wB - wB.dot(sparVecB) * sparVecB
		wnorm = np.sqrt(wB.dot(wB))
		lwnorm = np.sqrt(lwB.dot(lwB))
		lwpnorm = np.sqrt(lwpB.dot(lwpB))

		# Lift/drag directions
		eD = lwpB / lwpnorm
		eL = np.array([0,0,1]) #lwB / lwnorm
		# FIXME: the latter version needs some reversal for one half-stroke
		# the uncommented version agrees with Chen (2017) science robotics

		# Calculate aero force
		aoa = np.arccos(chordB.dot(wB) / wnorm)
		Cf = self.CF(aoa)
		# Cf *= 0.5 * rho * beta
		FaeroB = (Cf[0] * eL + Cf[1] * eD) * lwnorm**2

		# Body to world frame --
		pcom = q[4:7]
		Rb = Rotation.from_quat(q[7:11]) # scalar-last format
		pcopW = pcom + Rb.apply(pcopB)
		FaeroW = Rb.apply(FaeroB)

		return pcopW, FaeroW
