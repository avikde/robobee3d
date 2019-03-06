import numpy as np

class QuasiSteadySDAB:
	# Model with force control of the wing spar
	# Could make it modular so that lumped actuator models can be introduced as well

	def aerodynamics(bRight):
		pwingB = 0
		# pwing = pcom + 

		pcopW = np.array([0, 0, 0])
		FaeroW = np.array([0, 0, 0])

		return pcopW, FaeroW
