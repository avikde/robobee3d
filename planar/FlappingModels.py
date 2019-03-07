import numpy as np

class PlanarThrustStrokeDev:
	mb = 100e-6
	l = 12e-3
	g = 9.81
	ib = 1/12. * mb * l**2
	d = 2e-3

	def getLin(self, y, u, dt):
		phi = y[2]
		cphi0 = np.cos(phi)
		sphi0 = np.sin(phi)

		# derivatives of the continuous dynamics
		# dy2dotdq = np.array([
		# 	[0,0,(cphi0*self.g)],
		# 	[0,0,-(self.g*sphi0)],
		# 	[0,0,(-2*self.d*self.g*self.mb*np.cos(2*phi))/self.ib]
		# ])
		# dy2dotdu = np.array([
		# 	[(sphi0/self.mb),0],
		# 	[cphi0/self.mb,0],
		# 	[-((self.d*np.sin(2*phi))/self.ib),(self.g*self.mb*np.cos(2*phi))/self.ib]
		# ])
		dy2dotdq = np.array([
			[0,0,-(cphi0*self.g)],
			[0,0,-(self.g*sphi0)],
			[0,0,0]
		])
		dy2dotdu = np.array([
			[-(sphi0/self.mb),0],
			[cphi0/self.mb,0],
			[0,(self.g*self.mb)/self.ib]
		])
		nq = dy2dotdq.shape[0]

		Aupper = np.hstack([np.eye(nq), dt * np.eye(nq)])
		Alower = np.hstack([dt * dy2dotdq, np.eye(nq)])
		Ad = np.vstack([Aupper, Alower])
		Bd = np.vstack([np.zeros((3,2)), dy2dotdu])

		return Ad, Bd
	
	def getLimits(self):
		umin = np.array([-self.mb*self.g, -50e-3])
		umax = np.array([3*self.mb*self.g, 50e-3])
		xmin = np.array([-np.inf,-np.inf,-10*np.pi,-np.inf,-np.inf,-np.inf])
		xmax = -xmin
		return umin, umax, xmin, xmax

	nx = 6
	nu = 2
