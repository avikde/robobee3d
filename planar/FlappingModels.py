import numpy as np
import controlutils.py.ControlUtils as cu

class RefTraj:
	# generate a reference trajectory
	def __init__(self, N, dt):
		self.N = N
		self.dt = dt

	def generate(self, q0, twistDes):
		qtraj = cu.twistInt(q0, twistDes, self.N, self.dt)
		return {'q':qtraj}

class PlanarThrustStrokeDev:
	mb = 100e-6
	l = 12e-3
	g = 9.81
	ib = 1/12. * mb * l**2
	d = 2e-3
	# initial conditions
	y0 = np.array([0,0,0,0,0,0,g])
	u0 = np.array([0, 0])

	def getLinearDynamics(self, y, u):
		phi = y[2]
		cphi0 = np.cos(phi)
		sphi0 = np.sin(phi)

		# derivatives of the continuous dynamics
		# New: adding g as a state
		All = self.dt / self.mb * np.array([
			[0,0,-(self.g * self.mb + u[0]) * cphi0/self.mb],
			[0,0,-(self.g + u[0] / self.mb) * sphi0],
			[0,0,0]
		])
		Bl = self.dt * np.array([
			[-(sphi0/self.mb),0],
			[cphi0/self.mb,0],
			[u[1]/self.ib,(self.g * self.mb + u[0])/self.ib]
		])
		nq = All.shape[0]

		Aupper = np.hstack([np.eye(nq), self.dt * np.eye(nq), np.zeros((nq,1))])
		Alower = np.hstack([All, np.eye(nq), self.dt * np.array([[-sphi0],[-1 + cphi0],[self.mb * u[1]/self.ib]])])
		Ag = np.array([0,0,0,0,0,0,1])
		Ad = np.vstack([Aupper, Alower, Ag])
		Bd = np.vstack([np.zeros((nq,2)), Bl, np.zeros((1,2))])

		return Ad, Bd
	
	def getLimits(self):
		umin = np.array([-self.mb*self.g, -5e-3])
		umax = np.array([3*self.mb*self.g, 5e-3])
		# umin = np.array([-np.inf, -50e-3])
		# umax = np.array([np.inf, 50e-3])
		xmin = np.array([-np.inf,-np.inf,-10*np.pi,-np.inf,-np.inf,-np.inf,-np.inf])
		xmax = -xmin
		return umin, umax, xmin, xmax

	def dynamics(self, y, u, useLinearization=False):
		# Full nonlinear dynamics
		if useLinearization:
			Ad, Bd = self.getLinearDynamics(y, u)
			return Ad @ y + Bd @ u
		else:
			phi = y[2]
			thrust = self.g * self.mb + u[0]
			# accelerations
			y2dot = np.array([-thrust * np.sin(phi), 
			-self.g * self.mb + thrust * np.cos(phi), 
			thrust * u[1] / self.ib
			])
			y1dot = y[3:6]
			# zero order hold?
			yNog = y[0:6] + np.hstack((y1dot, y2dot)) * self.dt
			return np.hstack((yNog, self.g))

	nx = 7
	nu = 2
