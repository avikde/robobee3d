import numpy as np
import sys
sys.path.append('..')
import controlutils.py.geometry as geom

class RefTraj:
	# generate a reference trajectory
	def __init__(self, N, dt):
		self.N = N
		self.dt = dt

	def generate(self, q0, twistDes):
		qtraj = geom.twistInt(q0, twistDes, self.N, self.dt)
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
	# can be reset
	dt = 0.005

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

def visualizeTraj(ax, traj):
	# Plots what RefTraj.generate() returns
	from matplotlib.patches import Rectangle, Circle
	from matplotlib.collections import PatchCollection
	import controlutils.py.misc as misc
	
	N = traj['q'].shape[0]
	robotBodies = []
	toes = []
	for k in range(N):
		qk = traj['q'][k,:]
		robotBodies.append(misc.rectangle(qk[0:2], qk[2], 0.002, 0.005))

		# # wing
		# if k > 0:
		# 	pk = traj['p'][k,:]
		# 	for j in range(2):
		# 		pjk = pk[2*j:2*j+2]
		# 		toes.append(Circle(pjk, radius=toeRadius))
		# 		# add some text so we can tell k
		# 		ax.text(pjk[0], pjk[1], str(k))
	
	pc = PatchCollection(robotBodies, facecolor='r', edgecolor='k', alpha=0.3)
	ax.add_collection(pc)
	pc = PatchCollection(toes, facecolor='b', edgecolor='k', alpha=0.3)
	ax.add_collection(pc)

	ax.set_aspect(1)
	ax.set_xlim([-0.01,0.01])
	ax.set_ylim([-0.01,0.01])
	ax.grid(True)
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	
if __name__ == "__main__":
	import matplotlib.pyplot as plt
	# For visulatization
	fig, ax = plt.subplots(2)
	# TODO: plot traj
	model = PlanarThrustStrokeDev()
	Ndraw = 10
	Y = np.zeros((Ndraw, model.nx))
	U = np.zeros((Ndraw, model.nu))
	# initial conditions
	Y[0,:], U[0,:] = model.y0, model.u0
	U[0,:] = np.array([1e-3, 1e-5])
	for ti in range(1, Ndraw):
		Y[ti,:] = model.dynamics(Y[ti-1,:], U[ti-1,:])
		U[ti,:] = U[ti-1,:]
	visualizeTraj(ax[0], {'q':Y[:, 0:3]})

	plt.show()