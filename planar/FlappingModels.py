import numpy as np
import sys
sys.path.append('..')
import controlutils.py.geometry as geom
import controlutils.py.misc as misc

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
	l = 12e-3  # body length
	w = 2e-3  # body width (for visualization only)
	g = 9.81
	ib = 1/12. * mb * l**2
	d = 2e-3
	# initial conditions
	y0 = np.array([0,0,0,0,0,0,g])
	u0 = np.array([0, 0])
	# can be reset
	dt = 0.005
	STROKE_EXTENT = 5e-3

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
		umin = np.array([-self.mb*self.g, -self.STROKE_EXTENT])
		umax = np.array([3*self.mb*self.g, self.STROKE_EXTENT])
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

	# Non-standard model functions
	def visualizationInfo(self, y, u, Faeroscale=0.02):
		Ryaw = geom.rot2(y[2])
		pcop = y[0:2] + Ryaw @ (np.array([0,self.d]) + u[1])
		Faero = Ryaw @ np.array([0, self.mb * self.g + u[0]])
		strokeExtents = y[0:2] + np.vstack((Ryaw @ np.array([-self.STROKE_EXTENT, self.d]), Ryaw @ np.array([self.STROKE_EXTENT, self.d])))
		return misc.rectangle(y[0:2], y[2], self.w, self.l), pcop, Faeroscale * Faero, strokeExtents

	nx = 7
	nu = 2

def visualizeTraj(ax, traj, model):
	# Plots what RefTraj.generate() returns
	from matplotlib.patches import Rectangle, Circle
	from matplotlib.collections import PatchCollection
	import controlutils.py.misc as misc
	
	N = traj['q'].shape[0]
	robotBodies = []
	for k in range(N):
		qk = traj['q'][k,:]
		uk = traj['u'][k,:]

		# get info from model
		body, pcop, Faero, strokeExtents = model.visualizationInfo(qk, uk)
		
		robotBodies.append(body)
		ax.plot(strokeExtents[:,0], strokeExtents[:,1], 'k--', linewidth=1,  alpha=0.3)
		ax.arrow(pcop[0], pcop[1], Faero[0], Faero[1], width=0.0003, alpha=0.3)
	
	pc = PatchCollection(robotBodies, facecolor='r', edgecolor='k', alpha=0.3)
	ax.add_collection(pc)

	ax.set_aspect(1)
	ax.set_xlim([-0.03,0.03])
	ax.set_ylim([-0.03,0.03])
	ax.grid(True)
	ax.set_xlabel('x')
	ax.set_ylabel('y')
	
if __name__ == "__main__":
	import matplotlib.pyplot as plt
	# For visulatization
	fig, ax = plt.subplots(1)
	
	model = PlanarThrustStrokeDev()
	model.dt = 0.02
	Ndraw = 5
	Y = np.zeros((Ndraw, model.nx))
	U = np.zeros((Ndraw, model.nu))
	# initial conditions
	Y[0,:], U[0,:] = model.y0, model.u0
	# FIXME: these are just openloop inputs to test vis
	U[0,:] = np.array([1e-1, 1e-6])
	for ti in range(1, Ndraw):
		Y[ti,:] = model.dynamics(Y[ti-1,:], U[ti-1,:], useLinearization=False)
		U[ti,:] = U[ti-1,:]
	visualizeTraj(ax, {'q':Y[:, 0:3], 'u':U}, model)

	plt.show()