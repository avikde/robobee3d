import numpy as np
import sys
sys.path.append('..')
import controlutils.py.kinematics as kin
import controlutils.py.misc as misc

class PlanarThrustStrokeDev:
	mb = 100e-6
	l = 12e-3  # body length
	w = 2e-3  # body width (for visualization only)
	g = 9.81
	ib = 1/12. * mb * l**2
	d = 2e-3
	# initial conditions
	nx = 6  # originally 6 states for linearized (affine) model, but adding on another 6 to store the affine pieces including the gravity terms as well as fk the linearization residual
	nu = 2
	y0 = np.array([0,0,0,0,0,0])
	u0 = np.zeros(nu)
	# can be reset
	dt = 0.005
	STROKE_EXTENT = 3e-3

	def getLinearDynamics(self, y, u):
		'''Returns Ad, Bd[, fd]
		where fd is only there if the system is affine.
		x[k+1] = Ad @ x[k] + Bd @ u[k] + fd
		'''
		# Need the actual cts dynamics for the affine term
		fy = self.dydt(y, u)  # result is 6x1

		# Linearizations of the continuous dynamics. TODO: In the future can use autograd
		phi = y[2]
		cphi0 = np.cos(phi)
		sphi0 = np.sin(phi)

		thrust = self.g + u[0]/self.mb

		All = np.array([
			[0,0,-thrust * cphi0],
			[0,0,-thrust * sphi0],
			[0,0,0]
		])
		Bl = np.array([
			[-(sphi0/self.mb),0],
			[cphi0/self.mb,0],
			[u[1]/self.ib,thrust * self.mb / self.ib]
		])
		# Jac wrt state
		dfdy = np.vstack((
			np.hstack((np.zeros((3,3)), np.eye(3))),
			np.hstack((All, np.zeros((3,3))))
		))
		# Jac wrt inputs
		dfdu = np.vstack((
			np.zeros((3,2)),
			Bl
		))

		# Compute the affine term
		fd = self.dt * (fy - dfdy @ y[0:6] - dfdu @ u)

		# New linearization: see https://github.com/avikde/robobee3d/pull/29#issuecomment-475222003
		Ad = np.eye(6) + self.dt * dfdy
		Bd = self.dt * dfdu

		return Ad, Bd, fd
	
	def getLimits(self):
		umin = np.array([-self.mb*self.g, -self.STROKE_EXTENT])
		umax = np.array([3*self.mb*self.g, self.STROKE_EXTENT])
		# umin = np.array([-np.inf, -50e-3])
		# umax = np.array([np.inf, 50e-3])
		xmin = np.array([-np.inf,-np.inf,-10*np.pi,-np.inf,-np.inf,-np.inf])
		xmax = -xmin
		return umin, umax, xmin, xmax

	def dydt(self, y, u):
		# Full continuous nonlinear vector field
		phi = y[2]
		thrust = self.g + u[0]/self.mb
		# accelerations
		y2dot = np.array([-thrust * np.sin(phi), 
		-self.g + thrust * np.cos(phi), 
		thrust * u[1] * self.mb / self.ib
		])
		y1dot = y[3:6]
		return np.hstack((y1dot, y2dot))


	def dynamics(self, y, u, useLinearization=False):
		umin, umax, _, _ = self.getLimits()
		# FIXME: input constraints are not being satisfied. See #36
		# u = np.clip(u, umin, umax)
		u[1] = np.clip(u[1], umin[1], umax[1])
		# Full nonlinear dynamics
		if useLinearization:
			Ad, Bd, fd = self.getLinearDynamics(y, u)
			# print(Ad)
			return Ad @ y + Bd @ u + fd
		else:
			# 1st order integration
			return y[0:6] + self.dydt(y, u) * self.dt # np.hstack((yNog, self.g))

	# Non-standard model functions
	def visualizationInfo(self, y, u, ax, col='r', Faeroscale=1, rawxy=False):
		Ryaw = kin.rot2(y[2])
		pcop = y[0:2] + Ryaw @ np.array([u[1],self.d])
		Faero = Ryaw @ np.array([0, self.mb * self.g + u[0]])
		strokeExtents = np.vstack((y[0:2] + Ryaw @ np.array([-2*self.STROKE_EXTENT, self.d]), y[0:2] + Ryaw @ np.array([2*self.STROKE_EXTENT, self.d])))
		# Plot these custom things from here
		ax.plot(strokeExtents[:,0], strokeExtents[:,1], 'k--', linewidth=1,  alpha=0.3)
		ax.arrow(pcop[0], pcop[1], Faero[0], Faero[1], width=0.0002, alpha=0.3, facecolor=col)
		return misc.rectangle(y[0:2], y[2], self.w, self.l, rawxy)


class PlanarStrokeSpeed:
	# also trying rescaling the problem
	mb = 100e-6
	l = 12e-3  # body length
	w = 2e-3  # body width (for visualization only)
	g = 9.81
	ib = 1/12. * mb * l**2
	d = 2e-3

	nx = 7
	nu = 2
	y0 = np.array([0,0,0,0,0,0,0])

		# TODO:
	kat = 1e-3
	omegat = 1

	def getLinearDynamics(self, y, u):
		'''Returns Ad, Bd[, fd]
		where fd is only there if the system is affine.
		x[k+1] = Ad @ x[k] + Bd @ u[k] + fd
		'''
		raise NotImplementedError
	
	def getLimits(self):
		# return umin, umax, xmin, xmax
		raise NotImplementedError

	def dydt(self, y, u, avg=False):
		''' Full continuous nonlinear vector field when avg=False
		Otherwise return dydsigma
		'''
		phi = y[2]
		uss = u[0]
		v = y[-1]
		# v'(psi) is the stroke speed
		dv = u[0]
		sdv = np.sign(dv)
		cphi = np.cos(phi)
		sphi = np.sin(phi)

		if avg:
			dphi = y[5]
			u1 = u[0]
			psi0 = u[1]  # reversal point
			sigma0 = y[-1]  #initial stroke

			# from mathematica. for y = (phi,dx,dz,dphi)
			fav = np.array([dphi/self.omegat, 
			-self.kat * u1**2 * psi0 * self.omegat * ((1-2*psi0) * cphi + sphi) / (self.mb * (-1. + psi0)),
			-self.g / self.omegat + self.kat * u1**2 * psi0 * self.omegat * (-cphi + (1-2*psi0) * sphi) / (self.mb * (-1. + psi0)),
			-(self.kat * u1**2 * (2 * sigma0 * (1 + (-1+psi0)*psi0) + psi0*(u1 + self.d*(2-4*psi0) + u1*(-1+psi0)*psi0)) * self.omegat)/(2*self.ib*(-1+psi0))
			])
			
			# return 
			dxzphi = y[3:6]/self.omegat
			ddxzphi = fav[1:]
			return np.hstack((dxzphi, ddxzphi, dv / self.omegat))
		else:
			ddxzphi = np.array([dv**2 * self.kat * self.omegat * (cphi * sdv + sphi) / self.mb, 
			-self.g / self.omegat + dv**2 * self.kat * self.omegat * (cphi - sphi * sdv) / self.mb,
			dv**2 * self.kat * self.omegat * (v + self.d * sdv) / self.ib])
			# return np.hstack((y1dot, y2dot))

			return np.hstack((y[3:6], ddxzphi, dv))


	def dynamics(self, y, u, useLinearization=False):
		# Full nonlinear dynamics
		if useLinearization:
			Ad, Bd, fd = self.getLinearDynamics(y, u)
			# print(Ad)
			return Ad @ y + Bd @ u + fd
		else:
			# 1st order integration
			# FIXME: need to figure out time to stroke end
			return y + self.dydt(y, u) * dt

	# Non-standard model functions
	def visualizationInfo(self, y, u, ax, col='r', Faeroscale=1, rawxy=False):
		# Ryaw = kin.rot2(y[2])
		# pcop = y[0:2] + Ryaw @ np.array([u[1],self.d])
		# Faero = Ryaw @ np.array([0, self.mb * self.g + u[0]])
		# strokeExtents = np.vstack((y[0:2] + Ryaw @ np.array([-2*self.STROKE_EXTENT, self.d]), y[0:2] + Ryaw @ np.array([2*self.STROKE_EXTENT, self.d])))
		# # Plot these custom things from here
		# ax.plot(strokeExtents[:,0], strokeExtents[:,1], 'k--', linewidth=1,  alpha=0.3)
		# ax.arrow(pcop[0], pcop[1], Faero[0], Faero[1], width=0.0002, alpha=0.3, facecolor=col)
		return misc.rectangle(y[0:2], y[2], self.w, self.l, rawxy)
	

def visualizeTraj(ax, traj, model, col='r', xylim=None):
	# Plots what RefTraj.generate() returns
	from matplotlib.patches import Rectangle, Circle
	from matplotlib.collections import PatchCollection
	import controlutils.py.misc as misc
	
	N = traj['q'].shape[0]
	robotBodies = []
	for k in range(N):
		qk = traj['q'][k,:]
		uk = traj['u'][k,:] if traj['u'] is not None else None

		# get info from model
		body = model.visualizationInfo(qk, uk, ax, col)
		
		robotBodies.append(body)
	
	pc = PatchCollection(robotBodies, facecolor=col, edgecolor='k', alpha=0.3)
	ax.add_collection(pc)

	ax.set_aspect(1)
	if xylim is not None:
		ax.set_xlim(xylim[0:2])
		ax.set_ylim(xylim[2:4])
	ax.grid(True)
	ax.set_xlabel('x')
	ax.set_ylabel('z')
	
if __name__ == "__main__":
	import matplotlib.pyplot as plt
	import controlutils.py.lqr as lqr
	np.set_printoptions(suppress=True, linewidth=100, precision=3)
	
	model = PlanarThrustStrokeDev()
	model.dt = 0.01

	# For visualization
	fig, ax = plt.subplots(1)
	Ndraw = 100

	Y = np.zeros((Ndraw, model.nx))
	U = np.zeros((Ndraw, model.nu))
	# initial conditions
	Y[0,:], U[0,:] = model.y0, model.u0
	# Openloop
	U[0,:] = np.zeros(2)
	lqrgoal = np.array([-0.03, -0.02, 0, 0, 0, 0])
	for ti in range(1, Ndraw):
		# get linearization
		Ad, Bd, fd = model.getLinearDynamics(Y[ti-1,:], U[ti-1,:])
		# actually compute u
		K, X = lqr.dlqr(Ad, Bd, np.diag([100,100,10,1,1,1]), np.diag([0.1,0.1]))
		ulqr = K @ (lqrgoal - Y[ti-1,:])
		Y[ti,:] = Ad @ Y[ti-1,:] + Bd @ ulqr + fd  # actual dynamics
		U[ti,:] = ulqr  # for next linearization
	visualizeTraj(ax, {'q':Y[:, 0:3], 'u':U}, model)
	ax.plot(lqrgoal[0], lqrgoal[1], 'c*')
	print(Y)

	Ndraw = 5
	Y = np.zeros((Ndraw, model.nx))
	U = np.zeros((Ndraw, model.nu))
	# initial conditions
	Y[0,:], U[0,:] = model.y0, model.u0
	# Openloop
	Y[0,0] = 0.02
	U[0,:] = np.array([1e-3,1e-3])
	for ti in range(1, Ndraw):
		Y[ti,:] = model.dynamics(Y[ti-1,:], U[ti-1,:], useLinearization=True)
		U[ti,:] = U[ti-1,:]
	visualizeTraj(ax, {'q':Y[:, 0:3], 'u':U}, model, col='b')
	print(Y)

	# custom legend
	from matplotlib.lines import Line2D
	custom_lines = [Line2D([0], [0], color='r', alpha=0.3), 
		Line2D([0], [0], color='b', alpha=0.3)]
	ax.legend(custom_lines, ['LQR', 'Linearized'])

	plt.show()