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
	
	# rescaled system
	lscaled = 0.1

	def __init__(self, rescale=False):
		self.rescale = rescale
		self.ibmb = 1./12. * self.lscaled**2 if rescale else self.ib / self.mb

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

		thrust = self.g + u[0] if self.rescale else self.g + u[0]/self.mb

		All = np.array([
			[0,0,-thrust * cphi0],
			[0,0,-thrust * sphi0],
			[0,0,0]
		])
		if self.rescale:
			# inputs are u0tilde or "thrust", and u1
			Bl = np.array([
				[-(sphi0),0],
				[cphi0,0],
				[u[1]/self.ibmb, thrust/self.ibmb]
			])
		else:
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

		# Compute the affine term. Look at notes, this is right.
		fd = self.dt * (fy - dfdy @ y[0:6] - dfdu @ u)

		# New linearization: see https://github.com/avikde/robobee3d/pull/29#issuecomment-475222003
		Ad = np.eye(6) + self.dt * dfdy
		Bd = self.dt * dfdu

		return Ad, Bd, fd
	
	def getLimits(self):
		if self.rescale:
			umin = np.array([-self.g, -0.05 * self.lscaled])
			umax = np.array([3*self.g, 0.05 * self.lscaled])
		else:
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
		thrust = self.g + u[0] if self.rescale else self.g + u[0]/self.mb
		
		# accelerations
		y2dot = np.array([-thrust * np.sin(phi), 
		-self.g + thrust * np.cos(phi), 
		thrust * u[1] / self.ibmb
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
	def visualizationInfo(self, y, u, ax=None, col='r', Faeroscale=1, rawxy=False):
		Ryaw = kin.rot2(y[2])
		pcop = y[0:2] + Ryaw @ np.array([u[1],self.d])
		
		if self.rescale:
			Faero = Ryaw @ np.array([0, self.g + u[0]])
		else:
			Faero = Ryaw @ np.array([0, self.mb * self.g + u[0]])
		umin, umax, _, _ = self.getLimits()
		strokeExtents = np.vstack((y[0:2] + Ryaw @ np.array([umin[1], self.d]), y[0:2] + Ryaw @ np.array([umax[1], self.d])))
		if ax is not None:
			# plot onto ax
			ax.plot(strokeExtents[:,0], strokeExtents[:,1], 'k--', linewidth=1,  alpha=0.3)
			Faero *= Faeroscale
			ax.arrow(pcop[0], pcop[1], Faero[0], Faero[1], width=0.0002, alpha=0.3, facecolor=col)

			return misc.rectangle(y[0:2], y[2], self.w, self.l, rawxy)
		else:
			# return all the stuff necessary to draw
			return misc.rectangle(y[0:2], y[2], self.w, self.l, rawxy), pcop, Faero, strokeExtents


class PlanarStrokeSpeed:
	# also trying rescaling the problem
	mb = 100e-6
	l = 12e-3  # body length
	w = 2e-3  # body width (for visualization only)
	maxChord = 5e-3  # wing chord for vis
	g = 9.81
	ib = 1/12. * mb * l**2
	d = 2e-3
	STROKE_EXTENT = 3e-3

	nx = 7
	nu = 2
	y0 = np.array([0,0,0,0,0,0,0])
	
	rescale = True

	# Parameter values
	eps = 1e-2
	ka = 1e-3
	
	# TODO:
	omega = 100
	# Relate params by eps
	kat = ka/eps**2
	omegat = omega * eps

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

		#

		if avg:
			dphi = y[5]
			u1 = u[0]
			sigma0 = self.sigma0  #initial stroke
			# FIXME: this needs to be sampled at the section
			psi0 = 0.42 # (u[1] - sigma0)/u[0]  # reversal point in phase
			# print(psi0)

			# 
			mb = self.mb
			ib = self.ib
			kat = self.kat
			omegat = self.omegat
			u12 = u1**2
			g = self.g
			d = self.d
			w1 = u[0] / self.omega
			w12 = w1**2

			# from mathematica. for (dx,dz,dphi)
			fav = np.array([
				-((kat*omegat*(cphi*(1 - 2*psi0) + (1 - 2*psi0 + 2*psi0**2)*sphi)*w12)/(mb*(-1 + psi0)**2)),
				-(g/omegat) + (cphi*kat*omegat*(1 + psi0**2/(-1 + psi0)**2)* w12)/mb + (kat*omegat*(-1 + 2*psi0)*sphi*w12)/(mb*(-1 + psi0)**2),
				(kat*omegat*(d*(-2 + 4*psi0) + (1 - 2*psi0 + 2*(psi0)**2)*(2*sigma0 + psi0*w1))*w12)/(2.*ib*(-1 + psi0)**2)
			])
			
			# return 
			dxzphi = y[3:6]
			# fav = avg(dy/dpsi), so dy/dt ~= fav * dpsi/dt
			ddxzphi = fav * omegat
			return np.hstack((dxzphi, ddxzphi, dv))
		else:
			ddxzphi = np.array([-dv**2 * self.ka * (cphi * sdv + sphi) / self.mb, 
			-self.g + dv**2 * self.ka * (cphi - sphi * sdv) / self.mb,
			dv**2 * self.ka * (v - self.d * sdv) / self.ib])
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
	def visualizationInfo(self, y, u, ax, col='r', rawxy=False, Faeroscale=2 * STROKE_EXTENT):
		Rb = kin.rot2(y[2])
		strokeExtents = np.vstack((y[0:2] + Rb @ np.array([-2*self.STROKE_EXTENT, self.d]), y[0:2] + Rb @ np.array([2*self.STROKE_EXTENT, self.d])))
		# Plot these custom things from here
		ax.plot(strokeExtents[:,0], strokeExtents[:,1], 'k--', linewidth=1,  alpha=0.3)

		# plt the wing
		Rhinge = kin.rot2(0.25*np.pi if u[0] < 0 else -0.25*np.pi)
		wingStartB = np.array([y[-1], self.d])
		wingEndB = wingStartB + Rhinge @ np.array([0, -self.maxChord])
		wingExtents = np.vstack((y[0:2] + Rb @ wingStartB, y[0:2] + Rb @ wingEndB))
		ax.plot(wingExtents[:,0], wingExtents[:,1], 'k-', linewidth=2,  alpha=0.5)
		# stroke vel arrow
		strokeVelDirection = Rb @ np.array([u[0], 0]) * Faeroscale
		midWing = np.mean(wingExtents, axis=0)
		ax.arrow(midWing[0], midWing[1], strokeVelDirection[0], strokeVelDirection[1], width=0.0002, alpha=0.5, facecolor=col)
		
		return misc.rectangle(y[0:2], y[2], self.w, self.l, rawxy)
	

def visualizeTraj(ax, traj, model, col='r', Faeroscale=1, tplot=None):
	'''Plots a trajectory, with model info from model.visualizationInfo'''
	from matplotlib.patches import Rectangle, Circle
	from matplotlib.collections import PatchCollection
	import controlutils.py.misc as misc

	Faeroscale = 1e-4 if model.rescale else 1
	
	robotBodies = []

	if tplot is None:
		# plot the entire trajectory
		trajp = traj
	else:
		# use tplot as the keyframe times to display
		assert 't' in traj
		trajp = {}
		for keyi in ['q','u']:
			if traj[keyi] is not None:
				trajfunc = interp1d(traj['t'], traj[keyi], axis=0)
				trajp[keyi] = trajfunc(tplot)
			else:
				trajp[keyi] = None 

	for k in range(trajp['q'].shape[0]):
		qk = trajp['q'][k,:]
		uk = trajp['u'][k,:] if trajp['u'] is not None else None

		# get info from model
		# body = model.visualizationInfo(qk, uk, ax, col)
		body = model.visualizationInfo(qk, uk, ax, Faeroscale=Faeroscale)
		
		robotBodies.append(body)
	
	pc = PatchCollection(robotBodies, facecolor=col, edgecolor='k', alpha=0.3)
	ax.add_collection(pc)

	ax.set_aspect(1)
	ax.set_xlim([np.amin(traj['q'][:,0])-0.05,np.amax(traj['q'][:,0])+0.05])
	ax.set_ylim([np.amin(traj['q'][:,1])-0.05,np.amax(traj['q'][:,1])+0.05])
	ax.grid(True)
	ax.set_xlabel('x')
	ax.set_ylabel('z')
	
if __name__ == "__main__":
	import matplotlib.pyplot as plt
	import controlutils.py.lqr as lqr
	np.set_printoptions(suppress=True, linewidth=100, precision=3)
	
	rescale = True

	model = PlanarThrustStrokeDev(rescale=rescale)
	model.dt = 0.1 if rescale else 0.01
	Faeroscale = 1e-5 if rescale else 1

	# For visualization
	fig, ax = plt.subplots(1)
	Ndraw = 10

	# LQR
	wx = np.diag([100,100,10,1,1,1])
	wu = np.diag([0.1,0.1])
	Y = np.zeros((Ndraw, model.nx))
	U = np.zeros((Ndraw, model.nu))
	# initial conditions
	Y[0,:], U[0,:] = model.y0, model.u0
	# Openloop
	U[0,:] = np.array([0.01, 0.0]) if rescale else np.zeros(2)
	lqrgoal = np.array([-0.03, -0.02, 0, 0, 0, 0])
	for ti in range(1, Ndraw):
		# get linearization
		Ad, Bd, fd = model.getLinearDynamics(Y[ti-1,:], U[ti-1,:])
		# actually compute u
		try:
			K, X = lqr.dlqr(Ad, Bd, wx, wu)
		except np.linalg.linalg.LinAlgError:
			print("Ad =", Ad, "Bd =", Bd, "fd =", fd)
			raise
		ulqr = K @ (lqrgoal - Y[ti-1,:])
		Y[ti,:] = Ad @ Y[ti-1,:] + Bd @ ulqr + fd  # actual dynamics
		U[ti,:] = ulqr  # for next linearization
	visualizeTraj(ax, {'q':Y[:, 0:3], 'u':U}, model)
	ax.plot(lqrgoal[0], lqrgoal[1], 'c*')
	print(Y)

	Ndraw = 5
	Y = np.zeros((Ndraw, model.nx))
	Ynl = np.zeros_like(Y)
	U = np.zeros((Ndraw, model.nu))
	# initial conditions
	Y[0,:], U[0,:] = model.y0, model.u0
	# Openloop
	Y[0,0] = 0.02
	Ynl[0,0] = 0.03
	U[0,:] = np.array([0.1, 0.05 * model.lscaled]) if rescale else np.array([1e-3,1e-3])
	for ti in range(1, Ndraw):
		Y[ti,:] = model.dynamics(Y[ti-1,:], U[ti-1,:], useLinearization=True)
		Ynl[ti,:] = model.dynamics(Ynl[ti-1,:], U[ti-1,:], useLinearization=False)
		U[ti,:] = U[ti-1,:]
	visualizeTraj(ax, {'q':Y[:, 0:3], 'u':U}, model, col='b')
	visualizeTraj(ax, {'q':Ynl[:, 0:3], 'u':U}, model, col='g')
	if rescale:
		ax.set_xlim([-0.05, 0.05])
		ax.set_ylim([-0.05, 0.05])
	print(Y)

	# custom legend
	from matplotlib.lines import Line2D
	custom_lines = [Line2D([0], [0], color='r', alpha=0.3), 
		Line2D([0], [0], color='b', alpha=0.3), 
		Line2D([0], [0], color='g', alpha=0.3)]
	ax.legend(custom_lines, ['LQR', 'Lin OL', 'Nonlin OL'])

	plt.show()