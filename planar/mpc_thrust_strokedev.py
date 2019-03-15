import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
import MPCUtils
import FlappingModels

np.set_printoptions(precision=2, suppress=True, linewidth=100)

model = FlappingModels.PlanarThrustStrokeDev()

# Trajectory following?
def getXr(t):
	# xr = np.array([t, t, 0.,0.,0.,0.,model.g])

	# Sinusoidal
	xr = np.array([0.5 * np.sin(10 * t), 0.1 * t, 0.,0.,0.,0., model.g])
	
	# xr = np.array([0, 0, 3 * np.pi *t,0.,0.,0., model.g])
	
	# xr = np.array([0.5 * t,0.0, 0,0.,0.,0.])
	# if t < 0.1:
	# 	xr = np.array([0.0, 0.1, 0,0.,0.,0.])
	# else:
	# 	xr = np.array([0.0, 0.1, 0.5 * np.pi * (t - 0.1),0.,0.,0.])

	return xr

def runMPCSim(wx, wu, N=20, dt=0.002, epsi=1e-6):
	'''
	N = horizon (e.g. 20)
	dt = timestep (e.g. 0.002)

	Learnings about parameters:
	- tuning parameters: weights, N, dt, eps_abs, eps_rel
	- TODO: reason about units for wxi, wui and reduce those to two scalars?
	- longer N => problem harder to solve => infeasible result more likely
	- shorter N and/or longer dt => instability more likely
	- eps_i lower => infeasible more likely (TODO: needs more testing)
	'''
	# Sim parameters
	nsim = 200
	# control types
	CTRL_LIN_CUR = 0
	CTRL_LIN_HORIZON = 1
	CTRL_OPEN_LOOP = 2
	ctrlType = CTRL_LIN_CUR # could make this a param

	mpc = MPCUtils.MPCHelper(model, N, wx, wu, verbose=False, eps_abs=epsi, eps_rel=epsi)

	# Initial and reference states
	# x0 = 0.01 * np.random.rand(model.nx)
	x0, ctrl = model.getInit()

	# Simulate in closed loop
	X = np.zeros((nsim, model.nx))
	U = np.zeros((nsim, model.nu))
	firstA = 0
		
	t = np.arange(nsim) * dt
	desTraj = np.zeros((nsim,3))

	for i in range(nsim):
		# tgoal = (i + N) * dt
		tgoal = (i + 1) * dt
		xr = getXr(tgoal)
		# for logging
		desTraj[i,:] = xr[0:3]

		if ctrlType == CTRL_LIN_CUR:
			ctrl = mpc.update(x0, np.zeros(2), xr, dt)
		elif ctrlType == CTRL_LIN_HORIZON:
			# traj to linearize around
			x0horizon = np.zeros((N,model.nx))
			for xi in range(model.nx):
				x0horizon[:,xi] = np.linspace(x0[xi], xr[xi], N, endpoint=False)
			ctrl = mpc.update(x0horizon, np.tile(np.zeros(2),(N,1)), xr, dt)
		else: # openloop
			ctrl = np.array([1e-6,0])

		# simulate forward
		Ad, Bd = model.getLin(x0, ctrl, dt)
		x0 = Ad.dot(x0) + Bd.dot(ctrl)
		print(i, x0, ctrl)
		# print(x0horizon)
		X[i, :] = x0
		U[i, :] = ctrl

	return t, desTraj, X, U

t1, desTraj1, X1, U1 = runMPCSim([0.01, 0.01, 1, 0, 0, 0, 0], [1, 10000], N=20, dt=0.004)
t2, desTraj2, X2, U2 = runMPCSim([0.01, 0.01, 1, 0, 0, 0, 0], [1, 10000], N=15, dt=0.004)
t3, desTraj3, X3, U3 = runMPCSim([0.01, 0.01, 1, 0, 0, 0, 0], [1, 10000], N=10, dt=0.004)
labels = ['N=20','N=15','N=10']

# print(x0.shape)
plt.subplot(2,1,1)
plt.plot(t1, X1[:, 2],'.-', label=labels[0])
plt.plot(t2, X2[:, 2],'.-', label=labels[1])
plt.plot(t3, X3[:, 2],'.-', label=labels[2])
plt.plot(t1, desTraj1[:,2], 'k--', label='des traj')
plt.ylabel('phi')

plt.subplot(2,1,2)
plt.plot(X1[:,0], X1[:,1],'.-', label=labels[0])
plt.plot(X2[:,0], X2[:,1],'.-', label=labels[1])
plt.plot(X3[:,0], X3[:,1],'.-', label=labels[2])
plt.plot(desTraj1[:,0], desTraj1[:,1], 'k--', label='des traj')
plt.legend()
# plt.plot(t, X[:,3])
# plt.plot(t, X[:,4])

plt.ylabel('xz')

# plt.subplot(3,1,3)
# plt.plot(t1, 1000*U1,'.-')
# plt.ylabel('u')
plt.show()
