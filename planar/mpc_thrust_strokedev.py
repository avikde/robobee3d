import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
import MPCUtils
import FlappingModels

np.set_printoptions(precision=2, suppress=True, linewidth=100)

# Sim parameters
nsim = 100
N = 15 #horizon

model = FlappingModels.PlanarThrustStrokeDev()
mpc = MPCUtils.MPCHelper(model, N, [100, 100, 100, 0, 0, 0, 0], [1000, 1000], verbose=False)

# Initial and reference states
# x0 = 0.01 * np.random.rand(model.nx)
x0, ctrl = model.getInit()

# Simulate in closed loop
X = np.zeros((nsim, model.nx))
U = np.zeros((nsim, model.nu))
dt = 0.1
firstA = 0

# Trajectory following?
def getXr(t):
	xr = np.array([t, 0, 0,0.,0.,0.,model.g])
	# xr = np.array([0.5 * np.sin(1 * t), 0.05 * t, 0.,0.,0.,0., model.g])
	
	# xr = np.array([0.5 * t,0.0, 0,0.,0.,0.])
	# if t < 0.1:
	# 	xr = np.array([0.0, 0.1, 0,0.,0.,0.])
	# else:
	# 	xr = np.array([0.0, 0.1, 0.5 * np.pi * (t - 0.1),0.,0.,0.])

	return xr
	
t = np.arange(nsim) * dt
desTraj = np.zeros((nsim,3))

for i in range(nsim):
	tgoal = (i + N) * dt
	xr = getXr(tgoal)
	# for logging
	desTraj[i,:] = xr[0:3]

	# traj to linearize around
	x0horizon = np.zeros((N,model.nx))
	for xi in range(model.nx):
		x0horizon[:,xi] = np.linspace(x0[xi], xr[xi], N, endpoint=False)
	ctrl = mpc.update(x0horizon, np.zeros((N,2)), xr, dt)
	# ctrl = mpc.update(x0, xr, dt)
	# ctrl = np.array([1e-6,0])

	# simulate forward
	Ad, Bd = model.getLin(x0, ctrl, dt)
	x0 = Ad.dot(x0) + Bd.dot(ctrl)
	print(i, x0, ctrl)
	# print(x0horizon)
	X[i, :] = x0
	U[i, :] = ctrl

# print(x0.shape)
plt.subplot(3,1,1)
plt.plot(t, X[:, 2],'.-')
plt.plot(t, desTraj[:,2], 'k--', label='des traj')
plt.ylabel('phi')
plt.subplot(3,1,2)
plt.plot(X[:,0], X[:,1],'.-', label='traj')
plt.plot(desTraj[:,0], desTraj[:,1], 'k--', label='des traj')

# plt.plot(t, X[:,3])
# plt.plot(t, X[:,4])

plt.ylabel('xz')
plt.legend()
plt.subplot(3,1,3)
plt.plot(t, 1000*U,'.-')
plt.ylabel('u')
plt.show()
