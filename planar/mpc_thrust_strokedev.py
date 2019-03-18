import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
import controlutils.py.MPCUtils as MPCUtils
import FlappingModels

np.set_printoptions(precision=2, suppress=True, linewidth=100)

model = FlappingModels.PlanarThrustStrokeDev()

results = []  # List of sim results: is populated by runMPCSim  below

EXP_SINE = 0
EXP_SOMERSAULT = 1
EXP_11 = 2
exp = EXP_SINE

# Trajectory following?
def getXr(t):
	if exp == EXP_11:
		xr = np.array([t, t, 0.,0.,0.,0.,model.g])
	elif exp == EXP_SINE:
		# Sinusoidal
		xr = np.array([0.5 * np.sin(10 * t), 0.1 * t, 0.,0.,0.,0., model.g])
	elif exp == EXP_SOMERSAULT:
		xr = np.array([0, 0, 3*t,0.,0.,3., model.g])
	else:
		raise 'experiment not implemented'
	# xr = np.array([0.5 * t,0.0, 0,0.,0.,0.])
	# if t < 0.1:
	# 	xr = np.array([0.0, 0.1, 0,0.,0.,0.])
	# else:
	# 	xr = np.array([0.0, 0.1, 0.5 * np.pi * (t - 0.1),0.,0.,0.])

	return xr

def runMPCSim(wx, wu, N=20, dt=0.002, epsi=1e-6, label=''):
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

	model.dt = dt
	mpc = MPCUtils.MPCHelper(model, N, wx, wu, verbose=False, eps_abs=epsi, eps_rel=epsi)

	# Initial and reference states
	# x0 = 0.01 * np.random.rand(model.nx)
	x0, ctrl = model.y0, model.u0

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
		Ad, Bd = model.getLinearDynamics(x0, ctrl)
		x0 = Ad.dot(x0) + Bd.dot(ctrl)
		print(i, x0, ctrl)
		# print(x0horizon)
		X[i, :] = x0
		U[i, :] = ctrl

	results.append({'t':t, 'desTraj':desTraj, 'X':X, 'U':U, 'label':label})

if exp == EXP_SINE:
	# Sine experiments
	runMPCSim([0.01, 0.01, 1, 0, 0, 0, 0], [1, 10000], N=20, dt=0.004, epsi=1e-6, label='N=20')
	runMPCSim([0.01, 0.01, 1, 0, 0, 0, 0], [1, 10000], N=15, dt=0.005, epsi=1e-2, label='N=15')
	runMPCSim([0.01, 0.01, 1, 0, 0, 0, 0], [1, 10000], N=10, dt=0.005, epsi=1e-2, label='N=10')
elif exp == EXP_SOMERSAULT:
	# somersault experiments
	runMPCSim([0.01, 0.01, 10, 0, 0, 1, 0], [1, 10000], N=20, dt=0.004, epsi=1e-2, label='N=20')
	runMPCSim([0.01, 0.01, 10, 0, 0, 1, 0], [1, 10000], N=15, dt=0.004, epsi=1e-2, label='N=15')
	runMPCSim([0.01, 0.01, 10, 0, 0, 1, 0], [1, 10000], N=10, dt=0.004, epsi=1e-2, label='N=10')
else:
	raise 'experiment not implemented'
# print(x0.shape)

nplots = 3 if len(results) > 1 else 2

plt.subplot(nplots,1,1)
for res in results:
	plt.plot(res['t'], res['X'][:, 2],'.-', label=res['label'])
plt.plot(results[0]['t'], results[0]['desTraj'][:,2], 'k--', label='des')
plt.ylabel('phi')

plt.subplot(nplots,1,2)
for res in results:
	plt.plot(res['X'][:, 0], res['X'][:, 1],'.-', label=res['label'])
plt.plot(results[0]['desTraj'][:,0], results[0]['desTraj'][:,1], 'k--', label='des')
plt.legend()
# plt.plot(t, X[:,3])
# plt.plot(t, X[:,4])
plt.ylabel('xz')

# The following are only draw for a single trial
if nplots >= 3:
	plt.subplot(nplots,1,3)
	plt.plot(results[0]['t'], results[0]['U'] * 1000, '.-')
	plt.ylabel('u')

plt.show()
