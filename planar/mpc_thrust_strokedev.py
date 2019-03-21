import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
import controlutils.py.MPCUtils as MPCUtils
import controlutils.py.lqr as lqr
import FlappingModels

np.set_printoptions(precision=2, suppress=True, linewidth=100)

model = FlappingModels.PlanarThrustStrokeDev()

results = []  # List of sim results: is populated by runSim  below

EXP_SINE = 0
EXP_SOMERSAULT = 1
EXP_11 = 2
EXP_GOAL = 3
exp = EXP_GOAL

# control types
CTRL_LIN_CUR = 0
CTRL_LIN_HORIZON = 1
CTRL_OPEN_LOOP = 2
CTRL_LQR = 3

# Trajectory following?
def getXr(t):
	if exp == EXP_GOAL:
		xr = np.array([-0.03, -0.02, 0, 0, 0, 0])
	elif exp == EXP_11:
		xr = np.array([t, t, 0.,0.,0.,0.])
	elif exp == EXP_SINE:
		# Sinusoidal
		xr = np.array([0.5 * np.sin(10 * t), 0.1 * t, 0.,0.,0.,0.])
	elif exp == EXP_SOMERSAULT:
		xr = np.array([0, 0, 3*t,0.,0.,3.])
	else:
		raise 'experiment not implemented'
	# xr = np.array([0.5 * t,0.0, 0,0.,0.,0.])
	# if t < 0.1:
	# 	xr = np.array([0.0, 0.1, 0,0.,0.,0.])
	# else:
	# 	xr = np.array([0.0, 0.1, 0.5 * np.pi * (t - 0.1),0.,0.,0.])

	return xr

def runSim(wx, wu, N=20, dt=0.002, epsi=1e-2, label='', ctrlType=CTRL_LQR, nsim=200):
	'''
	N = horizon (e.g. 20)
	dt = timestep (e.g. 0.002)
	'''
	model.dt = dt
	if ctrlType in [CTRL_LIN_CUR, CTRL_LIN_HORIZON]:
		mpc = MPCUtils.LTVMPC(model, N, wx, wu, verbose=False, scaling=0, eps_abs=epsi, eps_rel=epsi)

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
		
		if ctrlType == CTRL_LQR:
			# print(x0)
			# get current linearization
			Ad, Bd, fd = model.getLinearDynamics(x0, ctrl)
			# actually compute u
			K = lqr.dlqr(Ad, Bd, np.diag(wx), np.diag(wu))[0]
			ctrl = K @ (xr - x0)
		elif ctrlType == CTRL_LIN_CUR:
			ctrl = mpc.update(x0, np.zeros(2), xr)
			# ctrl = mpc.update(x0, np.zeros(2), xr, trajMode=mpc.ITERATE_TRAJ_LIN)
		elif ctrlType == CTRL_LIN_HORIZON:
			# traj to linearize around
			# x0horizon = np.zeros((N,model.nx))
			# for xi in range(model.nx):
			# 	x0horizon[:,xi] = np.linspace(x0[xi], xr[xi], N, endpoint=False)
			# ctrl = mpc.update(x0horizon, np.tile(np.zeros(2),(N,1)), xr)
			# NOTE: left this here, but moved a lot of the traj stuff into MPCUtils through trajMode
			pass
		else: # openloop
			ctrl = np.array([1e-6,0])

		# simulate forward
		x0 = model.dynamics(x0, ctrl)
		print(i, x0, ctrl)
		# print(x0horizon)
		X[i, :] = x0
		U[i, :] = ctrl

	results.append({'t':t, 'desTraj':desTraj, 'X':X, 'U':U, 'label':label})

# if exp == EXP_SINE:
# 	# Sine experiments
# 	runSim([0.01, 0.01, 1, 0, 0, 0, 0], [1, 10000], N=20, dt=0.004, label='N=20')
# 	runSim([0.01, 0.01, 1, 0, 0, 0, 0], [1, 10000], N=15, dt=0.005, label='N=15')
# 	runSim([0.01, 0.01, 1, 0, 0, 0, 0], [1, 10000], N=10, dt=0.005, label='N=10')
# 	# # runSim([0.01, 0.01, 1, 0, 0, 0, 0], [1, 10000], N=2, dt=0.005, epsi=1e-2, label='N=20')
# 	# # runSim([0.01, 0.01, 1, 0, 0, 0, 0], [1, 10000], N=3, dt=0.005, epsi=1e-2, label='N=15')
# 	# runSim([0.01, 0.01, 1, 0, 0, 0, 0], [1, 10000], N=1, dt=0.005, epsi=1e-2, label='N=10')
# elif exp == EXP_SOMERSAULT:
# 	# somersault experiments
# 	runSim([0.01, 0.01, 10, 0, 0, 1, 0], [1, 10000], N=10, dt=0.004, label='N=10')
# 	runSim([0.01, 0.01, 10, 0, 0, 1, 0], [1, 10000], N=15, dt=0.004, label='N=15')
# else:
# 	raise 'experiment not implemented'
# # print(x0.shape)

# nplots = 3 if len(results) > 1 else 2

# plt.subplot(nplots,1,1)
# for res in results:
# 	plt.plot(res['t'], res['X'][:, 2],'.-', label=res['label'])
# plt.plot(results[0]['t'], results[0]['desTraj'][:,2], 'k--', label='des')
# plt.ylabel('phi')

# plt.subplot(nplots,1,2)
# for res in results:
# 	plt.plot(res['X'][:, 0], res['X'][:, 1],'.-', label=res['label'])
# plt.plot(results[0]['desTraj'][:,0], results[0]['desTraj'][:,1], 'k--', label='des')
# plt.legend()
# # plt.plot(t, X[:,3])
# # plt.plot(t, X[:,4])
# plt.ylabel('xz')

# # The following are only draw for a single trial
# if nplots >= 3:
# 	plt.subplot(nplots,1,3)
# 	plt.plot(results[0]['t'], results[0]['U'] * 1000, '.-')
# 	plt.ylabel('u')

# plt.show()

np.set_printoptions(suppress=True, linewidth=100, precision=3)

model = FlappingModels.PlanarThrustStrokeDev()
model.dt = 0.01

# For visualization
fig, ax = plt.subplots(1)
Ndraw = 100
# Use runSim
runSim([100,100,10,1,1,1], [0.1,0.1], dt=0.01, ctrlType=CTRL_LQR, label='LQR')
FlappingModels.visualizeTraj(ax, {'q':results[0]['X'][:, 0:3], 'u':results[0]['U']}, model)
lqrgoal = getXr(0)
ax.plot(lqrgoal[0], lqrgoal[1], 'c*')

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
FlappingModels.visualizeTraj(ax, {'q':Y[:, 0:3], 'u':U}, model, col='b')
print(Y)

# custom legend
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='r', alpha=0.3), 
	Line2D([0], [0], color='b', alpha=0.3)]
ax.legend(custom_lines, ['LQR', 'Linearized'])

plt.show()
