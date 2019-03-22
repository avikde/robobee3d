import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
import controlutils.py.mpc as mpc
import controlutils.py.lqr as lqr
import FlappingModels

np.set_printoptions(precision=2, suppress=True, linewidth=100)

model = FlappingModels.PlanarThrustStrokeDev()

results = []  # List of sim results: is populated by runSim  below

EXP_SINE = 0
EXP_SOMERSAULT = 1
EXP_11 = 2
EXP_GOAL = 3
EXP_VELDES = 4
exp = EXP_SINE

# control types
CTRL_LIN_CUR = 0
CTRL_LIN_HORIZON = 1
CTRL_OPEN_LOOP = 2
CTRL_LQR = 3

# Trajectory following?
def getXr(t):
	if exp in [EXP_GOAL, EXP_VELDES]:
		xr = np.array([-0.03, 0.02, 0, 0, 0, 0])
	elif exp == EXP_11:
		xr = np.array([t, t, 0.,0.,0.,0.])
	elif exp == EXP_SINE:
		# Sinusoidal
		xr = np.array([0.04 * np.sin(10 * t), 0.1 * t, 0.,0.,0.,0.])
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

def runSim(wx, wu, N=20, dt=0.002, epsi=1e-2, label='', ctrlType=CTRL_LQR, nsim=200, x0=None, u0=None):
	'''
	N = horizon (e.g. 20)
	dt = timestep (e.g. 0.002)
	'''
	model.dt = dt
	if ctrlType in [CTRL_LIN_CUR, CTRL_LIN_HORIZON]:
		ltvmpc = mpc.LTVMPC(model, N, wx, wu, verbose=False, scaling=0, eps_abs=epsi, eps_rel=epsi)

	# Initial and reference states
	# x0 = 0.01 * np.random.rand(model.nx)
	if x0 is None:
		x0 = model.y0
	ctrl = model.u0
	if u0 is not None:
		ctrl = u0

	# Simulate in closed loop
	X = np.zeros((nsim, model.nx))
	X[0,:] = x0
	U = np.zeros((nsim, model.nu))
	firstA = 0
		
	t = np.arange(nsim) * dt
	desTraj = np.zeros((nsim,3))

	for i in range(1,nsim):
		# tgoal = (i + N) * dt
		tgoal = i * dt
		if ctrlType == CTRL_LIN_CUR:
			tgoal = (i + N) * dt
		xr = getXr(tgoal)
		# for logging
		desTraj[i,:] = xr[0:3]

		# Try outer loop for the traj
		# velDes = 
		
		if ctrlType == CTRL_LQR:
			# print(x0)
			# get current linearization
			Ad, Bd, fd = model.getLinearDynamics(x0, ctrl)
			# actually compute u
			K = lqr.dlqr(Ad, Bd, np.diag(wx), np.diag(wu))[0]
			ctrl = K @ (xr - x0)
		elif ctrlType == CTRL_LIN_CUR:
			ctrl = ltvmpc.update(x0, ctrl, xr, costMode=mpc.TRAJ)#, trajMode=mpc.ITERATE_TRAJ)
		elif ctrlType == CTRL_LIN_HORIZON:
			# traj to linearize around
			# x0horizon = np.zeros((N,model.nx))
			# for xi in range(model.nx):
			# 	x0horizon[:,xi] = np.linspace(x0[xi], xr[xi], N, endpoint=False)
			# ctrl = mpc.update(x0horizon, np.tile(np.zeros(2),(N,1)), xr)
			# NOTE: left this here, but moved a lot of the traj stuff into MPCUtils through trajMode
			pass
		else: # openloop
			pass
			# ctrl = np.array([1e-3,1e-3])

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

# Run simulations
fig, ax = plt.subplots(ncols=3)

wx = [1000, 1000, 0.05, 5, 5, 0.005]
wu = [0.01,0.01]
nsimi = 200
if exp == EXP_SOMERSAULT:
	wx = [100, 100, 1, 1, 1, 1]
	wu = [0.01,0.01]
	nsimi = 5
if exp == EXP_VELDES:
	wx = [1,1,1, 1, 1, 1]
	wu = [0.01,0.01]

y0 = np.zeros(6)
if exp == EXP_VELDES:
	y0[3:6] = np.random.randn(3)
runSim(wx, wu, dt=0.01, ctrlType=CTRL_LQR, label='LQR', nsim=nsimi, x0=y0)
# MPC
runSim(wx, wu, dt=0.01, ctrlType=CTRL_LIN_CUR, label='MPC', nsim=nsimi, N=10, x0=y0)
if exp == EXP_SOMERSAULT:
	# Openloop
	y0[0] = 0.0
	runSim([], [], dt=0.01, ctrlType=CTRL_OPEN_LOOP, label='OL', nsim=10, x0=y0, u0=np.array([1e-3,1e-3]))

# some colors to draw with
results[0]['col'] = 'r'
results[1]['col'] = 'g'
# results[2]['col'] = 'b'

# Vis
for res in results:
	FlappingModels.visualizeTraj(ax[0], {'q':res['X'][:, 0:3], 'u':res['U']}, model, col=res['col'])
if exp == EXP_SINE:
	ax[0].plot(results[1]['desTraj'][:,0], results[1]['desTraj'][:,1], 'k--', label='des')
elif exp == EXP_GOAL:
	lqrgoal = getXr(0)
	ax[0].plot(lqrgoal[0], lqrgoal[1], 'c*')

# custom legend
from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color=res['col'], alpha=0.3) for res in results]
ax[0].legend(custom_lines, ['LQR', 'MPC', 'OL'])

# Plot time traces
# ax[1].plot(results[1]['X'][:, 0])
# ax[1].plot(results[1]['X'][:, 1])
# ax[1].plot(results[1]['X'][:, 2])
for res in results:
	ax[1].plot(res['X'][:, 3], color=res['col'])
	ax[1].plot(res['X'][:, 4], '--', color=res['col'])
ax[1].set_xlabel('Iters')
ax[1].set_ylabel('dxdz')

for res in results:
	ax[2].plot(res['X'][:, 5], color=res['col'])
ax[2].set_xlabel('Iters')
ax[2].set_ylabel('dphi')

plt.show()
