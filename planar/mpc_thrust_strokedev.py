import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
sys.path.append('..')
import controlutils.py.mpc as mpc
import controlutils.py.lqr as lqr
import FlappingModels
import matplotlib.animation as animation

np.set_printoptions(precision=4, suppress=True, linewidth=200)

# Select here
paramstr = 'hover'
y0 = np.zeros(6)

# Experiment params

paramsets = {
	'defaults': {'eps': 1e-2, 'wx': [1000, 1000, 0.05, 5, 5, 0.005], 'wu': [0.01,0.01], 'nsim': 200, 'dt': 0.01, 'N': 5},
	# Both of these seem to need two "stages" think about what that means
	# SOMERSAULT flip and hover
	'somersault': {'period': 0.2, 'num': 1, 'wx': [[1,1,1, 1, 1, 5], [100,100,1, 10, 10, 0.01]], 'trajMode': mpc.GIVEN_POINT_OR_TRAJ, 'eps': 1e-4, 'wu': [0.01,0.01], 'dt': 0.01},
	# SIDEPERCH position, orientation
	'sideperch': {'period': 0.2, 'periodO': 0.05, 'wx': [[100,100,1, 100, 100, 0.01], [10,10,1, 10, 10, 0.015]], 'wu':[0.01, 0.01], 'dt':0.003, 'N': 15},
	# hover
	'hover': {'wx': [1,1,1, 10, 10, 0.1], 'wu':[0.01, 0.01], 'dt':0.03, 'nsim': 50, 'eps':1e-5},
	# others
	'waypoint': {},
	'sine': {}
}
params = paramsets[paramstr]

if paramstr == 'somersault':
	params['nsim'] = int((params['period'] * params['num'] + 0.01)/ params['dt'])
elif paramstr == 'sideperch':
	params['nsim'] = int((params['period'] + params['periodO'])/params['dt'])
elif paramstr == 'hover':
	y0[3:6] = np.array([1,0,10])  # give it some initial velocity

# get parameters or default
def _pget(pname):
	# params from dict (DEFAULT PARAMS HERE)
	return params.get(pname, paramsets['defaults'][pname])

# control types
CTRL_MPC = 0
CTRL_LIN_HORIZON = 1
CTRL_OPEN_LOOP = 2
CTRL_LQR = 3

saveMovie = 0  #1 shows movie, 2 saves, 3 plots arrow
model = FlappingModels.PlanarThrustStrokeDev(rescale=True)
results = []  # List of sim results: is populated by runSim  below

# Trajectory following?
def goal(t):
	if paramstr == 'somersault':
		omega = 2 * np.pi / params['period']
		if t < params['num'] * params['period']:
			xr = np.array([0, 0, omega*t,0.,0.,omega])
		else:
			xr = np.array([0, 0, 2 * np.pi * params['num'],0.,0.,0])
	elif paramstr == 'sideperch':
		angDes = np.pi
		omega = angDes / params['period']
		if t < params['period']:
			xr = np.array([0.2 * t, 0, 0, 0.2,0.,0])
		elif t < params['periodO']:
			xr = np.array([0.2 * params['period'], 0, omega * (t - params['period']),0.,0.,omega])
		else:
			xr = np.array([0.2 * params['period'], 0, angDes,0.,0.,0])
	elif paramstr in ['hover', 'waypoint']:
		xr = np.array([-0.03, 0.02, 0, 0, 0, 0])
	# elif exp == EXP_11:
	# 	xr = np.array([t, t, 0.,0.,0.,0.])
	elif paramstr == 'sine':
		# Sinusoidal
		xr = np.array([0.04 * np.sin(10 * t), 0.1 * t, 0.,0.,0.,0.])
	else:
		raise 'experiment not implemented'
	# xr = np.array([0.5 * t,0.0, 0,0.,0.,0.])
	# if t < 0.1:
	# 	xr = np.array([0.0, 0.1, 0,0.,0.,0.])
	# else:
	# 	xr = np.array([0.0, 0.1, 0.5 * np.pi * (t - 0.1),0.,0.,0.])

	return xr

def runSim(params, label='', ctrlType=CTRL_LQR, x0=None, u0=None):
	'''
	N = horizon (e.g. 20)
	dt = MPC timestep (e.g. 0.002), different from sim timestep
	'''
	# params from dict (DEFAULT PARAMS HERE)
	peps = _pget('eps')
	pwx = _pget('wx')
	pwu = _pget('wu')
	nsim = _pget('nsim')
	dt = _pget('dt')
	N = _pget('N')

	# Check if a list of weights has been provided (multi-part behavior)
	weightsIdx = 0
	wx = pwx[weightsIdx] if isinstance(pwx[0], list) else pwx
	# wu = pwu[0] if isinstance(pwu, list) else pwu
	wu = pwu

	SIM_DT = 0.0005  # not the same as "dt" which is MPC dt
	model.dt = dt
	if ctrlType in [CTRL_MPC, CTRL_LIN_HORIZON]:
		# TODO: confirm this weight scaling
		wx = np.array(wx) / dt
		# wu = np.array(wu) / dt

		# default for eps_prim_inf = 1.0e-04, eps_dual_inf = 1.0e-04
		# print(peps)
		# sys.exit(0)
		ltvmpc = mpc.LTVMPC(model, N, wx, wu, verbose=False, polish=False, scaling=0, eps_rel=peps, eps_abs=peps, max_iter=1000000, kdamping=0)
		# ltvmpc.MAX_ULIM_VIOL_FRAC = 0.5
		# sys.exit(0)

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
		
	tvec = np.arange(nsim) * dt
	desTraj = np.zeros((nsim,3))

	for i in range(1,nsim):
		# tgoal = (i + N) * dt
		tgoal = tvec[i]
		if ctrlType == CTRL_MPC:
			tgoal = (i + N) * dt
		xr = goal(tgoal)
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
		elif ctrlType == CTRL_MPC:
			if ((paramstr == 'somersault' and tgoal > params['period'] * params['num'] * 1.5) or (paramstr == 'sideperch' and tgoal > params['period'])) and weightsIdx == 0:
				weightsIdx += 1
				print('Updating weights to', pwx[weightsIdx])
				ltvmpc.updateWeights(wx=np.array(pwx[weightsIdx])/dt)
			# elif exp == EXP_SIDEPERCH and tgoal > sidePerchParams['period'] + sidePerchParams['periodO']:
			# 	ltvmpc.updateWeights(wx=np.array(sidePerchParams['wx'][0])/dt)
			ctrl = ltvmpc.update(x0, xr, costMode=mpc.FINAL)
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

		# simulate forward -- faster than MPC rate
		for simiter in range(int(dt/SIM_DT)):
			ydot = model.dydt(x0, ctrl)
			x0 += ydot * SIM_DT  # euler integration
		# x0 = model.dynamics(x0, ctrl)
		print(i, x0, ctrl)
		# print(x0horizon)
		X[i, :] = x0
		U[i, :] = ctrl

	results.append({'t':tvec, 'desTraj':desTraj, 'X':X, 'U':U, 'label':label})

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
# runSim(params, ctrlType=CTRL_LQR, label='LQR', nsim=1, x0=y0)
# results[0]['col'] = 'r'
# MPC
runSim(params, ctrlType=CTRL_MPC, label='MPC', x0=y0)
results[0]['col'] = 'g'
# if exp == EXP_SOMERSAULT:
# 	# Openloop
# 	y0[0] = 0.0
# 	runSim([], [], dt=0.01, ctrlType=CTRL_OPEN_LOOP, label='OL', nsim=10, x0=y0, u0=np.array([1e-3,1e-3]))
# 	results[2]['col'] = 'b'

# Vis
if saveMovie > 0:
	
	from matplotlib.collections import PatchCollection
	from matplotlib.patches import Polygon, Arrow

	# Set up formatting for the movie files
	Writer = animation.writers['ffmpeg']
	writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

	# traj
	res = results[0]  # mpc
	traj = {'q':res['X'][:, 0:3], 'u':res['U']}
	N = traj['q'].shape[0]

	robotBodies = []

	fig, ax = plt.subplots()
	# ax.autoscale_view(True)
	seline, = ax.plot([], [], 'k--', linewidth=1, alpha=0.5)
	# aeroArrow = Arrow(0,0,0,0.01, width=0.0002, alpha=0.3, facecolor=res['col'])
	bodyPatch = Polygon([[0,0],[0,0],[0,0],[0,0]], alpha=0.5, facecolor=res['col'], edgecolor='k')
	if saveMovie == 3:
		Q = ax.quiver([0], [0], [0], [0.1], pivot='mid', color=res['col'], units='dots')

	def init():
		ax.set_aspect(1)
		ax.set_xlim([-0.05,0.05])
		ax.set_ylim([-0.05,0.05])
		ax.grid(True)
		ax.set_xlabel('x')
		ax.set_ylabel('z')
		# sets all the patches in the collection to white
		# pc.set_color(colors)
		ax.add_patch(bodyPatch)
		# arPatch = ax.add_patch(aeroArrow)
		# in the order body patch, stroke extents line
		return bodyPatch,

	def animate(k):
		qk = traj['q'][k,:]
		uk = traj['u'][k,:]

		# get info from model
		body, pcop, Faero, strokeExtents = model.visualizationInfo(qk, uk, rawxy=True)

		bodyPatch.set_xy(body)
		
		# 	ax.plot(, 'k--', linewidth=1,  alpha=0.3)
		# print(dir(aeroArrow))
		# sys.exit(-1)

		seline.set_data(strokeExtents[:,0], strokeExtents[:,1])
		# print(pcop)
		if saveMovie == 3:
			Q.set_offsets(np.array([pcop[0], pcop[1]]))
			Q.set_UVC(Faero[0], Faero[1])

		return bodyPatch, 

	ani = animation.FuncAnimation(fig, init_func=init, func=animate, frames=N, interval=_pget('dt'), blit=True)
	if saveMovie in [2,3]:
		ani.save('test.mp4', writer=writer)
	elif saveMovie == 1:
		plt.show()

else:
	# Regular plots
	fig, ax = plt.subplots(nrows=5)

	for res in results:
		FlappingModels.visualizeTraj(ax[0], {'q':res['X'][:, 0:3], 'u':res['U']}, model, col=res['col'])
	if paramstr == 'sine':
		ax[0].plot(results[1]['desTraj'][:,0], results[1]['desTraj'][:,1], 'k--', label='des')
	elif paramstr == 'waypoint':
		lqrgoal = goal(0)
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
		ax[1].plot(res['t'], res['X'][:, 3], color=res['col'])
		ax[1].plot(res['t'], res['X'][:, 4], '--', color=res['col'])
	ax[1].set_ylabel('dxdz')

	for res in results:
		ax[2].plot(res['t'], res['X'][:, 5], color=res['col'])
	ax[2].set_ylabel('dphi')
	
	# u stuff
	
	for res in results:
		ax[3].plot(res['t'], res['U'][:, 0], color=res['col'])
	umin, umax, _, _ = model.getLimits()
	ax[3].axhline(umin[0])
	ax[3].axhline(umax[0])
	ax[3].set_ylabel('u0')

	for res in results:
		ax[4].plot(res['t'], res['U'][:, 1], color=res['col'])
	umin, umax, _, _ = model.getLimits()
	ax[4].axhline(umin[1])
	ax[4].axhline(umax[1])
	ax[4].set_ylabel('u1')

	ax[-1].set_xlabel('t (sec)')

	plt.show()
