import numpy as np
import FlappingModels
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

np.set_printoptions(precision=4, suppress=True, linewidth=200)
model = FlappingModels.PlanarStrokeSpeed()

# crude sim
# SIM_DT = 0.0003

# tvec = np.arange(0, 0.1, SIM_DT)
# Nt = len(tvec)

# Y = np.zeros((Nt, model.nx))
# U = np.zeros((Nt, model.nu))
# Y[0,:] = model.y0
# U[0,:] = np.array([0,0])

def nonLinSim(y0, u0, tf, strokeExtents=[-1.5e-3,1.5e-3], avg=False):
	def dydt(t, y):
		dydt = model.dydt(y, np.asarray(u0), avg=avg)
		return dydt

	strokeEnd = strokeExtents[1]  # will be changed by the code

	def strokeEvent(t, y):
		return y[-1] - strokeEnd

	strokeEvent.terminal = True
	strokeEvent.direction = 1

	t0 = 0.

	sols = []

	tt = np.zeros(0)
	yy = np.zeros((model.nx,0))
	uu = np.zeros((model.nu,0))
	# store only at events times
	tev = np.zeros(0)
	yev = np.zeros((0,model.nx))
	uev = np.zeros((0,model.nu))
	if avg:
		model.sigma0 = y0[-1]
		model.strokeExtents = strokeExtents

	while True:
		# print(y0)
		sol = solve_ivp(dydt, [t0,tf], np.asarray(y0), events=strokeEvent, dense_output=True)
		# sols.append(sol)

		tt = np.hstack((tt, sol.t))
		yy = np.hstack((yy, sol.y))
		uu = np.hstack((uu, ))

		if sol.status == 1:  # termination event
			# restart
			t0 = sol.t_events[0][0]
			y0 = sol.sol(t0)
			if avg:
				model.sigma0 = y0[-1]
			# strokeReset()
			strokeEvent.direction = -strokeEvent.direction
			if strokeEvent.direction > 0:
				strokeEnd = strokeExtents[1]
			else:
				strokeEnd = strokeExtents[0]
			# strokeEnd = -strokeEnd
			u0[0] = -u0[0]
			# logging at event times
			tev = np.hstack((tev, t0))
			yev = np.vstack((yev, y0))
			uev = np.vstack((uev, np.asarray(u0)))
		elif sol.status == 0:
			# end of time
			break
	
	return {'t':tt, 'y':yy, 'u':uu, 'tev':tev, 'yev':yev, 'uev':uev}

def snapshotsPlot(ax, r, col='b'):
	FlappingModels.visualizeTraj(ax, {'t': r['tev'], 'q':r['yev'], 'u':r['uev']}, model, col=col)
			
tf = 0.06
y0 = np.array([0.02,0,0.3,0,0,0,0])
u0 = np.array([1.5,0.0])

# Run some sims
r1 = nonLinSim([0.02,0,0.3,0,0,0,0], u0, tf, strokeExtents=[-2e-3,1.5e-3])
print('Number of strokes =',len(r1['t']))
r2 = nonLinSim([0.0,0,0.3,0,0,0,0], u0, tf, strokeExtents=[-2e-3,1.5e-3], avg=True)

# u0[0] = 2
# r3 = nonLinSim([-0.02,0,0,0,0,0,0], u0, tf, strokeExtents=[-1.5e-3,2e-3])
	
# # Compare to averaged model
# tav = 0
# Yav = np.zeros((model.nx, 1))
# Yav[:,0] = model.y0
# Uav = np.zeros((model.nu, 1))
# Uav[:,0] = np.array([1.5, 1.5e-3])
# for aiter in range(69):
# 	# this is really valid near equilibrium
# 	dydsigma = model.dydt(Yav[:,-1], Uav[:,-1], avg=True)
# 	tmode = 0.001#Uav[1,-1]
# 	ynew = Yav[:,-1] + dydsigma * tmode
# 	Yav = np.hstack((Yav, ynew[:,np.newaxis]))
# 	Uav = np.hstack((Uav, np.array([1.5,1.5e-3])[:,np.newaxis]))

# print(sols[1])

fig, ax = plt.subplots(4)
ax[0].plot(r1['t'], r1['y'][-1,:], label='nl')
ax[0].plot(r2['t'], r2['y'][-1,:], label='avg')
ax[0].set_ylabel('sigma')
ax[0].legend()
ax[0].grid(True)

ax[1].plot(r1['t'], r1['y'][0,:])
ax[1].plot(r1['t'], r1['y'][1,:])
ax[1].set_ylabel('xz')

ax[2].plot(r1['t'], r1['y'][2,:], label='nl')
ax[2].plot(r2['t'], r2['y'][2,:], label='avg')
ax[2].set_ylabel('phi')
ax[0].legend()

ax[3].set_aspect(1)
ax[3].plot(r1['y'][0,:], r1['y'][1,:], '.-', label='nl')
ax[3].plot(r2['y'][0,:], r2['y'][1,:], '.-', label='avg')
# ax[3].plot(Yav[0,:], Yav[1,:], '.-', label='av')
ax[3].grid(True)
ax[3].legend()

ax[-1].set_xlabel('t')


# fig, ax = plt.subplots()
# snapshotsPlot(ax, r1, 'b')
# snapshotsPlot(ax, r2, 'g')
# snapshotsPlot(ax, r3, 'r')
# ax.set_title('Test controllability of nonlinear stroke/speed model')

fig, ax = plt.subplots()
snapshotsPlot(ax, r1, 'b')
snapshotsPlot(ax, r2, 'r')
# print(Yav.shape, Uav.shape)
# FlappingModels.visualizeTraj(ax, {'t': range(69), 'q':Yav.T, 'u':Uav.T}, model, col='r')

plt.show()

