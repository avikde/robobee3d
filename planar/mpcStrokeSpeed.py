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


def dydt(t, y):
	dydt = model.dydt(y, u0)
	return dydt

strokeEnd = 1.5e-3  # will be changed by the code

def strokeEvent(t, y):
	return y[-1] - strokeEnd

strokeEvent.terminal = True
strokeEvent.direction = 1

tf = 0.05
t0 = 0.
y0 = np.array([0,0,0.3,0,0,0,0])
u0 = np.array([1.5,0.0])

sols = []

tt = np.zeros(0)
yy = np.zeros((model.nx,0))
uu = np.zeros((model.nu,0))
# store only at events times
tev = np.zeros(0)
yev = np.zeros((0,model.nx))
uev = np.zeros((0,model.nu))

while True:
	print(y0)
	sol = solve_ivp(dydt, [t0,tf], y0, events=strokeEvent, dense_output=True)
	# sols.append(sol)

	tt = np.hstack((tt, sol.t))
	yy = np.hstack((yy, sol.y))
	uu = np.hstack((uu, ))

	if sol.status == 1:  # termination event
		# restart
		t0 = sol.t_events[0][0]
		y0 = sol.sol(t0)
		# strokeReset()
		strokeEnd = -strokeEnd
		strokeEvent.direction = -strokeEvent.direction
		u0[0] = -u0[0]
		# logging at event times
		tev = np.hstack((tev, t0))
		yev = np.vstack((yev, y0))
		uev = np.vstack((uev, u0))
	elif sol.status == 0:
		# end of time
		break
	
# Compare to averaged model
tav = 0
Yav = np.zeros((model.nx, 1))
Yav[:,0] = model.y0
Uav = np.zeros((model.nu, 1))
Uav[:,0] = np.array([1.0, 0.0])
for aiter in range(2):
	# this is really valid near equilibrium
	dydsigma = model.dydt(Yav[:,-1], Uav[:,-1], avg=True)
	tmode = 0.001#Uav[1,-1]
	ynew = Yav[:,-1] + dydsigma * tmode
	Yav = np.hstack((Yav, ynew[:,np.newaxis]))

# print(sols[1])

fig, ax = plt.subplots(4)
ax[0].plot(tt, yy[-1,:])
ax[0].set_ylabel('sigma')

ax[1].plot(tt, yy[0,:])
ax[1].plot(tt, yy[1,:])
ax[1].set_ylabel('xz')

ax[2].plot(tt, yy[2,:])
ax[2].set_ylabel('phi')

ax[3].set_aspect(1)
ax[3].plot(yy[0,:], yy[1,:], '.-', label='nl')
ax[3].plot(Yav[0,:], Yav[1,:], '.-', label='av')
ax[3].grid(True)
ax[3].legend()

ax[-1].set_xlabel('t')


fig, ax = plt.subplots()
FlappingModels.visualizeTraj(ax, {'t': tev, 'q':yev, 'u':uev}, model, col='b', xylim=[-0.05,0.05,-0.05,0.05])

plt.show()

