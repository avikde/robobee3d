import numpy as np
import FlappingModels
from scipy.integrate import solve_ivp

model = FlappingModels.PlanarStrokeSpeed()

# crude sim
SIM_DT = 0.001

tvec = np.arange(0, 0.1, SIM_DT)
Nt = len(tvec)

Y = np.zeros((Nt, model.nx))
U = np.zeros((Nt, model.nu))
Y[0,:] = model.y0
U[0,:] = np.array([0,0])

# for ti in range(1, len(tvec)):
# 	dydt = model.dydt(Y[ti-1,:], U[ti-1,:])
# 	# first order integration
# 	Y[ti,:] = Y[ti-1,:] + dydt * SIM_DT

# 	# integrate till event?
# 	if model.guard(Y[ti,:])

def upward_cannon(t, y): return [y[1], -0.5]
def hit_ground(t, y): return y[1]
hit_ground.terminal = True
hit_ground.direction = -1
sol = solve_ivp(upward_cannon, [0, 100], [0, 10], events=hit_ground)
print(sol.t_events)

print(sol)

