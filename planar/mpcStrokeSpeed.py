import numpy as np
import FlappingModels
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

model = FlappingModels.PlanarStrokeSpeed()

# crude sim
SIM_DT = 0.001

tvec = np.arange(0, 0.1, SIM_DT)
Nt = len(tvec)

Y = np.zeros((Nt, model.nx))
U = np.zeros((Nt, model.nu))
Y[0,:] = model.y0
U[0,:] = np.array([0,0])

u0 = np.array([1.0,0.0])

def dydt(t, y):
	dydt = model.dydt(y, u0)
	return dydt

def endStroke(t, y):
	return y[-1] - 1e-2

endStroke.terminal = True
endStroke.direction = 1

tf = 0.1
sol = solve_ivp(dydt, [0.0,tf], model.y0, t_eval=np.arange(0,tf,SIM_DT), events=endStroke, dense_output=True)

print(sol.t_events)
print(sol.status)

plt.plot(sol.t, sol.y[0,:])
plt.plot(sol.t, sol.y[0,:])

plt.show()

