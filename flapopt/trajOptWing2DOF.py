import autograd.numpy as np
from autograd import jacobian
import matplotlib.pyplot as plt
from matplotlib import animation
import sys
sys.path.append('..')
import planar.FlappingModels as FlappingModels
from scipy.integrate import solve_ivp
from matplotlib import animation
from matplotlib.collections import PatchCollection
np.set_printoptions(precision=4, suppress=True, linewidth=200)
from mosek.fusion import *
import mosek as msk

# MOSEK test ------------------------------------------------------

def primal_problem(P):

    print(msk.Env.getversion())
    
    k= len(P)
    if k==0: return -1,[]

    n= len(P[0])
    
    with Model("minimal sphere enclosing a set of points - primal") as M:

        r0 = M.variable(1    , Domain.greaterThan(0.))
        p0 = M.variable([1,n], Domain.unbounded())

        R0 = Var.repeat(r0,k)
        P0 = Var.repeat(p0,k)
       
        M.constraint( Expr.hstack( R0, Expr.sub(P0 , P) ), Domain.inQCone())

        M.objective(ObjectiveSense.Minimize, r0)
        M.setLogHandler(open('logp','wt'))

        M.solve()

        return r0.level()[0], p0.level()
# test problem
import random
n = 2
k = 500
p=  [ [random.gauss(0.,10.) for nn in range(n)] for kk in range(k)]
r0,p0 = primal_problem(p)

print ("r0^* = ", r0)
print ("p0^* = ", p0)

# ---------------------------------------------------

m = FlappingModels.Wing2DOF()

# discrete => do not need solve_ivp

dt = 1e-4
tf = 0.1
tvec = np.arange(0, tf, dt)
yi = np.zeros((m.nx, len(tvec)))
yi[:,0] = np.array([1e-2, 0, 0, 0])

# params
params = []

# Functions ---

def sigmades(t):
    return 15e-3 * np.sin(150 * 2 * np.pi * t)
def controller(t, y):
    return 1e2 * (sigmades(t) - y[0]) - 1e-2 * y[2]
def closedLoop(t, y):
    # u0 = 1e-3 * np.sin(100 * 2 * np.pi * tvec[ti])
    # pos servoing
    u0 = controller(t, y)
    # print(u0)
    return m.dydt(y, [u0], params)

def Jobjinst(y, u, params):
    _, Faero = m.aero(y, u, params)
    return -Faero[1] # minimization

def Jcostinst_dynpenalty(ynext, y, u, params):
    '''error on dynamics for penalty method'''
    dynErr = ynext - (y + m.dydt(y, u, params) * dt)
    return 1/2 * dynErr.T @ dynErr

# FIXME: need to find y that make one cycle
# as a first pass just average over the whole time

def Jcosttraj(yu, params):
    '''this is over a traj. yu = (nx+nu,Nt)-shaped'''
    Nt = yu.shape[1]
    c = 0
    PENALTY = 1e-6
    for i in range(Nt-1):
        c += Jobjinst(yu[:m.nx,i], yu[m.nx:,i], params) + PENALTY * Jcostinst_dynpenalty(yu[:m.nx,i+1], yu[:m.nx,i], yu[m.nx:,i], params)
    # TODO: any final cost?
    c += Jobjinst(yu[:m.nx,-1], yu[m.nx:,-1], params)
    return c

def Jcostsol(solt, soly, params):
    Nt = len(solt)
    Jcosti = 0
    for i in range(Nt):
        Jcosti += Jobjinst(soly[:,i], controller(solt[i], soly[:,i]), params)
    return Jcosti

Jgrad = jacobian(lambda yu : Jcosttraj(yu, params))

# ---

for ti in range(1, len(tvec)):
    yi[:,ti] = yi[:,ti-1] + dt * closedLoop(tvec[ti], yi[:,ti-1])


# compare to continuous
sol = solve_ivp(closedLoop, [0,tf], yi[:,0], dense_output=True, t_eval=tvec)

print('Avg cost =', Jcostsol(sol.t, sol.y, params))

# Trajectory to begin gradient descent from ------------
yu0 = sol.y.copy()
# yu0 = yi.copy()
utest = np.zeros(yu0.shape[1])
for i in range(yu0.shape[1]):
    utest[i] = controller(sol.t[i], sol.y[:,i])
yu0 = np.vstack((yu0, utest))

# Test compute gradient
print(Jcosttraj(yu0, params))
g1 = Jgrad(yu0)
# print(g1[:20])
# # Stupid "line search" for step size
# for i in range(-10,10):
#     print(i, Jcosttraj(yu0 - np.power(10.0,i) * g1, params))
# sys.exit(0)
yu1 = yu0 - 1e1 * g1
print(Jcosttraj(yu1, params))

# --------------------------------------------------------

# plots

fig, ax = plt.subplots(3)

# display


# ax[0].plot(tvec, yi[0,:])
ax[0].plot(tvec, yu0[0,:], label='0')
ax[0].plot(tvec, yu1[0,:], label='1')
ax[0].plot(tvec, sigmades(tvec), 'k--')
ax[0].legend()

# ax[1].plot(tvec, yi[1,:])
ax[1].plot(tvec, yu0[1,:], label='0')
ax[1].plot(tvec, yu1[1,:], label='1')
ax[1].legend()

# def makeAnim(_ax, _t, _y):

def flapkin(yui, xyoff):
    # wing extents
    wing1 = np.array([yui[0], 0]) + np.asarray(xyoff)
    c, s = np.cos(yui[1]), np.sin(yui[1])
    wing2 = wing1 + np.array([[c, -s], [s, c]]) @ np.array([0, -2*m.cbar])
    # aero arrow extents
    _, Faero = m.aero(yui[:m.nx], yui[m.nx:])
    pcop = (wing1 + wing2)/2
    aeroEnd = pcop + 0.3 * Faero
    return wing1, wing2, pcop, aeroEnd

_ax = ax[2]

p1, = _ax.plot([], [], 'b.-', linewidth=4)
paero, = _ax.plot([], [], 'r', linewidth=2)
p2, = _ax.plot([], [], 'b.-', linewidth=4)
paero2, = _ax.plot([], [], 'r', linewidth=2)
_ax.grid(True)
_ax.set_aspect(1)
_ax.set_ylim(m.cbar * np.array([-4, 2]))

def _init():
    return p1, paero, p2, paero2, 

def _animate(i):
    wing1, wing2, pcop, aeroEnd = flapkin(yu0[:,i], [-0.01,0])
    p1.set_xdata([wing1[0], wing2[0]])
    p1.set_ydata([wing1[1], wing2[1]])
    paero.set_xdata([pcop[0], aeroEnd[0]])
    paero.set_ydata([pcop[1], aeroEnd[1]])
    wing1, wing2, pcop, aeroEnd = flapkin(yu1[:,i], [0.01,0])
    p2.set_xdata([wing1[0], wing2[0]])
    p2.set_ydata([wing1[1], wing2[1]])
    paero2.set_xdata([pcop[0], aeroEnd[0]])
    paero2.set_ydata([pcop[1], aeroEnd[1]])
    return p1, paero, p2, paero2, 

anim = animation.FuncAnimation(fig, _animate, init_func=_init, frames=yu0.shape[1], interval=1e5*dt, blit=True)
# makeAnim(ax[2], sol.t, sol.y)

plt.tight_layout()
# print(sol.y.shape)

plt.show()
