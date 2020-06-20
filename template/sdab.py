import time, subprocess
import numpy as np
import pybullet as p
import robobee
import viewlog
np.set_printoptions(precision=2, suppress=True, linewidth=200)

# p.DIRECT for non-graphical
bee = robobee.RobobeeSim(p.GUI, slowDown=1, camLock=True, timestep=0.1, gui=1)
# load robot
startPos = [0,0,10]
startOrientation = p.getQuaternionFromEuler(np.zeros(3))
subprocess.call(["python", "../urdf/xacro.py", "../urdf/sdab.xacro", "-o", "../urdf/sdab.urdf"])
bid = bee.load("../urdf/sdab.urdf", startPos, startOrientation, useFixedBase=False)
data = viewlog.initLog()

# # Helper function: traj to track
# def traj(t):
#     return startPos + np.array([30 * np.sin(0.002*np.pi*t), 0, 0.5 * t])

# # draw traj
# T_END=1000
# tdraw = np.linspace(0, T_END, 20)
# for ti in range(1, len(tdraw)):
#     p.addUserDebugLine(traj(tdraw[ti-1]), traj(tdraw[ti]), lineColorRGB=[0,0,0], lifeTime=0)
    
# ---

# Params stored as (min, max, default) tuples
params = {'freq': (0, 0.3, 0.16), 'umean': (0, 200, 150), 'udiff': (-0.5, 0.5, 0), 'uoffs': (-0.5, 0.5, 0), 'testFL': (-10,10,0), 'testFR': (-10,10,0)}
# Params using pybullet GUI (sliders)
dbgIDs = {k : p.addUserDebugParameter(k, *params[k]) for k in params.keys()}
def P(k):
    try:
        return p.readUserDebugParameter(dbgIDs[k])
    except: # if in p.DIRECT mode, just return the default
        return params[k][-1]

def controller(t, q, dq):
    # Stroke kinematics
    omega = 2 * np.pi * P('freq') #ctrl['freq'] #
    ph = omega * bee.simt
    # force control
    umean = P('umean')
    udiff = P('udiff')
    uoffs = P('uoffs')
    return np.array([1 + udiff, 1 - udiff]) * umean * (np.sin(ph) + uoffs)

# --- Actual simulation ---
while True:
    try:
        # actual sim
        ss = bee.sampleStates()

        pdes = np.zeros(6)
        tau = controller(*ss)
        data = viewlog.appendLog(data, *ss, tau, pdes) # log
        
        bee.update(tau)#, testF=[P('testFL'), P('testFR')])
        time.sleep(bee._slowDown * bee.TIMESTEP * 1e-3)
    except:
        viewlog.saveLog('../logs/sdab', data)
        raise

