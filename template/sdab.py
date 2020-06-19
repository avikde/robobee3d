import time, subprocess
import numpy as np
import pybullet as p
import robobee
np.set_printoptions(precision=2, suppress=True, linewidth=200)


bee = robobee.RobobeeSim(slowDown=1, camLock=True, timestep=0.1, gui=1)
# load robot
startPos = [0,0,10]
startOrientation = p.getQuaternionFromEuler(np.zeros(3))
subprocess.call(["python", "../urdf/xacro.py", "../urdf/sdab.xacro", "-o", "../urdf/sdab.urdf"])
bid = bee.load("../urdf/sdab.urdf", startPos, startOrientation, useFixedBase=False)

# Helper function: traj to track
def traj(t):
    return startPos + np.array([30 * np.sin(0.002*np.pi*t), 0, 0.5 * t])

# draw traj
T_END=1000
tdraw = np.linspace(0, T_END, 20)
for ti in range(1, len(tdraw)):
    p.addUserDebugLine(traj(tdraw[ti-1]), traj(tdraw[ti]), lineColorRGB=[0,0,0], lifeTime=0)
    
# ---

idfreq = p.addUserDebugParameter("freq", 0, 0.3, 0.15)
idmean = p.addUserDebugParameter("umean [V]", 0, 100, 50)
iddiff = p.addUserDebugParameter("udiff [ ]", -0.5, 0.5, 0)
idoffs = p.addUserDebugParameter("uoffs [ ]", -0.5, 0.5, 0)	
idff1 = p.addUserDebugParameter("testFL", -10, 10, 0)	
idff2 = p.addUserDebugParameter("testFR", -10, 10, 0)

while True:
    try:
        # actual sim
        bee.sampleStates()
        
        # Stroke kinematics
        omega = 2 * np.pi * p.readUserDebugParameter(idfreq) #ctrl['freq'] #
        ph = omega * bee.simt
        # force control
        umean = p.readUserDebugParameter(idmean)
        udiff = p.readUserDebugParameter(iddiff)
        uoffs = p.readUserDebugParameter(idoffs)
        tau = np.array([1 + udiff, 1 - udiff]) * umean * (np.sin(ph) + uoffs)
        # np.full(2, ctrl['thrust'] * (np.sin(ph) + ctrl['strokedev']))

        # pass tau
        bee.update(tau)#, testF=[p.readUserDebugParameter(idff1), p.readUserDebugParameter(idff2)])

        time.sleep(bee._slowDown * bee.TIMESTEP * 1e-3)
    except:
        raise

