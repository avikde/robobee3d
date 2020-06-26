import time, subprocess
import numpy as np
import pybullet as p
import robobee
from robobee_test_controllers import OpenLoop, SimpleHover
import viewlog
np.set_printoptions(precision=2, suppress=True, linewidth=200)

# p.DIRECT for non-graphical
bee = robobee.RobobeeSim(p.GUI, slowDown=1, camLock=True, timestep=0.2, gui=0)
# load robot
startPos = [0,0,10]
startOrientation = p.getQuaternionFromEuler([0.5,0.5,0])
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

# controller = OpenLoop()
controller = SimpleHover()

# --- Actual simulation ---
while True:
    try:
        # actual sim
        ss = bee.sampleStates()

        pdes = np.zeros(6)
        tau = controller.update(*ss)
        data = viewlog.appendLog(data, *ss, tau, pdes) # log
        
        bee.update(tau)#, testF=[P('testFL'), P('testFR')])
        time.sleep(bee._slowDown * bee.TIMESTEP * 1e-3)
    except:
        viewlog.saveLog('../logs/sdab', data)
        raise

