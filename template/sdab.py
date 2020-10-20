import time, subprocess, argparse
import numpy as np
import pybullet as p
import robobee
from robobee_test_controllers import OpenLoop, WaypointHover
import viewlog
np.set_printoptions(precision=4, suppress=True, linewidth=200)

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('poptsFile', nargs='?', default='popts.npy')
parser.add_argument('-t', '--tend', type=float, default=np.inf, help='end time [ms]')
parser.add_argument('-d', '--direct', action='store_true', default=False, help='direct mode (no visualization)')
args = parser.parse_args()

# filtfreq is for the body velocity filter
bee = robobee.RobobeeSim(p.DIRECT if args.direct else p.GUI, slowDown=0, camLock=True, timestep=0.1, gui=0, filtfreq=0.16)
# load robot
startPos = [0,0,100]
startOrientation = p.getQuaternionFromEuler(np.zeros(3))#[0.5,-0.5,0])#
subprocess.call(["python", "../urdf/xacro.py", "../urdf/sdab.xacro", "-o", "../urdf/sdab.urdf"])
bid = bee.load("../urdf/sdab.urdf", startPos, startOrientation, useFixedBase=False)
data = viewlog.initLog()

# Helper function: traj to track
def traj(t):
    ph = 1*2*np.pi*(1e-3*t)
    return startPos + np.array([50 * np.sin(ph), 50 * (1 - np.cos(ph)), 0.0 * t])

# draw traj
tdraw = np.linspace(0, 500 if np.isinf(args.tend) else args.tend, 20)
for ti in range(1, len(tdraw)):
    p.addUserDebugLine(traj(tdraw[ti-1]), traj(tdraw[ti]), lineColorRGB=[0,0,0], lifeTime=0)
    
# ---

# controller = OpenLoop()
controller = WaypointHover(args.poptsFile, useh2=False)#, constPdes=[0.,0,10,0,0,0])

# --- Actual simulation ---
try:
    while bee.simt < args.tend:
        # actual sim
        ss = bee.sampleStates()

        controller.posdes = traj(bee.simt)
        tau = controller.update(*ss)
        # Also log the 4-dim u
        data = viewlog.appendLog(data, *ss, np.hstack((tau, controller.u4)), controller.accdes, controller.posdes) # log
        
        bee.update(tau)#, testF=[P('testFL'), P('testFR')])

        if not args.direct:
            time.sleep(bee._slowDown * bee.TIMESTEP * 1e-3)
finally:
    viewlog.saveLog('../logs/sdab', data)

