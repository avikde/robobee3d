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

# TEST
import sys
import scipy.linalg

# For x = (q position, p momentum)
Md = np.array([100, 100, 100, 3333, 3333, 1000])
M = np.diag(Md)
T0 = 9.81e-3 # current thrust
Z6 = np.zeros((6,6))
Z3 = np.zeros((3,3))
dpdotdphi = Md[0] * T0 * np.array([[0,1,0],
        [-1,0,0],
        [0,0,0]])
dy2dy = np.block([[Z3, dpdotdphi],
    [Z3, Z3]])

B = np.vstack((
    np.zeros((6,4)), 
    np.array([
        [0, 0, 0, 0],
        [0, 0, 0, 0],
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ])
    ))
A = np.vstack((
    np.hstack((Z6, np.linalg.inv(M))),
    np.hstack((Z6, dy2dy))
))
Q = np.diag(np.hstack((np.ones(6), np.array([1,1,1,0.1,0.1,0.1]))))
R = np.eye(4)
# P = scipy.linalg.solve_continuous_are(A, B, Q, R)
# Solved in mma
S = np.array([
    [148.63323830090116,0,0,0,1227.1305940665584,0,60.459197644011745,0,0,0,0.9999999997776448,0],
    [0,148.63323829662758,0,-1227.130592601852,0,2.78655025145515e-8,0,60.45919761098977,0,-1.0000000000976244,0,0],
    [0,0,100.99504938362075,0,0,0,0,0,1.,0,0,0],
    [0,-1227.130592578972,0,67852.6374376914,0,8.025985844048187e-7,0,-1790.593937328141,0,59.25001364938505,0,8.003670131039207e-9],
    [1227.1305926378923,0,0,0,67852.63756961962,0,1790.5939392762568, 0,0,0,59.25001366836082,0],
    [0,0,0,0,0,319.37438786323156,0,0,0,0,0,1.0000000001617644],
    [60.459197634119256,0,0,0,1790.593940128918,0,77.34072255317956,0,0,0,1.4863323826787476,0],
    [0,60.45919760722755,0,-1790.5939375232454,0,1.341900666357133e-8,0,77.34072249089046,0,-1.486332382820918,0,0],
    [0,0,0.9999999999999662,0,0,0,0,0,1.0099504938362078,0,0,0],
    [0,-1.0000000000850744,0,59.250013648479246,0,-3.7758848664609356e-10,0,-1.486332382600924,0,0.36817599534229906,0,0],
    [1.0000000001838987,0,0,0,59.250013701244946,0,1.4863323832950708,0,0,0,0.3681759956855034,0],
    [0,0,0,0,0,0.9999999997637417,0,0,0,0,0,0.31937438779944344]])
aa = A.T @ S + S @ A - S @ B @ np.linalg.inv(R) @ B.T @ S + Q
# Cc = np.hstack((B, A @ B, A @ A @ B))
print(aa)
sys.exit()

# filtfreq is for the body velocity filter
bee = robobee.RobobeeSim(p.DIRECT if args.direct else p.GUI, slowDown=0, camLock=True, timestep=0.1, gui=0, filtfreq=0.16)
# load robot
startPos = [0,0,100]
startOrientation = p.getQuaternionFromEuler(np.zeros(3))#[0.5,-0.5,0])
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
        data = viewlog.appendLog(data, *ss, np.hstack((tau, controller.u4)), controller.pdes) # log
        
        bee.update(tau)#, testF=[P('testFL'), P('testFR')])

        if not args.direct:
            time.sleep(bee._slowDown * bee.TIMESTEP * 1e-3)
finally:
    viewlog.saveLog('../logs/sdab', data)

