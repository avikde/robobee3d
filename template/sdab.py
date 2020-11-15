import time, subprocess, argparse
import numpy as np
import pybullet as p
import robobee
from robobee_test_controllers import OpenLoop, WaypointHover
import viewlog
import matplotlib.pyplot as plt
np.set_printoptions(precision=4, suppress=True, linewidth=200)

# Generate URDF
subprocess.call(["python", "../urdf/xacro.py", "../urdf/sdab.xacro", "-o", "../urdf/sdab.urdf"])

def runSim(poptsFile, direct, tend, **kwargs):
    # filtfreq is for the body velocity filter
    bee = robobee.RobobeeSim(p.DIRECT if direct else p.GUI, slowDown=0, camLock=True, timestep=0.1, gui=0, filtfreq=0.16)
    # load robot
    startPos = [0,0,100]
    startOrientation = p.getQuaternionFromEuler(np.zeros(3))#[0.5,-0.5,0])#
    bid = bee.load("../urdf/sdab.urdf", startPos, startOrientation, useFixedBase=False)
    data = viewlog.initLog()
    # controller = OpenLoop()
    controller = WaypointHover(poptsFile, startPos, useh2=False, **kwargs)#, constPdes=[0.,0,10,0,0,0])

    # --- Actual simulation ---
    try:
        while bee.simt < tend:
            # actual sim
            ss = bee.sampleStates()

            tau = controller.update(*ss)
            # Also log the 4-dim u
            data = viewlog.appendLog(data, *ss, np.hstack((tau, controller.u4)), controller.accdes, controller.posdes) # log
            
            bee.update(tau)#, testF=[P('testFL'), P('testFR')])

            if not args.direct:
                time.sleep(bee._slowDown * bee.TIMESTEP * 1e-3)
    finally:
        viewlog.saveLog('../logs/sdab', data)
    return data

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('poptsFile', nargs='?', default='popts.npy')
    parser.add_argument('-t', '--tend', type=float, default=np.inf, help='end time [ms]')
    parser.add_argument('-d', '--direct', action='store_true', default=False, help='direct mode (no visualization)')
    parser.add_argument('-m', '--mpc', action='store_true', default=False, help='Use MPC')
    parser.add_argument('-w', '--wlqp', action='store_true', default=False, help='Use WLQP')
    parser.add_argument('-e', '--exp', default='helix', help='Task')
    args = parser.parse_args()
    
    log = runSim(args.poptsFile, args.direct, args.tend, useMPC=args.mpc, useWLQP=args.wlqp, task=args.exp)
    if args.direct:
        viewlog.defaultPlots(log)
        plt.show()

