import time, subprocess, argparse
import numpy as np
import pybullet as p
import robobee
from robobee_test_controllers import OpenLoop, WaypointHover
import viewlog
np.set_printoptions(precision=4, suppress=True, linewidth=200)

# Generate URDF
subprocess.call(["python", "../urdf/xacro.py", "../urdf/sdab.xacro", "-o", "../urdf/sdab.urdf"])

def runSim(poptsFile, direct, tend, useMPC=True):
    # filtfreq is for the body velocity filter
    bee = robobee.RobobeeSim(p.DIRECT if direct else p.GUI, slowDown=0, camLock=True, timestep=0.1, gui=0, filtfreq=0.16)
    # load robot
    startPos = [0,0,100]
    startOrientation = p.getQuaternionFromEuler(np.zeros(3))#[0.5,-0.5,0])#
    bid = bee.load("../urdf/sdab.urdf", startPos, startOrientation, useFixedBase=False)
    data = viewlog.initLog()
    # controller = OpenLoop()
    controller = WaypointHover(poptsFile, useh2=False, useMPC=useMPC)#, constPdes=[0.,0,10,0,0,0])

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

def papPlots(poptsFile, tend=1000):
    import matplotlib.pyplot as plt
    l1 = runSim(poptsFile, True, tend, useMPC=True)
    l2 = runSim(poptsFile, True, tend, useMPC=False)
    qb = lambda data: data['q'][:,-7:]
    # dqb = data['dq'][:,-6:]
    # viewlog.defaultPlots(l1)
    fig, ax = plt.subplots(3,2)
    ax = ax.ravel()
    for i in range(3):
        ax[i].plot(l1['t'], qb(l1)[:,i], 'b')
        ax[i].plot(l2['t'], qb(l2)[:,i], 'r')
        ax[i].plot(l1['t'], l1['posdes'][:,i], 'k--', alpha=0.3)
    
    ax[3].plot(l1['t'], l1['u'][:,2], 'b') # Vmean
    ax[3].plot(l1['t'], l2['u'][:,2], 'r') # Vmean
    ax[3].set_ylabel('Vmean')
    ax[4].plot(l1['t'], l1['u'][:,3], 'b')
    ax[4].plot(l1['t'], l2['u'][:,3], 'r')
    ax[4].set_ylabel('uoffs')
    ax[5].plot(l1['t'], l1['u'][:,4], 'b')
    ax[5].plot(l1['t'], l2['u'][:,4], 'r')
    ax[5].set_ylabel('diff')
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('poptsFile', nargs='?', default='popts.npy')
    parser.add_argument('-t', '--tend', type=float, default=np.inf, help='end time [ms]')
    parser.add_argument('-d', '--direct', action='store_true', default=False, help='direct mode (no visualization)')
    args = parser.parse_args()
    
    # runSim(args.poptsFile, args.direct, args.tend)

    papPlots(args.poptsFile, tend=300)
