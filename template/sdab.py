import time, subprocess, argparse
import numpy as np
import pybullet as p
import robobee
from robobee_test_controllers import OpenLoop, WaypointHover
import viewlog
np.set_printoptions(precision=4, suppress=True, linewidth=200)
import matplotlib.pyplot as plt
from plot_helpers import *

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

def papExps(task, poptsFile, tend=1000):
    def doTask(s):
        l1 = runSim(poptsFile, True, tend, useMPC=True, task=s)
        l2 = runSim(poptsFile, True, tend, useMPC=False, task=s)
    doTask(task)

def papPlots(fmpc, freac, vscale):
    l1 = viewlog.readFile(fmpc)
    l2 = viewlog.readFile(freac)

    fig = plt.figure()
    ax3d = fig.add_subplot(1,1,1,projection='3d')
    traj3plot(ax3d, l1['t'], l1['p'], l1['s'], "Blues_r", vscale=vscale)
    aspectEqual3(ax3d, l1['p'])
    traj3plot(ax3d, l2['t'], l2['p'], l2['s'], "Reds_r", vscale=vscale)
    ax3d.plot(l1['posdes'][:,0], l1['posdes'][:,1], l1['posdes'][:,2], 'k--', alpha=0.5, zorder=9)
    ax3d.set_xlabel('x [mm]')
    ax3d.set_ylabel('y [mm]')
    ax3d.set_zlabel('z [mm]')
    
    fig, ax = plt.subplots(1,3,figsize=(7.5,2.5))
    ax=ax.ravel()
    ax[0].plot(l1['t'], l1['u'][:,2], 'b') # Vmean
    ax[0].plot(l1['t'], l2['u'][:,2], 'r') # Vmean
    ax[0].set_ylabel('Vmean')
    ax[1].plot(l1['t'], l1['u'][:,3], 'b')
    ax[1].plot(l1['t'], l2['u'][:,3], 'r')
    ax[1].set_ylabel('uoffs')
    ax[2].plot(l1['t'], l1['u'][:,4], 'b')
    ax[2].plot(l1['t'], l2['u'][:,4], 'r')
    ax[2].set_ylabel('diff')
    fig.tight_layout()
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('poptsFile', nargs='?', default='popts.npy')
    parser.add_argument('-t', '--tend', type=float, default=np.inf, help='end time [ms]')
    parser.add_argument('-d', '--direct', action='store_true', default=False, help='direct mode (no visualization)')
    args = parser.parse_args()
    
    # runSim(args.poptsFile, args.direct, args.tend, useMPC=True)

    # papExps('helix', args.poptsFile, tend=3000)
    papPlots('../logs/sdab_20201112190409.zip', '../logs/sdab_20201112190440.zip', 40)
    # papExps('line', args.poptsFile, tend=1000)
    papPlots('../logs/sdab_20201112190644.zip', '../logs/sdab_20201112190654.zip', 50)
