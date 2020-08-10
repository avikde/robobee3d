"""
This file generates numkins.py by running repeated simulations at different inputs and saving the kinematics features and wrench map as well.
"""
import subprocess, sys, progressbar
import numpy as np
import pybullet as p
import robobee

# Generate the data from simulation-----------------------------------------------------------------

def sweepFile(fname, Vmeans, uoffss, fs, udiffs, h2s):
    # p.DIRECT for non-graphical
    bee = robobee.RobobeeSim(p.DIRECT, slowDown=1, camLock=True, timestep=0.1, gui=0)
    # load robot
    startPos = [0,0,10]
    startOrientation = p.getQuaternionFromEuler(np.zeros(3))
    subprocess.call(["python", "../urdf/xacro.py", "../urdf/sdab.xacro", "-o", "../urdf/sdab.urdf"])
    bid = bee.load("../urdf/sdab.urdf", startPos, startOrientation, useFixedBase=True)

    # create a progress bar
    nrun = 0
    widgets = [
        'Progress: ', progressbar.Percentage(),
        ' ', progressbar.Bar(),
        ' ', progressbar.ETA(),
    ]
    bar = progressbar.ProgressBar(widgets=widgets, max_value=len(Vmeans)*len(uoffss)*len(fs)*len(udiffs)*len(h2s))

    def olAvgWrenchKinFeat(Vmean, uoffs, f, udiff, h2):
        """Incorporate both wings"""
        nonlocal nrun
        nrun += 1
        bar.update(nrun)
        # qw below contains wing kinematics as well as the wrench
        sw = bee.openLoop(Vmean * (1 + udiff), Vmean * (1 - udiff), uoffs, f, h2=h2, verbose=False)
        NptsInPeriod = int(1/(f * bee.TIMESTEP))
        # avg wrench (stored from sim for calibrating the kinfeat -> wrench analytical prediction)
        avgwrench = np.mean(sw[-NptsInPeriod:,-6:], axis=0)
        # Kinematics features:
        # amplitudes
        qw = sw[-NptsInPeriod:,:4]
        dqw = sw[-NptsInPeriod:,4:8]
        # to estimate alpha (fraction of period for upstroke)
        ratioPosStrokeVel = np.sum(dqw[:,0] >= 0, axis=0) / NptsInPeriod
        # for each wing, get max and min amplitude (corresponding to upstroke and downstroke)
        kins = np.hstack([np.hstack((np.amax(qw[:,2*i:2*i+2], axis=0), -np.amin(qw[:,2*i:2*i+2], axis=0))) for i in range(2)]) # size 8
        kins = np.hstack((kins, ratioPosStrokeVel)) # size 9
        return np.hstack((avgwrench, kins))

    res = np.array([
        np.hstack((Vmean, uoffs, f, udiff, h2, 
        olAvgWrenchKinFeat(Vmean, uoffs, f, udiff, h2))) 
        for Vmean in Vmeans for uoffs in uoffss for f in fs for udiff in udiffs for h2 in h2s])
    with open(fname, 'wb') as f:
        np.save(f, res)

    # xdata = np.array([
    #     np.hstack((Vmean, uoffs, f, udiff, h2)) 
    #     for Vmean in Vmeans for uoffs in uoffss for f in fs for udiff in udiffs for h2 in h2s])
    # cpus = multiprocessing.cpu_count()
    # pool = multiprocessing.Pool(processes=cpus)
    # res = pool.map(test, xdata)
    # TODO: multiprocessing but need different copies of the sim too

if __name__ == "__main__":
    # save to file
    Vmeans = np.linspace(90, 160, num=8)
    uoffss = np.linspace(-0.5, 0.5, num=8)
    fs = [0.17]
    udiffs = np.linspace(-0.2, 0.2, num=8)
    h2s = np.linspace(-0.2, 0.2, num=8)#[0.]
    sweepFile('numkinsaa.npy', Vmeans, uoffss, fs, udiffs, h2s)
