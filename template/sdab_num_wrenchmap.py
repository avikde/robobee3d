import subprocess, sys
import numpy as np
import matplotlib.pyplot as plt
import pybullet as p
import robobee

def sweepFile(fname, Vmeans, uoffss, fs, udiffs, h2s):
    # p.DIRECT for non-graphical
    bee = robobee.RobobeeSim(p.DIRECT, slowDown=1, camLock=True, timestep=0.1, gui=0)
    # load robot
    startPos = [0,0,10]
    startOrientation = p.getQuaternionFromEuler(np.zeros(3))
    subprocess.call(["python", "../urdf/xacro.py", "../urdf/sdab.xacro", "-o", "../urdf/sdab.urdf"])
    bid = bee.load("../urdf/sdab.urdf", startPos, startOrientation, useFixedBase=True)

    def olAvgWrench(Vmean, uoffs, f, udiff, h2):
        """Incorporate both wings"""
        qw = bee.openLoop(Vmean * (1 + udiff), Vmean * (1 - udiff), uoffs, f, h2=h2)
        # avg wrench
        NptsInPeriod = int(1/(f * bee.TIMESTEP))
        return np.mean(qw[-NptsInPeriod:,-6:], axis=0)

    res = np.array([
        np.hstack((Vmean, uoffs, f, udiff, h2, 
        olAvgWrench(Vmean, uoffs, f, udiff, h2))) 
        for Vmean in Vmeans for uoffs in uoffss for f in fs for udiff in udiffs for h2 in h2s])
    with open('test.npy', 'wb') as f:
        np.save(f, res)

if __name__ == "__main__":
    np.set_printoptions(precision=2, suppress=True, linewidth=200)

    if len(sys.argv) > 1:
        with open(sys.argv[1], 'rb') as f:
            dat = np.load(f)
        print(dat)
    else:
        # save to file
        Vmeans = np.linspace(120, 160, num=2)
        uoffss = np.linspace(-0.5, 0.5, num=1)
        fs = [0.15]
        udiffs = [0.]
        h2s = [0.]
        sweepFile('test.npy', Vmeans, uoffss, fs, udiffs, h2s)
