import time, subprocess
import numpy as np
import matplotlib.pyplot as plt
import pybullet as p
import robobee
np.set_printoptions(precision=2, suppress=True, linewidth=200)

# p.DIRECT for non-graphical
bee = robobee.RobobeeSim(p.DIRECT, slowDown=1, camLock=True, timestep=0.1, gui=0)
# load robot
startPos = [0,0,10]
startOrientation = p.getQuaternionFromEuler(np.zeros(3))
subprocess.call(["python", "../urdf/xacro.py", "../urdf/sdab.xacro", "-o", "../urdf/sdab.urdf"])
bid = bee.load("../urdf/sdab.urdf", startPos, startOrientation, useFixedBase=True)

def olAvgWrenchLeft(Vamp, uoffs, f, **kwargs):
    qw = bee.openLoopLeft(Vamp, uoffs, f, **kwargs)
    # avg wrench
    NptsInPeriod = int(1/(f * bee.TIMESTEP))
    return np.mean(qw[-NptsInPeriod:,-6:], axis=0)

def olAvgWrench(Vmean, uoffs, f, udiff, h2):
    """Incorporate both wings"""
    wl = olAvgWrenchLeft(Vmean * (1 + udiff), uoffs, f, h2=h2)
    wr = olAvgWrenchLeft(Vmean * (1 - udiff), uoffs, f, h2=-h2)
    # Need to mirror for wr
    wr[1] = -wr[1] # y
    wr[3] = -wr[3] # rx
    wr[5] = -wr[5] # rz
    return wl + wr
    
print(olAvgWrench(120, 0, 0.16, 0., -0.2))
