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

print(olAvgWrenchLeft(120, 0, 0.16))
