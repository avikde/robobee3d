import time, subprocess
import numpy as np
import pybullet as p
import robobee
np.set_printoptions(precision=2, suppress=True, linewidth=200)

# p.DIRECT for non-graphical
bee = robobee.RobobeeSim(p.DIRECT, slowDown=1, camLock=True, timestep=0.2, gui=1)
# load robot
startPos = [0,0,10]
startOrientation = p.getQuaternionFromEuler(np.zeros(3))
subprocess.call(["python", "../urdf/xacro.py", "../urdf/sdab.xacro", "-o", "../urdf/sdab.urdf"])
bid = bee.load("../urdf/sdab.urdf", startPos, startOrientation, useFixedBase=False)

def olSweepAndResult(Vamp, f, tendMS=100, h2=0, h3=0):
    qw = []
    omega = 2 * np.pi * f

    while bee.simt < tendMS:
        ss = bee.sampleStates()
        ph = omega * bee.simt
        V = Vamp * ((1+h3)*np.sin(ph) + h3*np.sin(3*ph) + h2*np.sin(2*ph))
        tau = [V,V] # same to both wings
        qw.append(np.copy(bee.q[:4]))
        bee.update(tau)
    
    bee.reset()
    qw = np.array(qw)
    return np.ptp(qw, axis=0) # get peak to peak

a = olSweepAndResult(50, 0.15)
print(a)
