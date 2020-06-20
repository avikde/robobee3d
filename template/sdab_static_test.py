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
bid = bee.load("../urdf/sdab.urdf", startPos, startOrientation, useFixedBase=False)

def olSweepAndResult(Vamp, f, tendMS=100, h2=0, h3=0):
    print(Vamp, f)
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
    # scale
    amps = np.rad2deg(np.ptp(qw, axis=0)) # get peak to peak
    amps[[1,3]] *= 0.5 # for hinge only want amplitude
    amps[[0,2]] /= Vamp # voltage-normalized
    return amps

fs = np.linspace(0.1, 0.2, num=15)
res1 = np.array([olSweepAndResult(30, f) for f in fs])
res2 = np.array([olSweepAndResult(40, f) for f in fs])

fig, ax = plt.subplots(2)

fs *= 1e3 # to Hz
ax[0].plot(fs, res1[:,0])
ax[0].plot(fs, res2[:,0])
ax[0].set_ylabel('Norm stroke [deg/V]')
ax[1].plot(fs, res1[:,1])
ax[1].plot(fs, res2[:,1])
ax[1].set_ylabel('Pitch amplitude [deg]')
ax[-1].set_xlabel('Freq [Hz]')

plt.show()