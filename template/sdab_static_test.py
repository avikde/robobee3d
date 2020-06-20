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
    """Drives the wings with an open-loop signal and the given amplitude and freq, and returns a 3-vector containing:
    - Voltage-normalized stroke amplitude [deg/V]
    - Wing pitch amplitude [deg]
    - Avg lift [mN]
    """
    print('Now testing:', np.array([Vamp, f]))
    qw = []
    omega = 2 * np.pi * f

    while bee.simt < tendMS:
        ss = bee.sampleStates()
        ph = omega * bee.simt
        V = Vamp * ((1+h3)*np.sin(ph) + h3*np.sin(3*ph) + h2*np.sin(2*ph))
        tau = [V,V] # same to both wings
        # call aerodynamics to calculate avg lift
        F, _ = robobee.aerodynamics(bee.q[0:2], bee.dq[0:2], 1, bee.urdfParams)
        qw.append(np.copy(np.hstack((bee.q[:2], F[2]))))
        bee.update(tau)
    
    bee.reset()
    qw = np.array(qw)
    # stroke and pitch amplitudes
    amps = np.rad2deg(np.ptp(qw[:,:2], axis=0)) # get peak to peak
    amps[1] *= 0.5 # for hinge only want amplitude
    amps[0] /= Vamp # voltage-normalized
    # avg lift
    NptsInPeriod = int(1/(f * bee.TIMESTEP))
    FLavg = np.mean(qw[-NptsInPeriod:,-1])

    return np.hstack((amps, FLavg))

Vamps = [120, 150, 180]
fs = np.linspace(0.1, 0.2, num=20)
res = [np.array([olSweepAndResult(V, f) for f in fs]) for V in Vamps]

fig, ax = plt.subplots(3)

fs *= 1e3 # to Hz
for i in range(3):
    ax[0].plot(fs, res[i][:,0], label=str(Vamps[i])+'V')
ax[0].set_ylabel('Norm stroke [deg/V]')
for i in range(3):
    ax[1].plot(fs, res[i][:,1])
ax[1].set_ylabel('Pitch amplitude [deg]')
for i in range(3):
    ax[2].plot(fs, 1000 / 9.81 * res[i][:,2])
ax[2].set_ylabel('Avg lift/wing [mg]')

ax[0].legend()
ax[-1].set_xlabel('Freq [Hz]')

plt.show()
