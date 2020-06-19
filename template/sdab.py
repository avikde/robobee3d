import time, subprocess
import numpy as np
import pybullet as p
import robobee
np.set_printoptions(precision=2, suppress=True, linewidth=200)

# Usage params
STROKE_FORCE_CONTROL = True # if false, use position control on the stroke

bee = robobee.RobobeeSim(slowDown=100, camLock=True, timestep=0.1)
# load robot
startPos = [0,0,10]
startOrientation = p.getQuaternionFromEuler(np.zeros(3))
subprocess.call(["python", "../urdf/xacro.py", "../urdf/sdab.xacro", "-o", "../urdf/sdab.urdf"])
bid = bee.load("../urdf/sdab.urdf", startPos, startOrientation, useFixedBase=True)
# print(jointId, urdfParams)

# Helper function: traj to track
def traj(t):
    return startPos + np.array([30 * np.sin(0.002*np.pi*t), 0, 0.5 * t])

# draw traj
T_END=1000
tdraw = np.linspace(0, T_END, 20)
for ti in range(1, len(tdraw)):
    p.addUserDebugLine(traj(tdraw[ti-1]), traj(tdraw[ti]), lineColorRGB=[0,0,0], lifeTime=0)
    
# ---


# conventional controller params
# ctrl = {'thrust': 40, 'strokedev': 0, 'ampl': np.pi/4, 'freq': 0.175}
tLastPrint = 0

idfreq = p.addUserDebugParameter("freq", 0, 0.3, 0.15)
idmean = p.addUserDebugParameter("umean [mN]", 0, 100, 40)
iddiff = p.addUserDebugParameter("udiff [ ]", -0.5, 0.5, 0)
idoffs = p.addUserDebugParameter("uoffs [ ]", -0.5, 0.5, 0)	
idff1 = p.addUserDebugParameter("testFL", -10, 10, 0)	
idff2 = p.addUserDebugParameter("testFR", -10, 10, 0)

while True:
    try:
        # actual sim
        bee.sampleStates()

        # # Conventional controller
        # # posErr = sim.q[4:7] - traj(sim.simt)
        # # FIXME: simple path
        # xDes = 0.03 * np.sin(2*np.pi*(sim.q[6] / 0.5))
        # zdDes = 0.05
        # ctrl['ampl'] = 0.9 - 0.1 * (sim.dq[6] - zdDes)
        # desPitch = np.clip(-(200 * (sim.q[4] - xDes) + 0.1 * sim.dq[4]), -np.pi/4.0, np.pi/4.0)
        # curPitch = Rotation.from_quat(sim.q[7:11]).as_euler('xyz')[1]
        # pitchCtrl = 2 * (curPitch - desPitch) + 0.1 * sim.dq[8]
        # ctrl['strokedev'] = np.clip(pitchCtrl, -0.4, 0.4)

        # Stroke kinematics
        omega = 2 * np.pi * p.readUserDebugParameter(idfreq) #ctrl['freq'] #
        ph = omega * bee.simt
        # th0 = ctrl['ampl'] * (np.sin(ph) + ctrl['strokedev'])
        # dth0 = omega * ctrl['ampl'] * np.cos(ph)
        # posdes = [th0,th0]

        # OR force control
        umean = p.readUserDebugParameter(idmean)
        udiff = p.readUserDebugParameter(iddiff)
        uoffs = p.readUserDebugParameter(idoffs)
        tau = np.array([1 + udiff, 1 - udiff]) * umean * (np.sin(ph) + uoffs)
        # np.full(2, ctrl['thrust'] * (np.sin(ph) + ctrl['strokedev']))

        # pass tau or posdes
        # bee.update(posdes, forceControl=False)
        bee.update(tau, forceControl=True, testF=[p.readUserDebugParameter(idff1), p.readUserDebugParameter(idff2)])

        time.sleep(bee._slowDown * bee.TIMESTEP * 1e-3)
    except:
        raise

