import time, subprocess
import numpy as np
import pybullet as p
import robobee
np.set_printoptions(precision=2, suppress=True, linewidth=200)

# Usage params
STROKE_FORCE_CONTROL = True # if false, use position control on the stroke

bee = robobee.RobobeeSim(slowDown=0.1, camLock=True, timestep=0.1)
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
ctrl = {'thrust': 40, 'strokedev': 0, 'ampl': np.pi/4, 'freq': 0.175}
tLastPrint = 0

idf = p.addUserDebugParameter("freq", 0, 0.3, 0.1)
idkh = p.addUserDebugParameter("khinge", 0, 0.1, 0.01)
idbh = p.addUserDebugParameter("bhinge", 0, 0.1, 0.02)
idrc = p.addUserDebugParameter("rcopnondim", 0, 2, 0.5)
idff = p.addUserDebugParameter("ff", -1, 1, 0)
idff2 = p.addUserDebugParameter("ff2", -1, 1, 0)

while True:
    try:
        # actual sim
        bee.sampleStates()
        robobee.rcopnondim = p.readUserDebugParameter(idrc)
        # p.setJointMotorControlArray(bid, [1,3], p.PD_CONTROL, targetPositions=[0,0], positionGains=p.readUserDebugParameter(idkh)*np.ones(2), velocityGains=p.readUserDebugParameter(idbh)*np.ones(2))
        p.setJointMotorControlArray(bid, [1,3], p.POSITION_CONTROL, targetPositions=[0,0], positionGains=p.readUserDebugParameter(idkh)*np.ones(2), velocityGains=p.readUserDebugParameter(idbh)*np.ones(2))

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

        # Stroke kinematics (not force controlled yet)
        omega = 2 * np.pi * p.readUserDebugParameter(idf) #ctrl['freq'] #
        ph = omega * bee.simt
        th0 = ctrl['ampl'] * (np.sin(ph) + ctrl['strokedev'])
        dth0 = omega * ctrl['ampl'] * np.cos(ph)
        posdes = [0.3,0.3]#[th0,th0]

        # OR force control
        tau = np.full(2, ctrl['thrust'] * (np.sin(ph) + ctrl['strokedev']))

        # pass tau or posdes
        bee.update(posdes, testF=[p.readUserDebugParameter(idff), p.readUserDebugParameter(idff2)], forceControl=False)
        # bee.update(tau, forceControl=True)
        time.sleep(bee._slowDown * bee.TIMESTEP)
    except:
        raise

