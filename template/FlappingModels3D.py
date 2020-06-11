import autograd.numpy as np
import sys
# trying to not rely on pybullet
from scipy.spatial.transform import Rotation
sys.path.append('..')
from controlutils.py.model import Model
import controlutils.py.kinematics as kinematics

g = 9.81

class ThrustStrokeDev(Model):
    '''4 inputs: left/right thrusts and each has a stroke deviation in the x direction. There is a fixed offset of the thrusts in the y direction (~COP offset along wing).'''
    m = 0.5
    # FIXME:
    Ib = np.diag([0.0005, 0.0005, 0.001])
    ycp = 0.5
    lwh = np.array([0.05,0.05,0.2])

    def dynamics(self, y, u):
        nq = 6
        # FIXME: does not work with autograd see https://github.com/avikde/controlutils/issues/10
        wRotb = Rotation.from_rotvec(y[3:6]).as_dcm()
        omega = y[9:]  # in world frame
        # Compute Fb
        FL = np.array([0, 0, u[0]])
        FR = np.array([0, 0, u[2]])
        # compute rb
        rL = np.array([u[1], self.ycp, 0])
        rR = np.array([u[3], -self.ycp, 0])
        # Compute dynamics (RHS of the newton euler eqs)
        mpdd = np.array([0, 0, -self.m * g]) + wRotb @ (FL + FR)
        Iomegadotb = np.cross(rL, FL) + np.cross(rR, FR) - np.cross(omega, self.Ib @ omega)
        omegadotb = np.linalg.inv(self.Ib) @ Iomegadotb
        omegadot = wRotb.T @ (omegadotb)
        # apply rotation
        # assemble vector
        ydot = np.hstack((y[nq:], mpdd / self.m, omegadot))
        return ydot
    
    def forcesW(self, y, u, Fscale=0.1):
        # Compute Fb
        FL = np.array([0, 0, u[0]])
        FR = np.array([0, 0, u[2]])
        # compute rb
        rL = np.array([u[1], self.ycp, 0])
        rR = np.array([u[3], -self.ycp, 0])
        wRotb = Rotation.from_rotvec(y[3:6]).as_dcm()
        return Fscale * wRotb @ FL, Fscale * wRotb @ FR, y[:3] + wRotb @ rL, y[:3] + wRotb @ rR

