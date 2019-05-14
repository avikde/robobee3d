import autograd.numpy as np
# trying to not rely on pybullet
from scipy.spatial.transform import Rotation
from controlutils.py.model import Model
import controlutils.py.kinematics as kinematics

g = 9.81

class ThrustStrokeDev(Model):
    '''4 inputs: left/right thrusts and each has a stroke deviation in the x direction. There is a fixed offset of the thrusts in the y direction (~COP offset along wing).'''
    m = 0.5
    # FIXME:
    Ib = np.diag([0.001, 0.001, 0.0001])
    ycp = 0.05

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
    

class QuasiSteadySDAB:
    # Model with force control of the wing spar
    # Could make it modular so that lumped actuator models can be introduced as well

    # Good parameters
    CD0 = 0.4
    CDmax = 3.4
    CLmax = 1.8

    # Modeling choices

    # True in Chen (2017) science robotics, but differently calculated in Osborne (1951)
    BODY_FRAME_FIXED_LIFT_DIRECTION = True
    RHO = 1.225 # density of air kg/m^3
    AERO_REGULARIZE_EPS = 1e-10 # stops undefined AoA when no wind
    
    def __init__(self, urdfParams):
        self.d = urdfParams['d']
        self.ycp = urdfParams['rcp']
        self.cbar = urdfParams['cbar']
    
    def CF(self, a):
        # in order lift,drag
        return np.array([self.CLmax * np.sin(2*a), (self.CDmax + self.CD0) / 2.0 - (self.CDmax - self.CD0) / 2.0 * np.cos(2*a)])

    def aerodynamics(self, q, dq, lrSign, worldFrame=True):
        # pass in current positions 11DOF: joints (4) + com (3 linear + 4 angular)

        # vector from center to distance along spar
        # joint angles
        if lrSign > 0:
            theta = -q[2:4]
            dtheta = -dq[2:4]
        else:
            theta = q[0:2]
            dtheta = dq[0:2]

        # These are all in the body frame
        Rspar = Rotation.from_euler('z', theta[0])
        Rhinge = Rotation.from_euler('y', -lrSign * theta[1])
        # vector along wing
        sparVecB = Rspar.apply(np.array([0, lrSign, 0]))
        # wing chord unit vector
        chordB = Rspar.apply(Rhinge.apply(np.array([0,0,-1])))
        # TODO: add external wind vector
        # NOTE: assuming spar velocity is dominant
        wB = -np.cross(dtheta[0] * np.array([0,0,1]), self.ycp * sparVecB)
        # print(dtheta[0] * np.array([0,0,1]), self.ycp * sparVecB)

        # COP: half od cbar down
        pcopB = np.array([0,0,self.d]) + self.ycp * sparVecB + 0.5 * self.cbar * chordB

        # Various directions in the notation of Osborne (1951)
        # l = vector along wing
        # w = relative wind vector
        # c = chord
        lwB = np.cross(sparVecB, wB)
        lwpB = wB - wB.dot(sparVecB) * sparVecB
        wnorm = np.sqrt(wB.dot(wB)) + self.AERO_REGULARIZE_EPS
        lwnorm = np.sqrt(lwB.dot(lwB)) + self.AERO_REGULARIZE_EPS
        lwpnorm = np.sqrt(lwpB.dot(lwpB)) + self.AERO_REGULARIZE_EPS

        # Lift/drag directions
        eD = lwpB / lwpnorm
        if self.BODY_FRAME_FIXED_LIFT_DIRECTION:
            eL = np.array([0,0,1])
        else:
            # FIXME: needs some reversal for one half-stroke
            eL = lwB / lwnorm
            raise 'Not implemented fully'

        # Calculate aero force
        aoa = np.arccos(chordB.dot(wB) / wnorm)
        Cf = self.CF(aoa)
        # Cf *= 0.5 * rho * beta
        FaeroB = 0.5 * self.RHO * self.cbar * self.ycp * (Cf[0] * eL + Cf[1] * eD) * lwnorm**2

        # Body to world frame --
        pcom = q[4:7]
        Rb = Rotation.from_quat(q[7:11]) # scalar-last format
        if worldFrame:
            pcopW = pcom + Rb.apply(pcopB)
            FaeroW = Rb.apply(FaeroB)
            # for external torque about wing hinge, use r X F
            hingeTorque = np.cross(0.5 * self.cbar * Rb.apply(chordB), FaeroW)
            return pcopW, FaeroW, hingeTorque
        else:
            # body frame
            hingeTorque = np.cross(0.5 * self.cbar * chordB, FaeroB)
            return pcopB, FaeroB, hingeTorque
