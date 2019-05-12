import numpy as np
# trying to not rely on pybullet
from scipy.spatial.transform import Rotation
from controlutils.py.model import Model

g = 9.81

class ThrustStrokeDev(Model):
    '''4 inputs: left/right thrusts and each has a stroke deviation in the x direction. There is a fixed offset of the thrusts in the y direction (~COP offset along wing).'''
    m = 0.5
    ib = 0.001

    def dynamics(self, y, u):
        nq = 3
        u1 = u[0]
        u2 = u[1]
        ydot = np.hstack((y[nq:], np.array([
            -u1 / self.m * np.sin(y[2]),
            u1 / self.m * np.cos(y[2]) - g,
            u1 * u2 / self.ib
        ])))
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
