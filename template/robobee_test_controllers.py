"""Each controller has a set of parameters and ranges as well as a function to evaluate it"""
import pybullet as p
import numpy as np

class RobobeeController(object):
    def __init__(self, params):
        self.params = params
        # Params using pybullet GUI (sliders)
        self.dbgIDs = {k : p.addUserDebugParameter(k, *params[k]) for k in params.keys()}

    def P(self, k):
        try:
            return p.readUserDebugParameter(self.dbgIDs[k])
        except: # if in p.DIRECT mode, just return the default
            return self.params[k][-1]
    
    def update(self, t, q, dq):
        raise NotImplementedError

class OpenLoop(RobobeeController):
    def update(self, t, q, dq):
        """Open-loop controller"""
        # Stroke kinematics
        omega = 2 * np.pi * self.P('freq') #ctrl['freq'] #
        ph = omega * t
        # force control
        umean = self.P('umean')
        udiff = self.P('udiff')
        uoffs = self.P('uoffs')
        return np.array([1 + udiff, 1 - udiff]) * umean * (np.sin(ph) + uoffs)
