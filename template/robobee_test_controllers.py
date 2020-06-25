"""Each controller has a set of parameters and ranges as well as a function to evaluate it"""
import pybullet as p
import numpy as np

class WaveformGenerator(object):
    """A simple waveform generator that keeps track of phase"""
    def __init__(self):
        self.ph = 0 # phase
        self.tprev = 0
    
    def update(self, t, f, h2=0., h3=0.1):
        self.ph += 2 * np.pi * f * (t - self.tprev)
        self.tprev = t
        return (1 + h3) * np.sin(self.ph) + h3 * np.sin(3 * self.ph) + h2 * np.sin(2 * self.ph)

class RobobeeController(object):
    def __init__(self, params):
        # Params stored as (min, max, default) tuples
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
    """Open-loop controller"""
    def __init__(self):
        super(OpenLoop, self).__init__({'freq': (0, 0.3, 0.16), 'umean': (0, 200, 150), 'udiff': (-0.5, 0.5, 0), 'uoffs': (-0.5, 0.5, 0), 'testFL': (-10,10,0), 'testFR': (-10,10,0)})
        self.wf = WaveformGenerator()

    def update(self, t, q, dq):
        # force control
        umean = self.P('umean')
        udiff = self.P('udiff')
        uoffs = self.P('uoffs')
        w = self.wf.update(t, self.P('freq'))
        return np.array([1 + udiff, 1 - udiff]) * umean * (w + uoffs)

# class SimpleHover(RobobeeController):
#     """Simplified version of Pakpong (2013)"""
#     def __init__(self):
#         super(SimpleHover, self).__init__({})
    
#     def update(self, t, q, dq):

