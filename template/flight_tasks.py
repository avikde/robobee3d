"""Tasks specified as desired pos, vel and nominal traj"""
import numpy as np

VERTICAL = np.array([0,0,1])

def helix(t, initialPos, trajAmp=80, trajFreq=1, dz=0.15, useY=True):
    """If trajAmp is small, can be used for hover"""
    posdes = np.copy(initialPos)
    dpdes = np.zeros(3)
    trajOmg = 2 * np.pi * trajFreq * 1e-3 # to KHz, then to rad/ms
    posdes[0] += trajAmp * np.sin(trajOmg * t)
    dpdes[0] = trajAmp * trajOmg * np.cos(trajOmg * t)
    if useY:
        posdes[1] += trajAmp * (1 - np.cos(trajOmg * t))
        dpdes[1] = trajAmp * trajOmg * np.sin(trajOmg * t)
    if trajAmp > 1e-3: # otherwise assume hover
        posdes[2] += dz*t
        dpdes[2] = dz

    return posdes, dpdes, VERTICAL

def straightAcc(t, initialPos, tduration=500, vdes=2):
    """duration in ms, vdes in m/s"""
    pdes = np.copy(initialPos)
    dpdes = np.zeros(3)
    dpdes[0] = vdes if t < tduration else 0
    pdes[0] += vdes * np.clip(t, 0, tduration)
    
    return pdes, dpdes, VERTICAL
