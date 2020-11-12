"""Tasks specified as desired pos, vel and nominal traj"""
import numpy as np

def helix(t, initialPos, trajAmp=80, trajFreq=1, dz=0.15):
    posdes = np.copy(initialPos)
    dpdes = np.zeros(3)
    trajOmg = 2 * np.pi * trajFreq * 1e-3 # to KHz, then to rad/ms
    posdes[0] += trajAmp * np.sin(trajOmg * t)
    dpdes[0] = trajAmp * trajOmg * np.cos(trajOmg * t)
    posdes[1] += trajAmp * (1 - np.cos(trajOmg * t))
    dpdes[1] = trajAmp * trajOmg * np.sin(trajOmg * t)
    posdes[2] += dz*t
    dpdes[2] = dz

    sdes = np.array([0,0,1])

    return posdes, dpdes, sdes

