import subprocess, sys, progressbar
import numpy as np
from scipy.interpolate import SmoothBivariateSpline
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pybullet as p
import robobee

def sweepFile(fname, Vmeans, uoffss, fs, udiffs, h2s):
    # p.DIRECT for non-graphical
    bee = robobee.RobobeeSim(p.DIRECT, slowDown=1, camLock=True, timestep=0.1, gui=0)
    # load robot
    startPos = [0,0,10]
    startOrientation = p.getQuaternionFromEuler(np.zeros(3))
    subprocess.call(["python", "../urdf/xacro.py", "../urdf/sdab.xacro", "-o", "../urdf/sdab.urdf"])
    bid = bee.load("../urdf/sdab.urdf", startPos, startOrientation, useFixedBase=True)

    # create a progress bar
    nrun = 0
    widgets = [
        'Progress: ', progressbar.Percentage(),
        ' ', progressbar.Bar(),
        ' ', progressbar.ETA(),
    ]
    bar = progressbar.ProgressBar(widgets=widgets, max_value=len(Vmeans)*len(uoffss)*len(fs)*len(udiffs)*len(h2s))

    def olAvgWrench(Vmean, uoffs, f, udiff, h2):
        """Incorporate both wings"""
        nonlocal nrun
        nrun += 1
        bar.update(nrun)
        qw = bee.openLoop(Vmean * (1 + udiff), Vmean * (1 - udiff), uoffs, f, h2=h2, verbose=False)
        # avg wrench
        NptsInPeriod = int(1/(f * bee.TIMESTEP))
        return np.mean(qw[-NptsInPeriod:,-6:], axis=0)

    res = np.array([
        np.hstack((Vmean, uoffs, f, udiff, h2, 
        olAvgWrench(Vmean, uoffs, f, udiff, h2))) 
        for Vmean in Vmeans for uoffs in uoffss for f in fs for udiff in udiffs for h2 in h2s])
    with open('test.npy', 'wb') as f:
        np.save(f, res)

unpackDat = lambda dat : (dat[:,0], dat[:,1], dat[:,2], dat[:,3], dat[:,4], dat[:,5:])

def splineContour(ax, xiu, yiu, Zfun, length=50, dx=0, dy=0):
    dimrow = lambda row : np.linspace(row.min(), row.max(), length)
    xi = dimrow(xiu)
    yi = dimrow(yiu)
    zi = Zfun(xi, yi, grid=True, dx=dx, dy=dy)
    return ax.contourf(xi, yi, zi.T, cmap='RdBu')

if __name__ == "__main__":
    np.set_printoptions(precision=2, suppress=True, linewidth=200)

    if len(sys.argv) > 1:
        with open(sys.argv[1], 'rb') as f:
            dat = np.load(f)
        Vmeans, uoffss, fs, udiffs, h2s, ws = unpackDat(dat)

        # scatter vis
        fig = plt.figure()
        ax = fig.add_subplot(221, projection='3d')
        ax.plot(Vmeans, uoffss, ws[:,2], '.')
        ax.set_xlabel('Vmean')
        ax.set_ylabel('uoffs')
        ax.set_zlabel('Fz')
        ax = fig.add_subplot(222, projection='3d')
        ax.plot(Vmeans, uoffss, ws[:,4], '.')
        ax.set_xlabel('Vmean')
        ax.set_ylabel('uoffs')
        ax.set_zlabel('Ry')

        # test spline
        Fzfun = SmoothBivariateSpline(Vmeans, uoffss, ws[:,2])
        Ryfun = SmoothBivariateSpline(Vmeans, uoffss, ws[:,4])
        # Fzs = Fzfun(np.linspace(), uoffss, grid=False)
        # dFzs = Fzfun(Vmeans, uoffss, dx=1, dy=1, grid=False)
        # Rys = Ryfun(Vmeans, uoffss, grid=False)

        ax = fig.add_subplot(223)#, projection='3d')
        c = splineContour(ax, Vmeans, uoffss, Fzfun)
        fig.colorbar(c, ax=ax)
        ax.set_xlabel('Vmean')
        ax.set_ylabel('uoffs')
        # ax.set_zlabel('Fz')
        ax = fig.add_subplot(224)
        c = splineContour(ax, Vmeans, uoffss, Ryfun)
        fig.colorbar(c, ax=ax)
        ax.set_xlabel('Vmean')
        ax.set_ylabel('uoffs')
        # ax.set_zlabel('dFz')
        plt.show()

    else:
        # save to file
        Vmeans = np.linspace(120, 160, num=10)
        uoffss = np.linspace(-0.5, 0.5, num=10)
        fs = [0.165]
        udiffs = np.linspace(-0.2, 0.2, num=10)
        h2s = [0.]
        sweepFile('test.npy', Vmeans, uoffss, fs, udiffs, h2s)
