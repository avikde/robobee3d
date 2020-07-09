import subprocess, sys, progressbar
import numpy as np
from scipy.interpolate import SmoothBivariateSpline
from scipy.optimize import curve_fit
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
    if isinstance(Zfun, SmoothBivariateSpline):
        zi = (Zfun(xi, yi, grid=True, dx=dx, dy=dy)).T
    else:
        # create a grid for the plot
        Xi, Yi = np.meshgrid(xi, yi)
        positions = np.vstack([Xi.ravel(), Yi.ravel()])
        zi = np.reshape(Zfun(positions.T), Xi.shape)
        # zi = Zfun(Xi, Yi, np.zeros_like(Xi))
    return ax.contourf(xi, yi, zi, cmap='RdBu')
    
class FunApprox:
    """Using https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html"""
    def __init__(self, k):
        # Rk to R1
        self.k = k
        self.A2 = np.zeros((self.k, self.k))
        self.xs, self.ys = np.triu_indices(self.k)
    
    def nparams(self):
        return int(1 + self.k + (self.k * (self.k+1))/2)
    
    def unpackp(self, p):
        a0 = p[0]
        a1 = np.array(p[1:self.k+1])
        a2 = np.array(p[self.k+1:])
        # create a symmetric matrix
        self.A2[self.xs, self.ys] = a2
        self.A2[self.ys, self.xs] = a2
        return a0, a1

    def f(self, xdata, *p):
        a0, a1 = self.unpackp(p)
        # print(a.shape, b, xdata.shape)
        y = a0 + xdata @ a1
        for i in range(xdata.shape[0]):
            xi = xdata[i,:]
            y[i] += xi.T @ self.A2 @ xi
        # print(y.shape)
        return y
    
    def df_dx(self, xi, *p):
        """This is dy/dx; NOT the jacobian wrt params required by the curve fit alg"""
        a0, a1 = self.unpackp(p)
        # y = a0 + a1 * x + x^T * A2 * x
        return a1 + self.A2 @ xi

def wrenchMap(xdata, popts):
    """popts = (6,k)-shaped array of optimized params for each wrench component.
    xdata = N,Nu
    Returns N,6"""
    Np = xdata.shape[0]
    return np.vstack([fa.f(xdata, *popts[i,:]) for i in range(6)]).T

if __name__ == "__main__":
    np.set_printoptions(precision=2, suppress=True, linewidth=200)

    if len(sys.argv) > 1:
        with open(sys.argv[1], 'rb') as f:
            dat = np.load(f)
        Vmeans, uoffss, fs, udiffs, h2s, ws = unpackDat(dat)
        
        fa = FunApprox(4) # k
        xdata = np.vstack((Vmeans, uoffss, udiffs, h2s)).T # k,M

        # Optimized param fits in each row for each component of the wrench
        popts = np.vstack([curve_fit(fa.f, xdata, ws[:,i], p0=np.ones(fa.nparams()))[0] for i in range(6)])

        def fitWi(i, ax3d, ax):
            # scatter
            ax3d.plot(Vmeans, uoffss, ws[:,i], '.')
            ax3d.set_xlabel('Vmean')
            ax3d.set_ylabel('uoffs')
            ax3d.set_zlabel('W'+str(i))

            # Spline2D
            ydata = ws[:,i] # M
            Sfun = SmoothBivariateSpline(Vmeans, uoffss, ydata)
            c = splineContour(ax[0], Vmeans, uoffss, Sfun)#Ryfun, dx=1,dy=1)
            fig.colorbar(c, ax=ax[0])
            ax[0].set_xlabel('Vmean')
            ax[0].set_ylabel('uoffs')

            # Custom fit
            ffit = lambda xdata2 : wrenchMap(np.hstack((xdata2, np.zeros((xdata2.shape[0], 2)))), popts)[:,i]
            c = splineContour(ax[1], Vmeans, uoffss, ffit)
            fig.colorbar(c, ax=ax[1])
            ax[1].set_xlabel('Vmean')
            ax[1].set_ylabel('uoffs')

        # scatter vis
        fig = plt.figure()

        ax3d1 = fig.add_subplot(2,3,1, projection='3d')
        ax1 = [fig.add_subplot(2,3,2), fig.add_subplot(2,3,3)]
        fitWi(2, ax3d1, ax1)
        ax3d2 = fig.add_subplot(2,3,4, projection='3d')
        ax2 = [fig.add_subplot(2,3,5), fig.add_subplot(2,3,6)]
        fitWi(4, ax3d2, ax2)
        fig.tight_layout()
        plt.show()

    else:
        # save to file
        Vmeans = np.linspace(120, 160, num=10)
        uoffss = np.linspace(-0.5, 0.5, num=10)
        fs = [0.165]
        udiffs = np.linspace(-0.2, 0.2, num=10)
        h2s = [0.]
        sweepFile('test.npy', Vmeans, uoffss, fs, udiffs, h2s)
