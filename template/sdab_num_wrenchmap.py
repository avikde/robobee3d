import subprocess, sys, progressbar
import autograd.numpy as np
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
    return ax.contourf(xi, yi, zi, cmap='RdBu_r')
    
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
        # print(a2, self.A2)
        return a0, a1

    def f(self, xdata, *p):
        a0, a1 = self.unpackp(p)
        # print(a.shape, b, xdata.shape)
        y = a0 + xdata @ a1
        for i in range(xdata.shape[0]):
            xi = xdata[i,:]
            y[i] += 0.5 * xi.T @ self.A2 @ xi
        # print(y.shape)
        return y

    def fsingle(self, xi, *p):
        # To try and use autograd FIXME:
        a0, a1 = self.unpackp(p)
        # print(a.shape, b, xdata.shape)
        return a0 + np.dot(xi, a1) + 0.5 * xi.T @ self.A2 @ xi
    
    def df_dx(self, xi, *p):
        """This is dy/dx; NOT the jacobian wrt params required by the curve fit alg"""
        a0, a1 = self.unpackp(p)
        # y = a0 + a1 * x + x^T * A2 * x
        return a1 + self.A2 @ xi

fa = FunApprox(4) # k

"""Function for the numerical versions of the wrench map and its Jacobian.
Here u = [Vmean,uoffs,udiff,h2] is the input,
and w in R^6 is the output
"""

def wrenchMap(xdata, popts):
    """popts = (6,k)-shaped array of optimized params for each wrench component.
    xdata = N,Nu
    Returns N,6"""
    if len(xdata.shape) > 1:
        return np.vstack([fa.f(xdata, *popts[i,:]) for i in range(6)]).T
    else:
        # xdata = np.reshape(xdata,(1,len(xdata)))
        return np.hstack([fa.fsingle(xdata, *popts[i,:]) for i in range(6)])

def dw_du(xdata, popts):
    """xdata must be (Nu,) shaped"""
    return np.vstack([fa.df_dx(xdata, *popts[i,:]) for i in range(6)])

if __name__ == "__main__":
    np.set_printoptions(precision=2, suppress=True, linewidth=200)

    if len(sys.argv) > 1:
        with open(sys.argv[1], 'rb') as f:
            dat = np.load(f)
        Vmeans, uoffss, fs, udiffs, h2s, ws = unpackDat(dat)
        print("Unique in data:", np.unique(Vmeans), np.unique(uoffss), np.unique(fs), np.unique(udiffs), np.unique(h2s))
        
        xdata = np.vstack((Vmeans, uoffss, udiffs, h2s)).T # k,M
        xlabels = ['Vmean', 'uoffs', 'udiff', 'h2']

        # Optimized param fits in each row for each component of the wrench
        popts = np.vstack([curve_fit(fa.f, xdata, ws[:,i], p0=np.ones(fa.nparams()))[0] for i in range(6)])
        np.save('popts.npy', popts)

        def plotFitWi(ui1, ui2, wi, ax3d, ax):
            def lbl(ax):
                ax.set_xlabel(xlabels[ui1])
                ax.set_ylabel(xlabels[ui2])
            lbl(ax3d)
            for i in range(3):
                lbl(ax[i])
            def cplot(ax, ffit, ttl):
                c = splineContour(ax, xdata[:,ui1], xdata[:,ui2], ffit)
                fig.colorbar(c, ax=ax)
                ax.set_title(ttl)

            # scatter
            ax3d.plot(xdata[:,ui1], xdata[:,ui2], ws[:,wi], '.')
            ax3d.set_zlabel('W'+str(wi))

            # Spline2D
            ydata = ws[:,wi] # M
            Sfun = SmoothBivariateSpline(xdata[:,ui1], xdata[:,ui2], ydata)
            cplot(ax[0], Sfun, 'W'+str(wi)+' spline2D')

            # Custom fit
            ffit2 = lambda xdata2 : wrenchMap(np.hstack((xdata2, np.zeros((xdata2.shape[0], 2)))), popts)[:,wi]
            cplot(ax[1], ffit2, 'W'+str(wi)+' fit')
            # FIXME: the filling out of xdata only works for 0,1 ui

            # Jac d/dVmean
            ffit3 = lambda xdata2 : np.hstack([dw_du(np.hstack((xdata2[j,:], np.zeros(2))), popts)[wi,0] for j in range(xdata2.shape[0])])
            cplot(ax[2], ffit3, 'dW'+str(wi)+'/dVmean')

        # scatter vis
        fig = plt.figure()

        ax3d1 = fig.add_subplot(2,4,1, projection='3d')
        ax1 = [fig.add_subplot(2,4,2), fig.add_subplot(2,4,3), fig.add_subplot(2,4,4)]
        plotFitWi(0, 1, 2, ax3d1, ax1)
        ax3d2 = fig.add_subplot(2,4,5, projection='3d')
        ax2 = [fig.add_subplot(2,4,6), fig.add_subplot(2,4,7), fig.add_subplot(2,4,8)]
        plotFitWi(0, 2, 2, ax3d2, ax2)
        # plotFitWi(0, 1, 4, ax3d2, ax2)
        # fig.tight_layout()
        plt.show()

    else:
        # save to file
        Vmeans = np.linspace(120, 180, num=8)
        uoffss = np.linspace(-0.5, 0.5, num=8)
        fs = [0.17]
        udiffs = np.linspace(-0.2, 0.2, num=8)
        h2s = np.linspace(-0.2, 0.2, num=8)#[0.]
        sweepFile('test2.npy', Vmeans, uoffss, fs, udiffs, h2s)
