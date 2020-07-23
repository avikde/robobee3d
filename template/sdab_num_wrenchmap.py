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
        # qw below contains wing kinematics as well as the wrench
        sw = bee.openLoop(Vmean * (1 + udiff), Vmean * (1 - udiff), uoffs, f, h2=h2, verbose=False)
        NptsInPeriod = int(1/(f * bee.TIMESTEP))
        # avg wrench
        avgwrench = np.mean(sw[-NptsInPeriod:,-6:], axis=0)
        # amplitudes
        qw = sw[-NptsInPeriod:,:4]
        dqw = sw[-NptsInPeriod:,4:8]
        # to estimate alpha (fraction of period for upstroke)
        ratioPosStrokeVel = np.sum(dqw[:,0] >= 0, axis=0) / NptsInPeriod
        # for each wing, get max and min amplitude (corresponding to upstroke and downstroke)
        kins = np.hstack([np.hstack((np.amax(qw[:,2*i:2*i+2], axis=0), -np.amin(qw[:,2*i:2*i+2], axis=0))) for i in range(2)]) # size 8
        kins = np.hstack((kins, ratioPosStrokeVel)) # size 9
        return np.hstack((avgwrench, kins))

    res = np.array([
        np.hstack((Vmean, uoffs, f, udiff, h2, 
        olAvgWrench(Vmean, uoffs, f, udiff, h2))) 
        for Vmean in Vmeans for uoffs in uoffss for f in fs for udiff in udiffs for h2 in h2s])
    with open(fname, 'wb') as f:
        np.save(f, res)

    # xdata = np.array([
    #     np.hstack((Vmean, uoffs, f, udiff, h2)) 
    #     for Vmean in Vmeans for uoffs in uoffss for f in fs for udiff in udiffs for h2 in h2s])
    # cpus = multiprocessing.cpu_count()
    # pool = multiprocessing.Pool(processes=cpus)
    # res = pool.map(test, xdata)
    # TODO: multiprocessing but need different copies of the sim too

unpackDat = lambda dat : (dat[:,0], dat[:,1], dat[:,2], dat[:,3], dat[:,4], dat[:,5:11], dat[:,11:])

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

def wrenchFromKinematics(kins, freq, params):
    """Analytical prediction of average wrench from kinematics features. See w2d_template.nb."""
    if len(kins.shape) > 1:
        # apply to each row and return result
        N = kins.shape[0]
        return np.array([wrenchFromKinematics(kins[i,:], freq[i], params) for i in range(N)])
    
    # unpack params
    CD0 = robobee.CD0
    CDmax = robobee.CDmax
    CLmax = robobee.CLmax
    f2 = freq**2
    kaero = 1/2 * robobee.RHO * params['Aw']**2 * params['AR'] * params['r2h']**2
    ycp = params['ycp']

    def iwrencha(Phim, Phid, Psi1, Psi2, alpha):
        """From w2d_template.nb"""
        c2Psi1, s2Psi1 = np.cos(2*Psi1), np.sin(2*Psi1)
        c2Psi2, s2Psi2 = np.cos(2*Psi2), np.sin(2*Psi2)
        cPhid, sPhid = np.cos(Phid), np.sin(Phid)
        cPhim, sPhim = np.cos(Phim), np.sin(Phim)
        Phid2 = Phid**2
        return np.array([
            (-4*((-2 + alpha)*c2Psi1*(CD0 - CDmax) + alpha*c2Psi2*(CD0 - CDmax) - 2*(-1 + alpha)*(CD0 + CDmax))*cPhim*f2*kaero*Phid*sPhid)/((-2 + alpha)*alpha),
            (-4*((-2 + alpha)*c2Psi1*(CD0 - CDmax) + alpha*c2Psi2*(CD0 - CDmax) - 2*(-1 + alpha)*(CD0 + CDmax))*f2*kaero*Phid*sPhid*sPhim)/((-2 + alpha)*alpha),
            (8*CLmax*f2*kaero*Phid2*((-2 + alpha)*s2Psi1 + alpha*s2Psi2))/((-2 + alpha)*alpha),
            (8*CLmax*cPhim*f2*kaero*Phid*((-2 + alpha)*s2Psi1 + alpha*s2Psi2)*sPhid*ycp)/((-2 + alpha)*alpha),
            (8*CLmax*f2*kaero*Phid*((-2 + alpha)*s2Psi1 + alpha*s2Psi2)*sPhid*sPhim*ycp)/((-2 + alpha)*alpha),
            (4*((-2 + alpha)*c2Psi1*(CD0 - CDmax) + alpha*c2Psi2*(CD0 - CDmax) - 2*(-1 + alpha)*(CD0 + CDmax))*f2*kaero*Phid2*ycp)/((-2 + alpha)*alpha)
            ])
    
    # unpack the stored numerical kinematics features
    ampls = kins[:4], kins[4:8]
    dalpha = kins[8] - 1.0

    # Symmetry mapping right to left
    Symmw = np.array([1,-1,1,-1,1,-1])

    return iwrencha(*ampls[0], 1 + dalpha) + Symmw * iwrencha(*ampls[1], 1 - dalpha)

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
        Vmeans, uoffss, fs, udiffs, h2s, ws0, kins = unpackDat(dat)
        params = robobee.wparams.copy()
        params.update({'ycp': 10, 'AR': 5})
        ws = wrenchFromKinematics(kins, fs, params)
        # TODO: compare ws0 to ws
        print(ws.shape)
        sys.exit()

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
        Vmeans = np.linspace(90, 160, num=8)
        uoffss = np.linspace(-0.5, 0.5, num=8)
        fs = [0.17]
        udiffs = np.linspace(-0.2, 0.2, num=8)
        h2s = np.linspace(-0.2, 0.2, num=8)#[0.]
        sweepFile('numkins.npy', Vmeans, uoffss, fs, udiffs, h2s)
