import sys, os, itertools
import autograd.numpy as np
from scipy.interpolate import SmoothBivariateSpline
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import robobee

# Load empirical data and convert to the same format as tested with sim ---------------------------------

def loadEmpiricalData(fnameCSV):
    # Loading the sim-generated data looks like this `Vmeans, uoffss, fs, udiffs, h2s, ws0, kins = unpackDat(dat)``
    # From the empirical data should be able to replicate all this except for ws0
    datcsv = np.genfromtxt(fnameCSV, delimiter=",", skip_header=2)
    # unpack https://docs.google.com/spreadsheets/d/1Sa6lT008fpYqdgcjl0M7Nx5MlD-OYx808FP6nt6We5A/edit#gid=0
    # fs, Vleft, Vright, drv_pch, h2 = datcsv[:,0], datcsv[:,1], datcsv[:,2], datcsv[:,3], datcsv[:,4]
    rawInp, rawT, rawB = datcsv[:,:5], datcsv[:,5:9], datcsv[:,9:13]
    Nrows = rawInp.shape[0]
    # print(datcsv)

    def convertVmean(Inp, T, B):
        """Convert Vleft/Vright from empirical trials to different combinations of Vmean/udiff"""
        outU = []
        outK = []
        # Produce rows from combinations https://github.com/avikde/robobee3d/pull/172#issuecomment-669311864
        Vleft, Vright = Inp[:,1], Inp[:,2]
        uVolts = np.unique(Vleft)
        combs = list(itertools.combinations(uVolts,2))
        inds = list(itertools.combinations(range(len(uVolts)), 2))
        # print('Vleft,Vright combinations', combs)

        def addRow(Vmean, udiff, kinT, kinB):
            uoffs = Inp[0,3]/Vmean # should be the same
            h2 =  Inp[0,4]
            alpha = 0.5 # TODO:
            outU.append(np.array([Vmean, uoffs, udiff, h2]))
            outK.append(np.hstack((kinT, kinB, 0.5)))

        for i in range(len(combs)):
            i0, i1 = inds[i] # i0,i1 are the indices in the original collected data
            # Get 4 different combinations
            # NOTE: assuming vleft = vright in each original data row
            V1, V2 = combs[i]
            Vm = 0.5*(V1+V2)
            addRow(V1, 0, T[i0,:], B[i0,:])
            # invert the CC [1+udiff;1-udiff]*Vm = [V1;V2]
            udiff = V2/Vm - 1# = 1-V1/Vm
            addRow(Vm, udiff, T[i0,:], B[i1,:])
            addRow(Vm, -udiff, T[i1,:], B[i0,:])
            addRow(V2, 0, T[i1,:], B[i1,:])
        # convert to 2D array
        outU = np.array(outU)
        outK = np.array(outK)
        return outU, outK

    def convertRawKins(rawKins, topSign):
        """Convert the empirically measured kinematics features to the ones used in the sim; working on one wing at a time
        *** In sim:
        For each wing, get max and min amplitude (corresponding to upstroke and downstroke)
        np.hstack((np.amax(qw[:,2*i:2*i+2], axis=0), -np.amin(qw[:,2*i:2*i+2], axis=0)))
        = [num max stroke, num max pitch, -num min stroke, -num min pitch] -- all should be positive numbers
        Sim data plot https://github.com/avikde/robobee3d/pull/172#issuecomment-671369920.
        
        *** Empirical:
        All positive numbers, have [stroke L, pitch L->R, stroke R, pitch R->L]. Would help to see this data for a sim trial for comparison.
        Comparing with the sim plot on github, stroke L, stroke R seem to match in both, except for deg to rad
        """
        scaledKins = np.zeros_like(rawKins)
        # NOTE: assuming top = left wing otherwise these need to be flipped
        scaledKins[:,[0,2]] = np.radians(rawKins[:,[0,2]])
        if topSign > 0:
            # ~ sin(t)
            scaledKins[:,1] = np.radians(rawKins[:,1])
            scaledKins[:,3] = np.radians(rawKins[:,3])
        else:
            scaledKins[:,1] = np.radians(rawKins[:,3])
            scaledKins[:,3] = np.radians(rawKins[:,1])
        return scaledKins

    # First find all the rows with the same uoffs, h2 (except Vmean, udiff)
    drv_pch = rawInp[:,3]
    uniqueuoffs, invinds = np.unique(drv_pch, return_inverse=True)
    outU = []
    outK = []
    # Next convert (Vleft/Vright) to Vmean/udiff from all co
    for i in range(len(uniqueuoffs)):
        originds = np.where(invinds == i)[0]
        outUi, outKi = convertVmean(rawInp[originds,:], convertRawKins(rawT[originds,:], 1), convertRawKins(rawB[originds,:], -1))
        outU.append(outUi)
        outK.append(outKi)
        # outU
    outU = np.vstack(outU)
    outK = np.vstack(outK)
    # print(uniqueuoffs, uinds)
    Nrows2 = outU.shape[0]
    fs = np.ones((Nrows2, 1)) * rawInp[0,0] * 1e-3 # to KHz
    # print(outU.shape, outK.shape)
    # outU = Vmean, uoffs, udiff, h2

    dat = np.hstack((outU[:,:2], fs, outU[:,2:], np.zeros((Nrows2, 6)), outK))
    # print(dat.shape)
    return dat

# Load data and fit a function --------------------------------------------------------------------------

# Vmeans, uoffss, fs, udiffs, h2s, ws0, kins = unpackDat(dat)
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

def wrenchFromKinematics(kins, freq, params, kaerox=1, strokex=1, wbias=np.zeros(6)):
    """Analytical prediction of average wrench from kinematics features. See w2d_template.nb."""
    if len(kins.shape) > 1:
        # apply to each row and return result
        N = kins.shape[0]
        return np.array([wrenchFromKinematics(kins[i,:], freq[i], params, kaerox=kaerox, strokex=strokex, wbias=wbias) for i in range(N)])
    
    # unpack params
    CD0 = robobee.CD0
    CDmax = robobee.CDmax
    CLmax = robobee.CLmax
    f2 = freq**2
    kaero = 1/2 * robobee.RHO * params['Aw']**2 * params['AR'] * params['r2h']**2 * kaerox
    ycp = params['ycp']
    R = params['R']

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
            (8*CLmax*f2*kaero*Phid*((-2 + alpha)*s2Psi1 + alpha*s2Psi2)*(Phid*R +cPhim*sPhid*ycp))/((-2 + alpha)*alpha),
            (8*CLmax*f2*kaero*Phid*((-2 + alpha)*s2Psi1 + alpha*s2Psi2)*sPhid*sPhim*ycp)/((-2 + alpha)*alpha),
            (4*((-2 + alpha)*c2Psi1*(CD0 - CDmax) + alpha*c2Psi2*(CD0 - CDmax) - 2*(-1 + alpha)*(CD0 + CDmax))*f2*kaero*Phid*(cPhim*R*sPhid + Phid*ycp))/((-2 + alpha)*alpha)
            ])
    
    def remapKins(q1max, q2max, q1nmin, q2nmin):
        Phim = 0.5 * (q1max - q1nmin)
        Phid = 0.5 * (q1max + q1nmin) * strokex
        Psi1 = q2max
        Psi2 = -q2nmin
        return Phim, Phid, Psi1, Psi2

    # unpack the stored numerical kinematics features
    ampls = remapKins(*kins[:4]), remapKins(*kins[4:8])
    dalpha = 2*kins[8] - 1.0

    # Symmetry mapping right to left
    Symmw = np.array([1,-1,1,-1,1,-1])

    return iwrencha(*ampls[0], 1 - dalpha) + Symmw * iwrencha(*ampls[1], 1 + dalpha) + np.asarray(wbias)

def wrenchCompare(ws, ws2):
    fig, ax = plt.subplots(6)
    for i in range(6):
        ax[i].plot(ws[:,i], '.')
        ax[i].plot(ws2[:,i], '.')
    plt.show()

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
    np.set_printoptions(precision=2, suppress=False, linewidth=100000)

    ext = os.path.splitext(sys.argv[1])[1]
    if ext == '.npy':
        with open(sys.argv[1], 'rb') as f:
            dat = np.load(f)
    elif ext == '.csv':
        dat = loadEmpiricalData(sys.argv[1])
    Vmeans, uoffss, fs, udiffs, h2s, ws0, kins = unpackDat(dat)

    params = robobee.wparams.copy()
    params.update({'ycp': 7.5, 'AR': 4.5, 'R': 3})
    # NOTE: check this bias
    ws = wrenchFromKinematics(kins, fs, params, kaerox=1, strokex=1.5)#, wbias=[0,0,0,0,-3,0])
    
    # wrenchCompare(ws0, ws) # compare ws0 to ws
    # sys.exit()

    print("Unique in data:", np.unique(Vmeans), np.unique(uoffss), np.unique(fs), np.unique(udiffs), np.unique(h2s))
    
    xdata = np.vstack((Vmeans, uoffss, udiffs, h2s)).T # k,M
    xlabels = ['Vmean', 'uoffs', 'udiff', 'h2']

    # Optimized param fits in each row for each component of the wrench
    popts = np.vstack([curve_fit(fa.f, xdata, ws[:,i], p0=np.ones(fa.nparams()))[0] for i in range(6)])
    print('popts row major =')
    print(np.array2string(np.ravel(popts, order='C'), separator=','))
    # np.save('poptsEmp2.npy', popts)

    def plotFitWi(ui1, ui2, wi, ax3d, ax):
        def lbl(ax):
            ax.set_xlabel(xlabels[ui1])
            ax.set_ylabel(xlabels[ui2])
        lbl(ax3d)
        for i in range(2):
            lbl(ax[i])
        def cplot(ax, ffit, ttl):
            c = splineContour(ax, xdata[:,ui1], xdata[:,ui2], ffit)
            fig.colorbar(c, ax=ax)
            ax.set_title(ttl)

        def fitSurface(ax, xiu, yiu, Zfun, length=50):
            dimrow = lambda row : np.linspace(row.min(), row.max(), length)
            xi = dimrow(xiu)
            yi = dimrow(yiu)
            
            # create a grid for the plot
            Xi, Yi = np.meshgrid(xi, yi)
            positions = np.vstack([Xi.ravel(), Yi.ravel()])
            zi = np.reshape(Zfun(positions.T), Xi.shape)
            Zi = zi.reshape(Xi.shape)
            # zi = Zfun(Xi, Yi, np.zeros_like(Xi))
            return ax.plot_surface(Xi, Yi, Zi, cmap='RdBu_r')
        
        if ui1 == 0 and ui2 in [1,2]:
            pass
        else:
            print("WARNING!!! the surface plots are picking the wrong slice for ui =", ui1, ui2)

        # scatter
        ax3d.plot(xdata[:,ui1], xdata[:,ui2], ws[:,wi], '.')
        ax3d.set_zlabel('W'+str(wi))

        # # Spline2D
        # ydata = ws[:,wi] # M
        # Sfun = SmoothBivariateSpline(xdata[:,ui1], xdata[:,ui2], ydata)
        # cplot(ax[0], Sfun, 'W'+str(wi)+' spline2D')

        # Custom fit
        def ffit2(xdata2):
            testxdataSection = np.zeros((xdata2.shape[0], 4))
            testxdataSection[:,ui1] = xdata2[:,0]
            testxdataSection[:,ui2] = xdata2[:,1]
            return wrenchMap(testxdataSection, popts)[:,wi]
        cplot(ax[0], ffit2, 'W'+str(wi)+' fit')
        fitSurface(ax3d, xdata[:,ui1], xdata[:,ui2], ffit2)
        # FIXME: the filling out of xdata only works for 0,1 ui

        # Jac d/dVmean
        ffit3 = lambda xdata2 : np.hstack([dw_du(np.hstack((xdata2[j,:], np.zeros(2))), popts)[wi,0] for j in range(xdata2.shape[0])])
        cplot(ax[1], ffit3, 'dW'+str(wi)+'/dVmean')

    # scatter vis
    fig = plt.figure()

    ax3d1 = fig.add_subplot(2,3,1, projection='3d')
    ax1 = [fig.add_subplot(2,3,2), fig.add_subplot(2,3,3)]
    plotFitWi(0, 1, 2, ax3d1, ax1)
    ax3d2 = fig.add_subplot(2,3,4, projection='3d')
    ax2 = [fig.add_subplot(2,3,5), fig.add_subplot(2,3,6)]
    # plotFitWi(0, 2, 2, ax3d2, ax2)
    # plotFitWi(0, 1, 4, ax3d2, ax2) # Vmean, uoffs -> pitch
    plotFitWi(0, 2, 3, ax3d2, ax2) # Vmean, udiff -> roll
    # fig.tight_layout()
    plt.show()
