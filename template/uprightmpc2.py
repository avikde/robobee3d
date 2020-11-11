
import numpy as np
from scipy.spatial.transform import Rotation
np.set_printoptions(precision=4, suppress=True, linewidth=200)
from genqp import quadrotorNLDyn
from template_controllers import createMPC, reactiveController
from time import perf_counter
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.cm as cm
import progressbar

def viewControlTestLog(log, log2=None, callShow=True, goal0=False, desTraj=False, vscale=0.4):
    def traj3plot(_ax, t, p, v, cmap, narrow=20):
        cnorm = t/t[-1]
        # _ax.scatter(p[:,0], p[:,1], p[:,2], c=cnorm, cmap=cmap, marker='.', label='_nolegend_')
        
        # https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/multicolored_line.html
        points = np.array([p[:,0], p[:,1], p[:,2]]).T.reshape(-1, 1, 3)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        # Create the 3D-line collection object
        lc = Line3DCollection(segments, cmap=cmap, linewidths=2, label='_nolegend_')
        lc.set_array(cnorm) 
        lc.set_linewidth(2)
        _ax.add_collection3d(lc, zs=p[:,-1], zdir='z')

        # Upright arrows
        ii = np.linspace(0, len(t), narrow, dtype=int, endpoint=False)
        v *= vscale
        _ax.quiver(p[ii,0], p[ii,1], p[ii,2], v[ii,0], v[ii,1], v[ii,2], color='b' if "Blue" in cmap else 'r')

    def aspectEqual3(_ax, xyz):
        X, Y, Z = xyz[:,0], xyz[:,1], xyz[:,2]
        # Create cubic bounding box to simulate equal aspect ratio
        max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
        Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
        Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
        Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
        # Comment or uncomment following both lines to test the fake bounding box:
        for xb, yb, zb in zip(Xb, Yb, Zb):
            _ax.plot([xb], [yb], [zb], 'w', label='_nolegend_')

    def posParamPlot(_ax):
        traj3plot(_ax, log['t'], log['y'][:,:3], log['y'][:,3:6], "Blues_r")
        aspectEqual3(_ax, log['y'][:,:3])
        if log2 is not None:
            traj3plot(_ax, log2['t'], log2['y'][:,:3], log2['y'][:,3:6], "Reds_r")
        # _ax.plot(log['t'], log['pdes'][:,0], 'k--', alpha=0.3)
        if goal0:
            _ax.plot([0], [0], [0], 'g*', markersize=10, zorder=10)
            _ax.legend(('MPC', 'Reactive', 'Goal'))
        else:
            _ax.legend(('MPC', 'Reactive'))
        if desTraj:
            _ax.plot(log['pdes'][:,0], log['pdes'][:,1], log['pdes'][:,2], 'k--', alpha=0.5, zorder=9)
        _ax.set_xlabel('x [mm]')
        _ax.set_ylabel('y [mm]')
        _ax.set_zlabel('z [mm]')

    def posPlot(_ax):
        _ax.plot(log['t'], log['y'][:,:3])
        if log2 is not None:
            _ax.plot(log2['t'], log2['y'][:,:3], '--')
        _ax.plot(log['t'], log['pdes'][:,0], 'k--', alpha=0.3)
        _ax.set_ylabel('p')
    def splot(_ax):
        _ax.plot(log['t'], log['y'][:,3:6])
        _ax.axhline(y=0, color='k', alpha=0.3)
        _ax.set_ylabel('s')
    def inputsPlot(_ax1, _ax2):
        _ax1.plot(log['t'], log['u'][:,0])
        _ax1.axhline(y=0, color='k', alpha=0.3)
        _ax1.set_ylabel('Sp. thrust')
        _ax2.plot(log['t'], log['u'][:,1:])
        _ax2.axhline(y=0, color='k', alpha=0.3)
        _ax2.set_ylabel('Moments')
    def velsPlot(_ax1, _ax2):
        _ax1.plot(log['t'], log['y'][:,6:9])
        _ax1.axhline(y=0, color='k', alpha=0.3)
        _ax1.set_ylabel('v')
        if _ax2 is not None:
            _ax2.plot(log['t'], log['y'][:,9:12])
            _ax2.axhline(y=0, color='k', alpha=0.3)
            _ax2.set_ylabel('omega')
    def accdesPlots(_ax1, _ax2):
        _ax1.plot(log['t'], log['accdes'][:,:3])
        _ax1.axhline(y=0, color='k', alpha=0.3)
        _ax1.set_ylabel('accdes pos')
        _ax2.plot(log['t'], log['accdes'][:,3:])
        _ax2.axhline(y=0, color='k', alpha=0.3)
        _ax2.set_ylabel('accdes ang')
    def wlqpuPlots(_ax):
        _ax.plot(log['t'], log['wlqpu'][:,:2])
        _ax.plot(log['t'], log['wlqpu'][:,2:],'--')
        _ax.legend(('0','1','2','3'))
        _ax.set_ylabel('wlqpu')
    
    # fig = plt.figure()
    # ax = [fig.add_subplot(3,3,i+1) for i in range(1,12)]
    # posPlot(ax[0])
    # splot(ax[1])
    # inputsPlot(ax[2], ax[3])
    # velsPlot(ax[4], ax[5])
    # accdesPlots(ax[6], ax[7])
    # wlqpuPlots(ax[8])
    # fig.tight_layout()
    
    fig = plt.figure()
    ax3d = fig.add_subplot(1,1,1,projection='3d')
    posParamPlot(ax3d)

    if callShow:
        plt.show()

def controlTest(mdl, tend, dtsim=0.2, useMPC=True, trajFreq=0, trajAmp=0, ascentIC=False, showPlots=True, tpert=None, speedTest=False, perchTraj=False, flipTask=False, taulim=100, **kwargs):
    """trajFreq in Hz, trajAmp in mm"""
    if mdl is None and useMPC:
        mdl, _ = createMPC(**kwargs)
    speedTestvdes = 2 # m/s
    # Initial conditions
    dq = np.zeros(6)
    if ascentIC or speedTest or perchTraj:
        p = np.array([0, 0, -50]) if ascentIC else np.array([-500*speedTestvdes, 0, 0])
        if perchTraj:
            p = np.array([-100, 0, 0])
        Rb = np.eye(3)
    else:
        p = np.array([0, 0, 0])
        Rb = np.eye(3) if flipTask else Rotation.from_euler('xyz', [0.5,-0.5,0]).as_matrix()
        dq[0] = 0.1
    pdes = np.zeros(3)
    dpdes = np.zeros(3)
    sdes = np.array([0,0,1])
    
    tt = np.arange(tend, step=dtsim)
    Nt = len(tt)

    # for logging
    log = {'t': tt, 'y': np.zeros((Nt, 12)), 'u': np.zeros((Nt, 3)), 'pdes': np.zeros((Nt, 3)), 'accdes': np.zeros((Nt,6)), 'wlqpu': np.zeros((Nt,4))}

    trajOmg = 2 * np.pi * trajFreq * 1e-3 # to KHz, then to rad/ms
    ddqdes = None # test integrate ddq sim below

    avgTime = 0.0

    for ti in range(Nt):
        # Traj to follow
        if flipTask:
            # rotation phase 0 to 1
            ph = np.clip((tt[ti] - 100) / 200, 0, 1)
            sdes = np.array([-np.sin(ph*2*np.pi), 0, np.cos(ph*2*np.pi)])
        elif perchTraj:
            if tt[ti] < 500:
                pdes[0] = -100 + 0.2 * tt[ti]
                dpdes[0] = 0.2
                # rotation phase 0 to 1
                ph = np.clip((tt[ti] - 450) / 100, 0, 1)
                sdes = np.array([-np.sin(ph*np.pi), 0, np.cos(ph*np.pi)])
            else:
                pdes[0] = dpdes[0] = 0
                sdes = np.array([-1,0,0])
        elif speedTest:
            if tt[ti] < 500:
                dpdes[0] = speedTestvdes
                pdes[0] = -500*speedTestvdes + speedTestvdes*(tt[ti])
            else:
                pdes[0] = 0
                dpdes[0] = 0
        else:
            pdes[0] = trajAmp * np.sin(trajOmg * tt[ti])
            dpdes[0] = trajAmp * trajOmg * np.cos(trajOmg * tt[ti])
            if trajAmp > 1e-3:
                pdes[2] = 0.1 * tt[ti]
                dpdes[2] = 0.1
                # Add perturbation for this traj
                if tpert is not None and tt[ti] > tpert:
                    dq[1] += 2
                    tpert = None

        # Call controller
        if useMPC:
            t1 = perf_counter()
            # what "u" is depends on w(u). Here in python testing with w(u) = [0,0,u0,u1,u2,u3]
            w0 = np.hstack((0,0,mdl.u0))
            dwdu0 = np.vstack((np.zeros((2,4)), np.eye(4)))
            u, log['accdes'][ti,:], uwlqp = mdl.update(p, Rb, dq, pdes, dpdes, sdes)
            log['wlqpu'][ti,:] = uwlqp
            avgTime += 0.01 * (perf_counter() - t1 - avgTime)
            # # Alternate simulation by integrating accDes
            # ddqdes = accdess[ti,:]
        else:
            u = reactiveController(p, Rb, dq, pdes, **kwargs)
        # u = np.array([1,0.1,0])
        # Input limit
        for i in range(2):
            u[i+1] = np.clip(u[i+1], -taulim, taulim)

        p, Rb, dq = quadrotorNLDyn(p, Rb, dq, u, dtsim, ddq=ddqdes)
        log['y'][ti,:] = np.hstack((p, Rb[:,2], dq))
        log['u'][ti,:] = u
        log['pdes'][ti,:] = pdes
    if useMPC and showPlots:
        print("Time (ms):", avgTime * 1e3)
    if showPlots:
        viewControlTestLog(log)
    return log

def logMetric(log):
    # A metric to plot about how good the tracking was
    Nt = len(log['t'])
    perr = log['y'][:,:3]
    tau = log['u'][:,1:3]
    serr = log['y'][:,3:6]
    serr[:,2] -= 1.0
    err = 0
    eff = 0
    for i in range(Nt):
        err += np.dot(perr[i,:], perr[i,:])# + 10 * np.dot(serr[i,:], serr[i,:])
        eff += np.dot(tau[i,:], tau[i,:])# + 10 * np.dot(serr[i,:], serr[i,:])
    err /= Nt
    eff /= Nt
    return err, eff

def papPlots(bmpc):
    """Baseline mpc as argument"""
    def flipTask():
        l1 = controlTest(bmpc, 1000, useMPC=True, showPlots=False, flipTask=True)
        viewControlTestLog(l1, desTraj=True, vscale=10)
        fig, ax = plt.subplots(1,3, figsize=(7.5,2.5))
        for i in range(0,3,2):
            ax[i].plot(1e-3*l1['t'], l1['y'][:,i], 'b')
            ax[i].plot(1e-3*l1['t'], l1['pdes'][:,i], 'k--', alpha=0.3)
            ax[i].set_xlabel('t [s]')
            
        ax[1].plot(1e-3*l1['t'], 180/np.pi*np.arctan2(l1['y'][:,3], l1['y'][:,5]), 'b')
        ax[1].plot([0, 0.1, 0.2], [0, 0, -180], 'k--', alpha=0.3)
        ax[1].plot([0.2, 0.3, 1], [180, 0, 0], 'k--', alpha=0.3)
        ax[1].set_ylabel('Angle [deg]')
        ax[0].set_ylabel('x [mm]')
        ax[2].set_ylabel('z [mm]')
        fig.tight_layout()
        plt.show()

    def perchTask():
        l1 = controlTest(bmpc, 550, useMPC=True, showPlots=False, perchTraj=True)
        viewControlTestLog(l1, desTraj=True, vscale=10)
        fig, ax = plt.subplots(1,3, figsize=(7.5,2.5))
        for i in range(0,3,2):
            ax[i].plot(1e-3*l1['t'], l1['y'][:,i], 'b')
            ax[i].plot(1e-3*l1['t'], l1['pdes'][:,i], 'k--', alpha=0.3)
            ax[i].set_xlabel('t [s]')
            
        ax[1].plot(1e-3*l1['t'], 180/np.pi*np.arctan2(l1['y'][:,3], l1['y'][:,5]), 'b')
        ax[1].plot([0, 0.45, 0.55], [0, 0, -90], 'k--', alpha=0.3)
        ax[1].set_ylabel('Angle [deg]')
        ax[0].set_ylabel('x [mm]')
        ax[2].set_ylabel('z [mm]')
        fig.tight_layout()
        plt.show()

    def hoverTask(show3d, reactiveArgs1, reactiveArgs2=None):
        l1 = controlTest(bmpc, 300, useMPC=True, showPlots=False)
        l2 = controlTest(bmpc, 1000, useMPC=False, showPlots=False, **reactiveArgs1)
        if reactiveArgs2 is not None:
            l3 = controlTest(bmpc, 1000, useMPC=False, showPlots=False, **reactiveArgs2)
        if show3d:
            viewControlTestLog(l1, log2=l2, goal0=True)
        else:
            fig, ax = plt.subplots(1,2, figsize=(5,2.5))
            for i in range(2):
                ax[i].plot(1e-3*l1['t'], l1['y'][:,i], 'b')
                ax[i].plot(1e-3*l2['t'], l2['y'][:,i], 'r')
                if reactiveArgs2 is not None:
                    ax[i].plot(1e-3*l3['t'], l3['y'][:,i], 'r--')
                ax[i].plot(1e-3*l2['t'], l2['pdes'][:,i], 'k--', alpha=0.3)
                ax[i].set_xlabel('t [s]')
            ax[0].set_ylabel('x [mm]')
            ax[1].set_ylabel('y [mm]')
            fig.tight_layout()
            plt.show()

    def sTask(show3d, reactiveArgs):
        l1 = controlTest(bmpc, 2000, useMPC=True, showPlots=False, trajAmp=50, trajFreq=1, tpert=1000)
        l2 = controlTest(bmpc, 2000, useMPC=False, showPlots=False, trajAmp=50, trajFreq=1, tpert=1000, **reactiveArgs)
        if show3d:
            viewControlTestLog(l1, log2=l2, desTraj=True, vscale=20)
        else:
            fig, ax = plt.subplots(1,2, figsize=(5,2.5))
            for i in range(2):
                ax[i].plot(1e-3*l1['t'], 1e-3*l1['y'][:,i], 'b')
                ax[i].plot(1e-3*l2['t'], 1e-3*l2['y'][:,i], 'r')
                ax[i].plot(1e-3*l2['t'], 1e-3*l2['pdes'][:,i], 'k--', alpha=0.3)
                ax[i].set_xlabel('t [s]')
            ax[0].set_ylabel('x [m]')
            ax[1].set_ylabel('y [m]')
            fig.tight_layout()
            plt.show()

    def accTask(reactiveArgs):
        l1 = controlTest(bmpc, 1000, useMPC=True, showPlots=False, speedTest=True)
        l2 = controlTest(bmpc, 1000, useMPC=False, showPlots=False, speedTest=True, **reactiveArgs)
        viewControlTestLog(l1, log2=l2, desTraj=True, vscale=100)
        fig, ax = plt.subplots(1,2, figsize=(5,2.5))
        ax[0].plot(1e-3*l1['t'], 1e-3*l1['y'][:,0], 'b')
        ax[0].plot(1e-3*l2['t'], 1e-3*l2['y'][:,0], 'r')
        ax[0].plot(1e-3*l2['t'], 1e-3*l2['pdes'][:,0], 'k--', alpha=0.3)
        ax[0].set_ylabel('x [m]')
        ax[0].set_xlabel('t [s]')
        ax[1].plot(1e-3*l1['t'], l1['y'][:,6], 'b')
        ax[1].plot(1e-3*l2['t'], l2['y'][:,6], 'r')
        ax[1].set_ylabel('xdot [m/s]')
        ax[1].set_xlabel('t [s]')
        fig.tight_layout()
        plt.show()

    # Hover tuning ---------
    def gainTuningSims(useMPC, kwgain, k1range, k2range, kwfixedn, kwfixedv, npts=10):
        k1s = np.linspace(*k1range,num=npts)
        k2s = np.linspace(*k2range,num=npts)
        xv, yv = np.meshgrid(k1s, k2s, indexing='ij') # treat xv[i,j], yv[i,j]
        costs = np.zeros_like(xv)
        efforts = np.zeros_like(xv)
        # create a progress bar
        widgets = [
            'Progress: ', progressbar.Percentage(),
            ' ', progressbar.Bar(),
            ' ', progressbar.ETA(),
        ]
        bar = progressbar.ProgressBar(widgets=widgets, max_value=np.prod(costs.shape))
        nrun = 0
        for i in range(len(k1s)):
            for j in range(len(k2s)):
                nrun += 1
                bar.update(nrun)
                try:
                    if useMPC:
                        kwgain = 'mpc_wpos' # overwrite filename
                        # Assume same ratio of final to running cost
                        kwargs = {'wpr': xv[i,j], 'wvr': yv[i,j], 'wpf': 5*xv[i,j], 'wvf': 2*yv[i,j]}
                    else: # Reactive
                        kwargs = {kwgain: [xv[i,j],yv[i,j]], kwfixedn: kwfixedv}
                    l2 = controlTest(None, 1000, useMPC=useMPC, showPlots=False, taulim=10, **kwargs)
                    costs[i,j], efforts[i,j] = logMetric(l2)
                except KeyboardInterrupt:
                    raise
                except:
                    costs[i,j] = efforts[i,j] = np.nan
        np.savez(kwgain+str('.npz'), xv=xv, yv=yv, costs=costs, efforts=efforts)

    def gainTuningPlots(maxcost=10):
        lmpc = controlTest(bmpc, 1000, useMPC=True, showPlots=False)
        empc, effmpc = logMetric(lmpc)
                        
        def plot1(ax, dat, cbar=False):
            costs = np.clip(dat['costs'] / empc, 0, maxcost)
            im = ax.pcolormesh(dat['xv'], dat['yv'], costs, cmap='gray_r', shading='auto', vmin=0, vmax=maxcost)
            if cbar:
                fig.colorbar(im)
        
        fig, ax = plt.subplots(1,3,figsize=(12,4))
        plot1(ax[0], np.load('ks.npz'))
        ax[0].plot([15], [100], 'r*', ms=20)
        plot1(ax[1], np.load('kpos.npz'))
        ax[1].plot([0.01, 0.04], [1.0, 1.25], 'r*', ms=20)
        plot1(ax[2], np.load('mpc_wpos.npz'), cbar=True)
        ax[2].plot([1], [1e3], 'b*', ms=20)
        plt.show()

    def trackingEffortPlot(files):
        # Baseline
        lmpc = controlTest(bmpc, 1000, useMPC=True, showPlots=False)
        empc, effmpc = logMetric(lmpc)
        costs2 = []
        effs2 = []
        costs2mpc = []
        effs2mpc = []
        for fname in files:
            dat = np.load(fname)
            costs = dat['costs'].ravel() / empc
            effs = dat['efforts'].ravel() / effmpc
            ii = np.where(costs < 10)[0]
            if 'mpc' in fname:
                costs2mpc.append(costs[ii])
                effs2mpc.append(effs[ii])
            else:
                costs2.append(costs[ii])
                effs2.append(effs[ii])
        fig, ax = plt.subplots(1, figsize=(4,4))
        if len(costs2mpc) > 0:
            ax.scatter(costs2mpc, effs2mpc, color='b', label='MPC')
        ax.scatter(costs2, effs2, color='r', label='Reactive')
        ax.axhline(1, color='k', linestyle='dashed', alpha=0.3)
        ax.axvline(1, color='k', linestyle='dashed', alpha=0.3)
        ax.set_xlim((0,10))
        ax.set_ylim((0,10))
        ax.set_aspect('equal')
        ax.set_xlabel('Relative tracking error [ ]')
        ax.set_ylabel('Relative actuator effort [ ]')
        ax.legend()
        plt.show()

    # Tuning ---------------
    # # Run and save data
    # # defaults kpos=[5e-3,5e-1], kz=[1e-1,1e0], ks=[10e0,1e2]
    # gainTuningSims(False, 'ks', [5e0,2e1], [2e1,2e2], 'kpos', [5e-3,5e-1])
    # gainTuningSims(False, 'kpos', [1e-3,8e-2], [1e-1,2e0], 'ks', [15,100])
    # # defaults wpr=1, wvr=1e3, wpf=5, wvf=2e3
    # gainTuningSims(True, 'wpos', [0.5,10], [0.5e3, 10e3], None, None)
    
    # hoverTask(False, {'ks':[15,100], 'kpos':[0.01,1]}, {'ks':[15,100], 'kpos':[0.04,1.25]})
    # sTask(False, ks=[15,100], kpos=[0.01,1])

    # gainTuningPlots()
    # trackingEffortPlot(['mpc_wpos.npz','kpos.npz'])

    # sim1hover -------------------
    # hoverTask(False, {'ks':[15,100], 'kpos':[0.01,1]})
    # sTask(False, {'ks':[15,100], 'kpos':[0.01,1]})
    accTask({'ks':[15,100], 'kpos':[0.01,1]})
    # sim1perch --------------
    # flipTask()
    # perchTask()

if __name__ == "__main__":
    up, upc = createMPC()

    # # Hover
    # controlTest(up, 500, useMPC=True)
    # # Ascent
    # controlTest(up, 500, useMPC=True, ascentIC=True)
    # # Traj
    # controlTest(up, 2000, useMPC=True, trajAmp=50, trajFreq=1)

    papPlots(up)
