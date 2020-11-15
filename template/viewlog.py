
import numpy as np
import gzip, pickle, glob, os, time
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from ca6dynamics import M
from plot_helpers import *

def initLog():
    return {'t': [], 'q': [], 'dq': [], 'u': [], 'accdes': [], 'posdes': []}

def appendLog(data, t, q, dq, u, accdes, posdes):
    """rot can be any scipy rotation type"""
    lastT = -np.inf if len(data['t']) == 0 else data['t'][-1]
    if t - lastT >= 0.999:
        data['t'].append(t)
        data['q'].append(np.copy(q))
        data['dq'].append(np.copy(dq))
        data['u'].append(np.copy(u))
        data['accdes'].append(np.copy(accdes))
        data['posdes'].append(np.copy(posdes))
    return data

def saveLog(f1, data):
    for k in data.keys():
        data[k] = np.array(data[k])

    t = time.localtime()
    timestamp = time.strftime('%Y%m%d%H%M%S', t)
    fname = f1 + '_' + timestamp + '.zip'
    zfile = gzip.GzipFile(fname, 'wb')
    pickle.dump(data, zfile)
    zfile.close()
    print('Saved log in', fname)

def addKeys(data):
    # Convert and add more keys for easier processing
    qb = data['q'][:,-7:]
    dqb = data['dq'][:,-6:]
    data['s'] = np.zeros((len(data['t']), 3))
    for i in range(len(data['t'])):
        data['s'][i,:] = Rotation.from_quat(qb[i,3:]).as_matrix()[:,2]
    data['eul'] = Rotation.from_quat(qb[:,3:]).as_euler('xyz')
    # Split up q
    data['p'] = qb[:,:3]
    data['dp'] = dqb[:,:3]
    data['omega'] = dqb[:,3:]
    return data

def readFile(fname, dataFmtList=None):
    f1, ext = os.path.splitext(fname)
    zfile = gzip.GzipFile(fname, 'rb')
    data = pickle.load(zfile)
    # print(data)
    print('Opened '+fname+'; average data rate = ' + str(1000.0/np.mean(np.diff(data['t'])))+'Hz')
    return addKeys(data)

def getData(fname):
    #get file name
    if not fname:
        # open the latest file
        files = glob.glob('../logs/' + '*')
        files.sort(key=os.path.getmtime, reverse=True)
        fname = files[0]
    return readFile(fname), 'ca6' in fname

def defaultPlots(data, ca6log=False):
    t = data['t']
    # SDAB log also has wings

    if 'p' not in data.keys():
        data = addKeys(data)

    fig, ax = plt.subplots(4,2)
    ax = ax.ravel()
    ax[0].plot(data['t'], data['p'])
    ax[0].plot(data['t'], data['posdes'][:,0], 'b--')
    ax[0].plot(data['t'], data['posdes'][:,1], 'r--')
    ax[0].plot(data['t'], data['posdes'][:,2], 'g--')
    ax[0].set_ylabel('pos [mm]')

    ax[1].plot(data['t'], data['s'])
    ax[1].set_ylabel('s')
    ax[2].plot(data['t'], data['eul'])
    ax[2].set_ylabel('orn [rad]')

    ax[3].plot(data['t'], data['accdes'][:,:3])
    ax[3].set_ylabel('Accdes pos')
    ax[4].plot(data['t'], data['accdes'][:,3:])
    ax[4].set_ylabel('Accdes ang')
    # actMom = (M @ dqb.T).T
    # ax[2].plot(data['t'], actMom[:,2], label='act')
    # ax[2].legend()

    if ca6log:
        ax[5].plot(data['t'], data['u'][:,[0,3]])
        ax[5].set_ylabel('u1')
        # ca6 log
        ax[6].plot(data['t'], data['u'][:,[1,4]])
        ax[6].set_ylabel('u2')
        ax[7].plot(data['t'], data['u'][:,[2,5]])
        ax[7].set_ylabel('u3')
    else:
        ax[5].plot(data['t'], data['omega'])
        ax[5].set_ylabel('Omega')
        # Inputs
        ax[6].plot(t, data['u'][:,2]) # Vmean
        ax[6].set_ylabel('Vmean')
        ax[7].plot(t, data['u'][:,3], label='offs')
        ax[7].plot(t, data['u'][:,4], label='diff')
        ax[7].plot(t, data['u'][:,5], label='h2')
        ax[7].set_ylabel('u')
        ax[7].legend()
        # # plot wing states
        # qw = data['q'][:,:4]
        # dqw = data['dq'][:,:4]
        # ax[5].plot(t, qw[:,[0,2]])
        # ax[5].set_ylabel('Stroke')
        # ax[6].plot(t, qw[:,[1,3]])
        # ax[6].set_ylabel('Pitch')

    ax[-1].set_xlabel('Time [ms]')
    fig.tight_layout()

def papPlots(l1, l2, vscale=50, traj3d=False, plotset='helix', purple=False):
    if isinstance(l1, str):
        l1 = readFile(l1)
        l2 = readFile(l2)

    if traj3d:
        fig = plt.figure()
        ax3d = fig.add_subplot(1,1,1,projection='3d')
        traj3plot(ax3d, l1['t'], l1['p'], l1['s'], "Blues_r", vscale=vscale)
        aspectEqual3(ax3d, l1['p'])
        traj3plot(ax3d, l2['t'], l2['p'], l2['s'], "Purples_r" if purple else "Reds_r", vscale=vscale)
        ax3d.plot(l1['posdes'][:,0], l1['posdes'][:,1], l1['posdes'][:,2], 'k--', alpha=0.5, zorder=9)
        ax3d.set_xlabel('x [mm]')
        ax3d.set_ylabel('y [mm]')
        ax3d.set_zlabel('z [mm]')
    def plott(ax, key, comp, lbl, dkey=None, dcomp=None):
        ax.plot(l1['t'], l1[key][:,comp], 'b')
        ax.plot(l2['t'], l2[key][:,comp], 'purple' if purple else 'r')
        if dkey is not None:
            ax.plot(l1['t'], l1[dkey][:,dcomp], 'k--', alpha=0.3)
        ax.set_ylabel(lbl)
    
    # pick 4 plots to show
    fig, ax = plt.subplots(2,2,figsize=(5,4))
    ax = ax.ravel()
    if plotset=='line':
        plott(ax[0], 'p', 0, 'x [mm]', 'posdes', 0)
        plott(ax[1], 'dp', 0, 'dx [mm]')
        plott(ax[2], 'omega', 1, 'omgy')
        plott(ax[3], 'eul', 1, 'euly')
    elif plotset=='helix_wlqp':
        plott(ax[0], 'p', 0, 'x [mm]', 'posdes', 0)
        plott(ax[1], 'p', 2, 'z [mm]', 'posdes', 2)
        plott(ax[2], 'omega', 1, 'omegay [rad/ms]')
        plott(ax[3], 'eul', 2, 'eulz [rad]')
    elif plotset=='hover_wlqp':
        plott(ax[0], 'p', 0, 'x [mm]', 'posdes', 0)
        plott(ax[1], 'p', 2, 'z [mm]', 'posdes', 2)
        plott(ax[2], 'omega', 1, 'omegay [rad/ms]')
        plott(ax[3], 'eul', 2, 'eulz [rad]')
    else:
        plott(ax[0], 'p', 0, 'x [mm]', 'posdes', 0)
        plott(ax[1], 'p', 2, 'z [mm]', 'posdes', 2)
        plott(ax[2], 'omega', 1, 'omegay [rad/ms]')
        plott(ax[3], 'eul', 2, 'eulz [rad]')
    fig.tight_layout()

if __name__ == "__main__":
    # data, ca6log = getData("")
    # defaultPlots(data, ca6log=ca6log)
    # plt.show()

    # For this one
    # papPlots('../logs/sdab_20201113124801.zip', '../logs/sdab_20201113124828.zip')
    # papPlots('../logs/sdab_20201115093555.zip', '../logs/sdab_20201115093651.zip', traj3d=True, vscale=100, plotset='line')

    # WLQP vs. not (same non-WLQP as above). second is WLQP 1,1.5
    papPlots('../logs/sdab_20201113124801.zip', '../logs/sdab_20201115120846.zip', purple=True, plotset='helix_wlqp')
    # papPlots('../logs/sdab_20201115141407.zip', '../logs/sdab_20201115141137.zip', purple=True, vscale=20, plotset='hover_wlqp') # force bias 3,0. No WLQP stopped at t=1800 (fails)
    plt.show()
