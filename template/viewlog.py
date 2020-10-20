
import numpy as np
import gzip, pickle, glob, os, time
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from ca6dynamics import M

def initLog():
    return {'t': [], 'q': [], 'dq': [], 'u': [], 'accdes': []}

def appendLog(data, t, q, dq, u, accdes):
    """rot can be any scipy rotation type"""
    lastT = -np.inf if len(data['t']) == 0 else data['t'][-1]
    if t - lastT >= 0.999:
        data['t'].append(t)
        data['q'].append(np.copy(q))
        data['dq'].append(np.copy(dq))
        data['u'].append(np.copy(u))
        data['accdes'].append(np.copy(accdes))
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

def readFile(fname, dataFmtList=None):
    f1, ext = os.path.splitext(fname)
    zfile = gzip.GzipFile(fname, 'rb')
    data = pickle.load(zfile)
    # print(data)
    print('Opened '+fname+'; average data rate = ' + str(1000.0/np.mean(np.diff(data['t'])))+'Hz')
    return data

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
    qb = data['q'][:,-7:]
    dqb = data['dq'][:,-6:]

    fig, ax = plt.subplots(3,2)
    ax = ax.ravel()
    ax[0].plot(data['t'], qb[:,:3])
    ax[0].set_ylabel('pos [mm]')

    s = np.zeros((len(t), 3))
    for i in range(len(t)):
        s[i,:] = Rotation.from_quat(qb[i,3:]).as_matrix()[:,2]
    ax[1].plot(data['t'], s)
    ax[1].set_ylabel('s')
    # eul = Rotation.from_quat(qb[:,3:]).as_euler('xyz')
    # ax[1].plot(data['t'], eul)
    # ax[1].set_ylabel('orn [rad]')

    # ax[2].plot(data['t'], data['accdes'][:,:3])
    ax[2].plot(data['t'], data['accdes'][:,3:])
    # actMom = (M @ dqb.T).T
    # ax[2].plot(data['t'], actMom[:,2], label='act')
    ax[2].set_ylabel('Accdes ang')
    # ax[2].legend()

    if ca6log:
        ax[3].plot(data['t'], data['u'][:,[0,3]])
        ax[3].set_ylabel('u1')
        # ca6 log
        ax[4].plot(data['t'], data['u'][:,[1,4]])
        ax[4].set_ylabel('u2')
        ax[5].plot(data['t'], data['u'][:,[2,5]])
        ax[5].set_ylabel('u3')
    else:
        ax[3].plot(data['t'], dqb[:,3:])
        ax[3].set_ylabel('Omega')
        # Inputs
        ax[4].plot(t, data['u'][:,2]) # Vmean
        ax[4].set_ylabel('Vmean')
        ax[5].plot(t, data['u'][:,3], label='offs')
        ax[5].plot(t, data['u'][:,4], label='diff')
        ax[5].plot(t, data['u'][:,5], label='h2')
        ax[5].set_ylabel('u')
        ax[5].legend()
        # # plot wing states
        # qw = data['q'][:,:4]
        # dqw = data['dq'][:,:4]
        # ax[5].plot(t, qw[:,[0,2]])
        # ax[5].set_ylabel('Stroke')
        # ax[6].plot(t, qw[:,[1,3]])
        # ax[6].set_ylabel('Pitch')

    ax[-1].set_xlabel('Time [ms]')
    fig.tight_layout()

if __name__ == "__main__":
    data, ca6log = getData("")
    defaultPlots(data, ca6log=ca6log)
    plt.show()
