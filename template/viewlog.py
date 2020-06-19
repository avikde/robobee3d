
import numpy as np
import gzip, pickle, glob, os, time
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from ca6dynamics import M

def initLog():
    return {'t': [], 'q': [], 'dq': [], 'u': [], 'pdes': []}

def appendLog(data, t, q, dq, u, pdes):
    """rot can be any scipy rotation type"""
    lastT = -np.inf if len(data['t']) == 0 else data['t'][-1]
    if t - lastT >= 0.999:
        data['t'].append(t)
        data['q'].append(np.copy(q))
        data['dq'].append(np.copy(dq))
        data['u'].append(np.copy(u))
        data['pdes'].append(np.copy(pdes))
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
    return readFile(fname)

def defaultPlots(data):
    t = data['t']
    # SDAB log also has wings
    qb = data['q'][:,-7:]
    dqb = data['dq'][:,-6:]

    fig, ax = plt.subplots(6)
    ax[0].plot(data['t'], qb[:,:3])
    ax[0].set_ylabel('pos [mm]')

    eul = Rotation.from_quat(qb[:,3:]).as_euler('xyz')
    ax[1].plot(data['t'], eul)
    ax[1].set_ylabel('orn [rad]')

    ax[2].plot(data['t'], data['pdes'][:,2], 'k--', alpha=0.3, label='des')
    actMom = (M @ dqb.T).T
    ax[2].plot(data['t'], actMom[:,2], label='act')
    ax[2].set_ylabel('Momentum')
    ax[2].legend()

    if data['u'].shape[1] > 2:
        ax[3].plot(data['t'], data['u'][:,[0,3]])
        ax[3].set_ylabel('u1')
        # ca6 log
        ax[4].plot(data['t'], data['u'][:,[1,4]])
        ax[4].set_ylabel('u2')
        ax[5].plot(data['t'], data['u'][:,[2,5]])
        ax[5].set_ylabel('u3')
    else:
        # plot wing states
        qw = data['q'][:,:4]
        dqw = data['dq'][:,:4]
        ax[3].plot(t, qw[:,[0,2]])
        ax[3].set_ylabel('Stroke')
        ax[4].plot(t, qw[:,[1,3]])
        ax[4].set_ylabel('Pitch')

    ax[-1].set_xlabel('Time [ms]')

if __name__ == "__main__":
    data = getData("")
    defaultPlots(data)
    plt.show()
