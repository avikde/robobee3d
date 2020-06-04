
import numpy as np
import gzip, pickle, glob, os, time
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation

def initLog():
    return {'t': [], 'q': [], 'v': [], 'u': []}

def appendLog(data, t, pw, Rb, dq, u):
    data['t'].append(t)
    data['q'].append(np.hstack((pw, Rb.as_euler('xyz'))))
    data['v'].append(dq)
    data['u'].append(u)
    return data

def saveLog(f1, data):
    for k in data.keys():
        data[k] = np.array(data[k])

    t = time.localtime()
    timestamp = time.strftime('%Y%m%d%H%M%S', t)
    fname = f1 + 'ca6_' + timestamp + '.pkz'
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
    fig, ax = plt.subplots(2)
    ax[0].plot(data['t'], data['q'][:,:3])
    ax[0].set_ylabel('pos [mm]')

    Nt = len(data['t'])
    # eul = np.zeros((Nt,3))
    # for ti in range(Nt):
    #     eul[ti,:] = Rotation.from_quat(data['q'][ti,:]).as_euler()
    #     # ax[1].plot(data['t'])

    ax[1].plot(data['t'], data['q'][:,3:])
    ax[1].set_ylabel('orn [rad]')
    ax[1].set_xlabel('Time [ms]')

if __name__ == "__main__":
    data = getData("")
    defaultPlots(data)
    plt.show()
