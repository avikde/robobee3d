
import numpy as np
import gzip, pickle, glob, os
import matplotlib.pyplot as plt

def readFile(fname, dataFmtList=None):
    f1, ext = os.path.splitext(fname)
    zfile = gzip.GzipFile(fname, 'rb')
    data = pickle.load(zfile)
    # print(data)
    print('Opened '+fname+'; average data rate = ' + str(1.0/np.mean(np.diff(data['t'])))+'Hz')
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
    fig, ax = plt.subplots(3)
    ax[0].plot(data['t'], data['q'][:,:3])
    ax[0].set_ylabel('pos')

if __name__ == "__main__":
    data = getData("")
    defaultPlots(data)
    plt.show()
