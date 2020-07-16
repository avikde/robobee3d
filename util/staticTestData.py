import csv, sys, itertools
import numpy as np
import matplotlib.pyplot as plt

def loadStaticTestData(fname):
    csvfile = open(fname, newline='')
    dat = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    # create dict with voltage as key
    freq = {}
    strokes = {}
    # create keys
    for row in dat:
        freq[row['Voltage']] = []
        strokes[row['Voltage']] = []
    # populate the data
    csvfile.seek(0)
    dat = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    for row in dat:
        freq[row['Voltage']].append(float(row['Frequency']))
        strokes[row['Voltage']].append(np.hstack((float(row['Stroke bot']), float(row['Stroke top']))))
    # convert to np array from list
    for key in freq.keys():
        freq[key] = np.array(freq[key])
        strokes[key] = np.array(strokes[key])
    return freq, strokes

def normStrokePlot(ax, freq, strokes, avg=False):
    marker = itertools.cycle(('^', '+', '.', 'o', '*', 'v')) 
    for key in freq.keys():
        V = float(key)
        if avg:
            ax.plot(freq[key], np.mean(strokes[key], axis=1) / V, marker=next(marker), label=key+'V bot')
        else:
            ax.plot(freq[key], strokes[key][:,0] / V, marker=next(marker), label=key+'V bot')
            ax.plot(freq[key], strokes[key][:,1] / V, marker=next(marker), label=key+'V top')
    ax.set_xlabel('Freq [Hz]')
    ax.set_ylabel('Norm stroke [deg/V]')

aa = loadStaticTestData(sys.argv[1])
fig, ax = plt.subplots()
normStrokePlot(ax, *aa)
ax.legend()
plt.show()


