import csv, sys
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
        strokes[row['Voltage']].append(np.hstack((float(row['Stroke bot']), float(row['Stroke bot']))))
    # convert to np array from list
    for key in freq.keys():
        freq[key] = np.array(freq[key])
        strokes[key] = np.array(strokes[key])
    return freq, strokes

# def normStrokePlot(ax, V, f, S):
#     ax.

f,S = loadStaticTestData(sys.argv[1])
print(f, S)


