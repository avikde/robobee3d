import csv, sys
import numpy as np

def loadStaticTestData(fname):
    csvfile = open(fname, newline='')
    dat = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    # convert to lists for each key
    voltage = []
    freq = []
    strokes = []
    for row in dat:
        voltage.append(float(row['Voltage']))
        freq.append(float(row['Frequency']))
        strokes.append(np.hstack((float(row['Stroke bot']), float(row['Stroke bot']))))
    return np.array(voltage), np.array(freq), np.array(strokes)

V, f, S = loadStaticTestData(sys.argv[1])
print(V, f, S)