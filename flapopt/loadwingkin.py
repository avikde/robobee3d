
import sys
from scipy.io import loadmat
import matplotlib.pyplot as plt

x = loadmat(sys.argv[1])
t = x['data'][:,0]
print(x['data'].shape)
print(x['currTest'])
print(x.keys())

fig, ax = plt.subplots(3)
ax[0].plot(t, x['data'][:,1])
ax[1].plot(t, x['data'][:,2])
ax[2].plot(t, x['data'][:,3])

plt.show()
