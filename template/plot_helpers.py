
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection

def traj3plot(_ax, t, p, v, cmap, narrow=20, vscale=0.4):
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
