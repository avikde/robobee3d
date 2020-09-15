
from wlqppy import WLController
import matplotlib.pyplot as plt

wl = WLController()

pdotdes = [0, 0, 10, 0, 0, 0]
mb = 100
g = 9.81e-3
h0 = [0, 0, mb * g, 0, 0, 0]
u0 = [140.0, 0., 0., 0.]

u = wl.update(u0, h0, pdotdes)
print(u)
