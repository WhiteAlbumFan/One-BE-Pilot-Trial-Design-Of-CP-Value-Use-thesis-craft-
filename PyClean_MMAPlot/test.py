from matplotlib.colors import LogNorm

from pylab import *


x = randn(100000)
y = randn(100000)+5 

hist2d(x, y, bins=40, norm=LogNorm())

colorbar()

show()