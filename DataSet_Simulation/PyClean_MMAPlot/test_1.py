import numpy as np
'''
import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde

x = np.random.normal(size=1000)
y = x *3+ np.random.normal(size=1000)


xy = np.vstack([x,y])

z = gaussian_kde(xy)(xy)

fig, ax = plt.subplots()

ax.scatter(x, y, c=z, s=100, edgecolor='')

plt.show()
'''
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
 
# 创建 3D 图形对象
fig = plt.figure()
ax = Axes3D(fig)
 
# 生成数据
X = np.arange(-2, 2, 0.1)
Y = np.arange(-2, 2, 0.1)
X, Y = np.meshgrid(X, Y)
Z = np.sqrt(X ** 2 + Y ** 2)
 
# 绘制曲面图，并使用 cmap 着色
ax.plot_surface(X, Y, Z, cmap=plt.cm.winter)
 
plt.show()

