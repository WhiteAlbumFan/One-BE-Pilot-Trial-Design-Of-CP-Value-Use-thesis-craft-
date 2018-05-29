'''
    import numpy as np
    import pickle
    import matplotlib.pyplot as plt
    import csv
    import pylab as pyb
    from scipy.stats import gaussian_kde
    from matplotlib.colors import LogNorm

    from mpl_toolkits.mplot3d import Axes3D
'''

import matplotlib.pyplot as plt
import numpy as np
import pickle
#from mpl_toolkits.mplot3d import Axes3D
import math

path="D:/Young_Lee_Current/ZXH_PD/R-strike/GS/"
savepath="D:/Young_Lee_Current/ZXH_PD/R-strike/GS/"

for ks in [1]:
    f=open(path+"T_x.pkl", "rb")
    x1=pickle.load(f)
    f.close()

    f=open(path+"T_y.pkl", "rb")
    y1=pickle.load(f)
    f.close()


    f=open(path+"F_x.pkl", "rb")
    x2=pickle.load(f)
    f.close()

    f=open(path+"F_y.pkl", "rb")
    y2=pickle.load(f)
    f.close()

'''
    f=open(path+"z_kk.pkl", "rb")
    z1=pickle.load(f)
    f.close()

    f=open(path+"z_kk_fake.pkl", "rb")
    z2=pickle.load(f)
    f.close()


    plt.figure(figsize=(25, 15), dpi=320 )
    plt.scatter(x2, y2, color="b", marker="+", s=5)
    plt.scatter(x1, y1, color="r", marker="x", s=5)


    #plt.show()
    plt.savefig(savepath+"trues.png")
'''






'''
def density_olds(x, y, slice=200):
    length_s=min( len(x), len(y))
    wait=[]
    rec=np.zeros(shape=(slice+1,slice+1), dtype=int )
    for k in range(length_s):
        x1=int(x[k]*slice)
        y1=int(y[k]*slice)
        rec[x1, y1]+=1
    for k in range(length_s):
        if (math.fabs(x[k]-1)<1e-3) & (math.fabs(y[k])<1e-3):
            wait.append(0)
        else:
            x1=int(x[k]*slice)
            y1=int(y[k]*slice)
            temp=float(rec[x1, y1])
            if temp>50:
                temp=0
            wait.append(temp)
    return wait

z1=density_olds(x1, y1, slice=500)
z2=density_olds(x2, y2, slice=500)

f=open(savepath+"T_z.pkl", "wb")
pickle.dump(z1, f)
f.close()

f=open(savepath+"F_z.pkl", "wb")
pickle.dump(z2, f)
f.close()

theta=""
f=open(savepath+"T.csv", "w")
for k in range(len(z1)):
    theta=theta+str(x1[k])+","+str(y1[k])+","+str(z1[k])+"\n"
f.write(theta)
f.close()


theta=""
f=open(savepath+"F.csv", "w")
for k in range(len(z2)):
    theta=theta+str(x2[k])+","+str(y2[k])+","+str(z2[k])+"\n"
f.write(theta)
f.close()
'''




