
import matplotlib.pyplot as plt
import numpy as np
import pickle
import math


path="D:/Young_Lee_Current/ZXH_PD/R-strike/"
filename1="GS_T.csv"
filename2="GS_F.csv"
savepath="D:/Young_Lee_Current/ZXH_PD/R-strike/GS/"


#   http://xsxwsq.ustc.edu.cn/

'''
with open(path+filename1) as f:
    line=f.readline()
    line=f.readline()
    x1=[]
    y1=[]
    while(line):
        line=line.strip()
        line=line.split(",")
        x1.append(float(line[0]))
        y1.append(float(line[1]))
        line=f.readline()

with open(path+filename2) as f:
    line=f.readline()
    line=f.readline()
    x2=[]
    y2=[]
    while(line):
        line=line.strip()
        line=line.split(",")
        x2.append(float(line[0]))
        y2.append(float(line[1]))
        line=f.readline()



f=open(savepath+"T_x.pkl", "wb")
pickle.dump(x1, f)
f.close()

f=open(savepath+"T_y.pkl", "wb")
pickle.dump(y1, f)
f.close()

f=open(savepath+"F_x.pkl", "wb")
pickle.dump(x2, f)
f.close()

f=open(savepath+"F_y.pkl", "wb")
pickle.dump(y2, f)
f.close()
'''












