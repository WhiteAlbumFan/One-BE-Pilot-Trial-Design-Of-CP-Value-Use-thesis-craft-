

'''
fig = plt.figure() 
axss= fig.add_subplot(111, projection='3d') 
#axss.plot_trisurf(x1, y1, z1) 
axss.plot_trisurf(x1, y1, z1, color="r")
axss.plot_trisurf(x2, y2, z2, color="b")
plt.savefig(savepath+'3d_fig.png',dpi=400)


plt.show()
'''




'''
print(z1[0:100])
fig, ax = plt.subplots()

ax.scatter(x1, y1, c=z1, s=1, edgecolor='')

plt.show()
'''

'''
def fun(x,y):  
    return np.power(x,2)+np.power(y,2)  


x1=x1[0:100]
y1=y1[0:100]
fig1=plt.figure()           #创建一个绘图对象  
ax=Axes3D(fig1)             #用这个绘图对象创建一个Axes对象(有3D坐标)  

x1,y1=np.meshgrid(x1,y1)  
z=fun(x1,y1)

ax.plot_surface(x1, y1, z, rstride=1, cstride=1, cmap="rainbow")  #用取样点(x,y,z)去构建曲面  
ax.set_xlabel('x label', color='r')  
ax.set_ylabel('y label', color='g')  
ax.set_zlabel('z label', color='b')         #给三个坐标轴注明  
plt.savefig(savepath+"GS.png")
plt.show()                                  #显示模块中的所有绘图对象  
'''

'''
xy = np.vstack([x1,y1])

z = gaussian_kde(xy)(xy)

fig, ax = plt.subplots()

ax.scatter(x1, y1, c=z, s=1, edgecolor='', marker=",")

plt.show()

'''





'''

with open(path+filename) as f:
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
filename="kk_fake.csv"
with open(path+filename) as f:
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



f=open(path+"kk_x.pkl", "wb")
pickle.dump(x1, f)
f.close()

f=open(path+"kk_y.pkl", "wb")
pickle.dump(y1, f)
f.close()

f=open(path+"kk_fake_x.pkl", "wb")
pickle.dump(x2, f)
f.close()

f=open(path+"kk_fake_y.pkl", "wb")
pickle.dump(y2, f)
f.close()
'''





