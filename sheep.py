import numpy as np
import h5py
import matplotlib.pyplot as plt
import math

FILEs = "sheep.h5"
FILEw = "wolve.h5"
DATASET = "DS1"

error = 0.0001
hfs = h5py.File(FILEs, 'r')
hfw = h5py.File(FILEw, 'r')
sheep = hfs[DATASET]
wolve = hfw[DATASET]
sheep = np.array(sheep)
wolve = np.array(wolve)
grass = np.zeros_like(wolve)

N, T = sheep.shape
grid_size = int(np.sqrt(N))
t=1
half = math.sqrt(N)

plt.rcParams['animation.html'] = 'jshtml'
fig2 = plt.figure(2)
ax=fig2.add_subplot(111)
plt.figure(3)
# print(sheep[:][1].sum())
xx = []
yys = []
yyw = []
for t in range(0,T):
    counts = 0
    countw = 0
    countg = 0
    xx.append(t)
    
    marker = [0]*N 
    for i in range(N):
        y = i//half
        x = i - y*half
        if sheep[i,t] > error:
            counts += 1
            marker[i] = plt.scatter(x, y, color='blue', alpha=1, s=10**2)
            plt.annotate(round(sheep[i,t]), (x, y))
        if wolve[i,t] > error:
            countw += 1
            marker[i] = plt.scatter(x, y, color='red', alpha=1, s=10**2)
            plt.annotate(round(wolve[i,t]), (x, y))
        if grass[i,t] > error:
            countg += 1
            marker[i] = plt.scatter(x, y, color='green', alpha=1, s=10**2,marker='s')
            plt.annotate(round(wolve[i,t]), (x, y))
    plt.xlim(-1, grid_size)
    plt.ylim(-1, grid_size)
    plt.grid(b=True, which='major', color='k', linestyle='--',alpha=0.1)
    plt.title('t: {}, sheep num: {}, wolve num: {},  grass num: {}'.format(str(t),str(counts),str(countw),str(countg)))
    
    yys.append(counts)
    yyw.append(countw)
    ax.plot(xx,yys,color='b')
    ax.plot(xx,yyw,color='r')
    ax.set_xlim(-1, T)
    ax.set_ylim(-1, N)
    ax.grid(b=True, which='major', color='k', linestyle='--',alpha=0.1)
    fig2.canvas.draw()

    plt.pause(1)
    plt.clf()

    

plt.show()