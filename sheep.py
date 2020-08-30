import numpy as np
import h5py
import matplotlib.pyplot as plt
import math

N = 25
T = 100
FILE = "animal.h5"
DATASETs = "Ds_sheep"
DATASETw = "Ds_wolve"
DATASETg = "Ds_grass"
error = 0
hf = h5py.File(FILE, 'r')
sheep = hf[DATASETs]
wolve = hf[DATASETw]
grass = hf[DATASETg]
sheep = np.array(sheep)
wolve = np.array(wolve)
grass = np.array(grass)

grid_size = int(np.sqrt(N))
half = math.sqrt(N)

plt.rcParams['animation.html'] = 'jshtml'
fig2 = plt.figure(2)
ax=fig2.add_subplot(111)
plt.figure(3)
xx = []
yys = []
yyw = []
yyg = []

print(np.all(sheep==1))
for t in range(0,T):
    counts = 0
    countw = 0
    countg = 0
    xx.append(t)
    
    marker = [0]*N 
    base = t*N
    for i in range(N):
        y = i//half
        x = i - y*half
        if sheep[i+base] > error:
            counts += 1
            marker[i] = plt.scatter(x, y, color='blue', alpha=1, s=10**2)
            plt.annotate(round(sheep[i+base],2), (x, y))
        if wolve[i+base] > error:
            countw += 1
            marker[i] = plt.scatter(x, y, color='red', alpha=1, s=10**2)
            plt.annotate(round(wolve[i+base],2), (x, y))
        if grass[i+base] > 0:
            countg += 1
            marker[i] = plt.scatter(x, y, color='green', alpha=0.3, s=10**2,marker='s')
        else:
            marker[i] = plt.scatter(x, y, color='gray', alpha=0.3, s=10**2,marker='s')
    plt.xlim(-1, grid_size)
    plt.ylim(-1, grid_size)
    plt.grid(b=True, which='major', color='k', linestyle='--',alpha=0.1)
    plt.title('t: {}, sheep num: {}, wolve num: {},  grass num: {}'.format(str(t),str(counts),str(countw),str(countg)))
    
    yys.append(counts)
    yyw.append(countw)
    yyg.append(countg)
    ax.plot(xx,yys,color='b')
    ax.plot(xx,yyw,color='r')
    ax.plot(xx,yyg,color='g')
    ax.set_xlim(-1, T)
    ax.set_ylim(-1, N)
    ax.grid(b=True, which='major', color='k', linestyle='--',alpha=0.1)
    fig2.canvas.draw()

    plt.pause(1)
    plt.clf()

    

plt.show()