import numpy as np
import h5py
import matplotlib.pyplot as plt
import math
from util import *
path = './animal.h5'
sheep, wolve, grass, count, setting = read_hdf(path)
N = setting[0]
T = setting[1]
# setting
realW = int(math.sqrt(N))
fig2 = plt.figure(2)
ax=fig2.add_subplot(111)
subsample = 1
xx = []
yys =[]
yyw = []
yyg = []
for t in range(0,T,subsample):
    xx.append(t)
    yys.append(count[t*3])
    yyw.append(count[1+t*3])
    yyg.append(count[2+t*3]/4)
    
ax.plot(xx,yys,color='b')
ax.plot(xx,yyw,color='r')
ax.plot(xx,yyg,color='g')
ax.legend(['sheep','wolves','grass'])
ax.set_xlim(-1, T)
ax.set_ylim(-1, N/4)
ax.set_xlabel('Time (T)')
ax.set_ylabel('Population (N)')
ax.grid(b=True, which='major', color='k', linestyle='--',alpha=0.1)
plt.savefig('population_dynamic.png')