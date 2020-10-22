import numpy as np
import matplotlib.pyplot as plt

o0 = '9.2 26.3 72.8 149.0 281.6 493.1 800.4 1253.6 1871.6 2750.3'
o1 = '7.7 30.9 73.4 154.8 284.7 492.2 800.0 1264.6 1888.5 2719.4'
o2 = '3.8 7.3 14.1 24.3 38.3 56.8 81.5 113.2 152.6 201.7'
o3 = '3.7 6.2 11.6 19.7 31.8 47.1 67.9 94.9 130.2 174.7'
o0 = [float(i) for i in o0.split(' ')]
o1 = [float(i) for i in o1.split(' ')]
o2 = [float(i) for i in o2.split(' ')]
o3 = [float(i) for i in o3.split(' ')]
xx = [ i*i for i in range(10,101,10)]
fig2 = plt.figure(2)
ax=fig2.add_subplot(111)
ax.plot(xx,o0,color='b',marker='o')
ax.plot(xx,o1,color='r',marker='o')
ax.plot(xx,o2,color='g',marker='o')
ax.plot(xx,o3,color='y',marker='o')
ax.legend(['O0','O1','O2','O3'])
# ax.set_xlim(-1, T)
# ax.set_ylim(-1, N/4)
ax.set_xlabel('Problem size (N)')
ax.set_ylabel('Mean running time(ms)')
plt.savefig('performance.png')