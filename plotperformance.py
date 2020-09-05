import numpy as np
import matplotlib.pyplot as plt

o0 = '31.5	123.7	307.6	492.8	820	1114	1429.5	1938.2	2406.4	3035.7'
o1 = '44.9	116.4	286.4	479.6	654.3	932.4	1350.5	1754.2	2382.6	2966.7'
o2 = '43	166.3	276.9	444.6	761.3	1118.3	1399.9	1884.2	2354.1	2756.2'
o3 = '25.8	101.3	260.8	438.7	711.7	1009.8	1449.1	1849.9	2084.1	2644.8'
o0 = [float(i) for i in o0.split('\t')]
o1 = [float(i) for i in o1.split('\t')]
o2 = [float(i) for i in o2.split('\t')]
o3 = [float(i) for i in o3.split('\t')]
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