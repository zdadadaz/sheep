from os import listdir
from os.path import isfile, join
from collections import defaultdict
import pandas as pd
outpath = './'

dictout = defaultdict(list)
dictout_num = defaultdict(list)
onlyfiles = [f for f in listdir(outpath) if isfile(join(outpath, f)) and f[-3:] == "out" and f[:5]=='slurm']
for f in onlyfiles:
    if f == 'slurm-6708450.out':
        continue
    with open(join(outpath, f), 'r') as f:
        script = f.read().split('\n')
        key = [script[1]] + script[0].split(' ')  + script[2].replace(',','').split(' ')
        key = ' '.join(key)

        if len(dictout[key])==0:
            for tidx in range(3,9):
                dictout[key].append(float(script[tidx].split(' ')[-2]))
                dictout_num[key].append(1)
        else:
            for tidx in range(3,9):
                dictout[key][tidx-3] += float(script[tidx].split(' ')[-2])
                dictout_num[key][tidx-3] += 1

# average 
dictlist = [i for i in dictout.keys()]
avgs = []
for key in dictlist:
    avg = [dictout[key][i]/dictout_num[key][i] for i in range(6)]
    # avg = ' '.join(str(avg)[1:-1].split(','))
    avgs.append(avg)
df = pd.DataFrame(columns=['exe', 'n', 'tpn', 't', 'ct', 'N','T','iS', 'iW','iG','totaltime','totalinit','act','eat','renew','time'],index=range(len(dictout.keys())))
for i in range(len(dictout.keys())):
    key = dictlist[i].split(' ')
    avg = avgs[i]
    df.iloc[i,0] = key[0]
    for j in range(1, 10):
        if key[(j-1)*2+2][-1] == ')':
            df.iloc[i,j] = key[(j-1)*2+2]
        else:
            df.iloc[i,j] = int(key[(j-1)*2+2])
    for j in range(10,16):
        df.iloc[i,j] = avg[j-10]
# print(df.sort_values(by='N'))
print(df[df['exe']=='sheep_hyb'].sort_values(by=['N','tpn','ct']))

