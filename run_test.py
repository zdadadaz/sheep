
import subprocess
from os import listdir
from os.path import isfile, join

set_dict = {'name':1,
            'nodes':2, 
            'taskpernode':3,
            'ntask':4,
            'omp':5,
            'mem':6
            }
exefile_dict = {
    0:'sheep',
    1:'sheep_omp',
    2:'sheep_mpi',
    3:'sheep_hyb'
}
mpi_line = 'time mpirun -mca btl ^openib ./'
omp_line = 'time ./'
template_path = './template.sh'
outpath = './scripts/'

def rw_template(path, outpath, input_setting, input_para,exefile):
    with open(path, 'r') as f:
        temp = f.read()
    temp = temp.split('\n')
    # print(temp.split('\n'))
    str_input_setting = [str(i) for i in input_setting]
    str_input_para = [str(i) for i in input_para]
    temp[1] = temp[1]+ '-'.join(str_input_setting) +'_'+ '_'.join(str_input_para)
    outpath = outpath+exefile +'-'+ '-'.join(str_input_setting) +'_'+ '_'.join(str_input_para) +'.sh'
    for i in range(2, 6):
        temp[i] = temp[i] + str(input_setting[i-2])
    temp[15] = "echo " + "\""+ exefile +"\"" 
    #parameters
    run_line = mpi_line
    if exefile == 'sheep_omp' or exefile == 'sheep':
        run_line = omp_line
    temp[16] = run_line + exefile +' ' + ' '.join(str_input_para)
    temp = '\n'.join(temp)
    # print(temp)
    with open(outpath, 'w') as f:
        f.write(temp)

def strong_test(input_para):
    nodes = [1 for i in range(5)] + [1 for i in range(5)] + [1 for i in range(5)]
    tpn = [1 for i in range(5)] + [1<<i for i in range(5)] +[1<<i for i in range(4)] + [12]
    nt = [1 for i in range(5)] + [1<<i for i in range(5)] + [1<<i for i in range(4)] + [12]
    omp = [1<<i for i in range(5)] + [1 for i in range(5)] +  [1<<i for i in range(3)] + [3,2]
    for exe in range(2,4):
        exefile = exefile_dict[exe]
        for i in range(len(nodes)):
            if exe <= 1 and i >4:
                continue
            if exe >1 and i <5:
                continue
            input_setting = [nodes[i], tpn[i], nt[i], omp[i]]
            rw_template(template_path, outpath, input_setting, input_para, exefile)
            if exe == 0:
                break

def soft_test():
    input_para = [[10000, 1000, 400, 200, 5000],  [40000, 1000, 1600, 800, 20000], \
                [90000, 1000,  4000, 2000, 45000],  [160000, 1000, 6400, 3200, 80000],  [250000, 1000, 8000, 4000, 125000]]
    nodes = [1 for i in range(5)] + [1 for i in range(5)]
    tpn = [1 for i in range(5)] + [1<<i for i in range(5)]
    nt = [1 for i in range(5)] + [1<<i for i in range(5)]
    omp = [1<<i for i in range(5)] + [1 for i in range(5)]
    for exe in range(0,4):
        exefile = exefile_dict[exe]
        for i in range(len(nodes)):
            _input = input_para[i%5]
            if exe <= 1 and i >4:
                continue
            if exe >1 and i <5:
                continue
            input_setting = [nodes[i], tpn[i], nt[i], omp[i]]
            rw_template(template_path, outpath, input_setting, _input, exefile)
            if exe == 0:
                break

# input_para = [10000, 1000, 400, 200, 5000]
# strong_test(input_para)

input_para = [160000, 1000, 25600, 12800, 80000]
strong_test(input_para)

soft_test()

# exefile = "sheep"
# n, tpn, nt, omp, 
# input_setting = [1,4,4,2]
# N, T, is, iw, ig 
# input_para = [10000, 1000, 400, 200, 5000]
# rw_template(template_path, outpath, input_setting, input_para, exefile)

onlyfiles = [f for f in listdir(outpath) if isfile(join(outpath, f))]
for f in onlyfiles:
    subprocess.run(["sbatch", outpath+f])
