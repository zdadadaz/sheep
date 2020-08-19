#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#define FILE "sheep.h5"
#define DATASET "DS1"
#define N 100
#define T 100
#define initSheepNum 10
#define sheepGainFromFood 4
#define sheepReproduce 4 //%
#define initWolveNum 5
#define wolveGainFromFood 20
#define wolveReproduce 5 //%
#define Grass 0

int half = 10;
int tot_sheep=0;
int tot_wolve=0;

int mat2hdf5(double wdata[N][T])
{
    hid_t file, space, dset; /* Handles */
    herr_t status;
    hsize_t dims[2] = {N, T};

    /*
     * Create a new file using the default properties.
     */
    file = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Create dataspace.  Setting maximum size to NULL sets the maximum
     * size to be the current size.
     */
    space = H5Screate_simple(2, dims, NULL);

    /*
     * Create the dataset and write the floating point data to it.  In
     * this example we will save the data as 64 bit little endian IEEE
     * floating point numbers, regardless of the native type.  The HDF5
     * library automatically converts between different floating point
     * types.
     */
    dset = H5Dcreate(file, DATASET, H5T_IEEE_F64LE, space, H5P_DEFAULT,
                     H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      wdata[0]);

    /*
     * Close and release resources.
     */
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Fclose(file);

    return 0;
}

void gen_surroundxy(int curX, int curY, int xlist[8], int ylist[8]){
    int xl[8] = {-1, -1 ,0 ,1, 1,1,0,-1};
    int yl[8] = {0  ,-1,-1,-1, 0,1,1, 1};
    for (int i = 0; i < 8; i++){
        xlist[i] = curX + xl[i];
        ylist[i] = curY + yl[i];
    }
}

bool check_move_valid(double animal[N][T], int curX, int curY, int x, int y, int t, int t2)
{
    if ((x >= 0) && (x < half) && (y >= 0) && (y < half) && animal[x + y * half][t]==0)
        animal[x + y * half][t2] = animal[curX + curY * half][t];
        if (t2 == t)
            animal[curX + curY * half][t] = 0;
        return true;
    return false;
}
bool check_sheep_exist(double animal[N][T], int x, int y, int t)
{
    if ((x >= 0) && (x < half) && (y >= 0) && (y < half) && animal[x + y * half][t] > 0)
        return true;
    return false;
}

void move(double animal[N][T], int curX, int curY, int t,int t2, int* newX, int* newY)
{
    int flag = false;
    int tmpX = curX, tmpY = curY;
    int count = 0;
    int xlist[8] = {0}, ylist[8] = {0};
    gen_surroundxy(curX, curY, xlist, ylist);
    for (int i = 0; i < 8; i++)
    {
        int randint = (float)rand() / (float)(RAND_MAX)*7;
        tmpX = xlist[randint];
        tmpY = ylist[randint];
        flag = check_move_valid(animal, curX, curY, tmpX, tmpY, t, t2);
        if (flag == true){
            *newX = tmpX;
            *newY = tmpY;
            break;
        }
    }
    if(flag == false)
    {
        *newX = curX;
        *newY = curY;
    }
}
void death(double animal[N][T], int curX, int curY, int t, int flag){
    if (animal[curX + curY* half][t]<0){
        animal[curX + curY * half][t] = 0;
        if (flag == 0)
            tot_sheep -= 1;
        else
            tot_wolve -= 1;
    }
}
void create_animal(double animal[N][T], int curX, int curY, int t, int flag){
    int gainFood = flag == 0 ? sheepGainFromFood : wolveGainFromFood;
    int Energy = 2 * gainFood * (float)rand() / (float)(RAND_MAX);
    int curIdx = curX + curY *half;
    int newX, newY;
    int tmpEnergy;

    tmpEnergy = animal[curIdx][t];
    animal[curIdx][t] = Energy;
    move(animal, curX, curY, t, t, &newX, &newY);
    animal[curIdx][t] = tmpEnergy;
    if (flag == 0)
        tot_sheep += 1;
    else
        tot_wolve += 1;
}

void reproduce(double animal[N][T], int curX, int curY, int t, int flag)
{
    int randint = (float)rand() / (float)(RAND_MAX)*100;
    int cond_reproduce_rate = (flag == 0) ? sheepReproduce : wolveReproduce;
    if (randint < cond_reproduce_rate){
        animal[curX + curY * half][t] /= 2;
        create_animal(animal, curX, curY, t, flag);
    }
}

void init_sheep_wolve(double sheep[N][T], double wolve[N][T])
{
    int sheepNum = 0;
    while (sheepNum < initSheepNum){
        int randint = (float)rand() / (float)(RAND_MAX) * (N - 1);
        if (sheep[randint][0] == 0){
            sheep[randint][0] = 2 * sheepGainFromFood * (float)rand() / (float)(RAND_MAX);
            sheepNum++;
            tot_sheep++;
        }
    }
    int wolveNum = 0;
    while (wolveNum < initWolveNum)
    {
        int randint = (float)rand() / (float)(RAND_MAX) * (N - 1);
        if (wolve[randint][0] == 0){
            wolve[randint][0] = 2 * wolveGainFromFood * (float)rand() / (float)(RAND_MAX);
            wolveNum++;
            tot_wolve++;
        }
    }
}
void catch_sheep(double wolve[N][T], double sheep[N][T], int *curX, int *curY, int tw, int ts)
{
    int count = 0;
    int xlist[8] = {0}, ylist[8] = {0};
    int count_list[8] = {0};
    int idx = (*curX) + (*curY) * half;
    gen_surroundxy(*curX, *curY, xlist, ylist);
    for (int i = 0; i < 8; i++)
    {
        bool check = check_sheep_exist(sheep, xlist[i], ylist[i], ts);
        if (check == true){
            count_list[i] = 1;
            count+=1;
        }
    }
    if (count > 0){
        int randint = (float)rand() / (float)(RAND_MAX) * (count - 1);
        sheep[xlist[randint] + half * ylist[randint]][ts] = 0; //kill
        tot_sheep += 1;
        wolve[idx][tw] += wolveGainFromFood;
    }
}

void ask_sheep(double sheep[N][T], int i, int t){
    int curY = i/half;
    int curX = i - curY * half;
    int newX, newY;
    move(sheep, curX, curY, t, t+1, &newX, &newY);
    // if Grass == 1{
    //     sheep[newX + newY * half][t] = sheep[newX + newY * half][t]-1
    // }
    // death(sheep, newX, newY, t+1, 0);
    // reproduce(sheep, newX, newY, t+1, 0);
}
void ask_wolve(double wolve[N][T], double sheep[N][T], int i, int t)
{
    int curY = i / half;
    int curX = i - curY * half;
    int newX, newY;
    move(wolve, curX, curY, t, t+1, &newX, &newY);
    wolve[i][t+1] -= 1;
    // catch_sheep(wolve, sheep, &newX, &newY, t+1, t+1);
    death(wolve, newX, newY, t+1, 1);
    reproduce(wolve, newX, newY, t+1, 0);
}

int main(void)
{
    double sheep[N][T] = {0}, wolve[N][T] = {0};
    init_sheep_wolve(sheep, wolve);
    for (int t = 0; t < T-1; t++)
    {
        printf("sheep = %d\n", tot_sheep);
        printf("wolve = %d\n", tot_wolve);
        for (int i = 0; i < N; i++)
        {
                ask_sheep(sheep, i, t);
    //             // ask_wolve(wolve, sheep, i, t);
    //             // if (Grass == 1)
    //             //     ask_patch(grass,i,t)
                
        }
    }

    mat2hdf5(sheep);
    // mat2hdf5(wolve);
    return 0;
}
