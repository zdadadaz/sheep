#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#define FILE "sheep.h5"
#define FILEw "wolve.h5"
#define FILEg "grass.h5"
#define DATASET "DS1"
#define N 100
#define T 100
#define initSheepNum 10
#define sheepGainFromFood 4
#define sheepReproduce 4 //%
#define initWolveNum 4
#define wolveGainFromFood 20
#define wolveReproduce 5 //%
#define Grass 1
#define initGrass 25
#define grassRegrowth 50 //%
#define max(a, b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define error 0.00001

int half = sqrt(N);
int tot_sheep=0;
int tot_wolve = 0;
int tot_grass = 0;

int mat2hdf5(double wdata[N][T], char *filename) 
// int mat2hdf5(double **wdata, char *filename)
{
    hid_t file, space, dset; /* Handles */
    herr_t status;
    hsize_t dims[2] = {N, T};

    file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    space = H5Screate_simple(2, dims, NULL);
    dset = H5Dcreate(file, DATASET, H5T_IEEE_F64LE, space, H5P_DEFAULT,
                     H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      wdata[0]);
    status = H5Dclose(dset);
    status = H5Sclose(space);
    status = H5Fclose(file);

    return 0;
}
int mat2hdf5int(int wdata[N][T], char *filename)
// int mat2hdf5(double **wdata, char *filename)
{
    hid_t file, space, dset; /* Handles */
    herr_t status;
    hsize_t dims[2] = {N, T};

    file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    space = H5Screate_simple(2, dims, NULL);
    dset = H5Dcreate(file, DATASET, H5T_IEEE_F64LE, space, H5P_DEFAULT,
                     H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      wdata[0]);
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
//  bool check_move_valid(double **animal, int curX, int curY, int x, int y, int t, int t2)
{
    if ((x >= 0) && (x < half) && (y >= 0) && (y < half) && (animal[x + y * half][t2] < error) && (animal[x + y * half][t] < error))
    {
        animal[x + y * half][t2] = animal[curX + curY * half][t];
        if (t2 == t)
            animal[curX + curY * half][t] = 0;
        return true;
    }
    return false;
}
bool check_sheep_exist(double animal[N][T], int x, int y, int t) 
// bool check_sheep_exist(double **animal, int x, int y, int t)
{
    if ((x >= 0) && (x < half) && (y >= 0) && (y < half) && animal[x + y * half][t] > error)
        return true;
    return false;
}

void move(double animal[N][T], int curX, int curY, int t, int t2, int *newX, int *newY) 
// void move(double **animal, int curX, int curY, int t, int t2, int *newX, int *newY)
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
        // printf("%d,%d\n", tmpX, tmpY);
        flag = check_move_valid(animal, curX, curY, tmpX, tmpY, t, t2);
        if (flag == true){
            // printf("tmpX %d,tmpY %d\n", tmpX, tmpY);
            *newX = tmpX;
            *newY = tmpY;
            break;
        }
    }
    // printf("t1 %d, t2 %d, newX %d,newY %d\n", t, t2, curX, curY);
    // stay the same place
    if (flag == false)
    {
        if (animal[curX + curY * half][t2] <= error)
        {
            // printf(" t1 %d, t2 %d, stay\n", t, t2);
            *newX = curX;
            *newY = curY;
            animal[curX + curY * half][t2] = animal[curX + curY * half][t];
        }
        else
        {   
            int count = 0, cnt = 0;
            while ((count < 1) && cnt < N)
            {
                int randint = (float)rand() / (float)(RAND_MAX) * (N - 1);
                if (animal[randint][t2] <= 0)
                {
                    animal[randint][t2] = animal[randint][t];
                    count++;
                    printf("fly\n");
                }
                cnt++;
            }
            if (count <1)
                printf("all occupied\n");
        }
    }
}
void death(double animal[N][T], int curX, int curY, int t, int flag){
// void death(double** animal, int curX, int curY, int t, int flag){
    if (animal[curX + curY * half][t] <= error)
    {
        animal[curX + curY * half][t] = 0;
        if (flag == 0)
            tot_sheep -= 1;
        else
            tot_wolve -= 1;
    }
}
void create_animal(double animal[N][T], int curX, int curY, int t, int flag){
// void create_animal(double** animal, int curX, int curY, int t, int flag){
    int gainFood = flag == 0 ? sheepGainFromFood : wolveGainFromFood;
    int Energy = 2 * gainFood * (float)rand() / (float)(RAND_MAX);
    int count = 0, cnt = 0;
    int idx = curX + curY *half;
    while ((count < 1) && cnt < N)
    {
        if (animal[cnt][t] <= error)
        {
            animal[cnt][t] = max(1.0, Energy);
            count++;
        }
        cnt++;
    }
    if (count>0){
        if (flag == 0)
            tot_sheep += 1;
        else
            tot_wolve += 1;
    }else{
        printf("no place for reproduce flag %d\n",flag);
    }
}

void reproduce(double animal[N][T], int curX, int curY, int t, int flag) 
// void reproduce(double **animal, int curX, int curY, int t, int flag)
{
    float randint = (float)rand() / (float)(RAND_MAX)*100;
    int cond_reproduce_rate = (flag == 0) ? sheepReproduce : wolveReproduce;
    if (randint < (float)cond_reproduce_rate){
        animal[curX + curY * half][t] /= 2;
        create_animal(animal, curX, curY, t, flag);
    }
}

void init_sheep_wolve(double sheep[N][T], double wolve[N][T], int grass[N][T])
// void init_sheep_wolve(double **sheep, double **wolve, double **grass)
{
    int sheepNum = 0;
    while (sheepNum < initSheepNum){
        int randint = (float)rand() / (float)(RAND_MAX) * (N-1);
        if (sheep[randint][0] < error){
            sheep[randint][0] = max(1.0, 2 * sheepGainFromFood * (float)rand() / (float)(RAND_MAX));
            // printf("randint %d, sheepNum %d, sheep %f, out %f \n", randint, sheepNum, sheep[randint][0], 2 * sheepGainFromFood * (float)rand() / (float)(RAND_MAX));
            sheepNum++;
            tot_sheep++;
        }
    }
    int wolveNum = 0;
    while (wolveNum < initWolveNum)
    {
        int randint = (float)rand() / (float)(RAND_MAX) * (N-1);
        if (wolve[randint][0] < error)
        {
            wolve[randint][0] = max(1.0, 2 * wolveGainFromFood * (float)rand() / (float)(RAND_MAX));
            wolveNum++;
            tot_wolve++;
        }
    }
    if (Grass == 0){
        for (int i =0; i<N; i++)
            grass[i][0] = 1;
    }else{
        int grassNum = 0;
        while (grassNum < initGrass)
        {
            int randint = (float)rand() / (float)(RAND_MAX) * (N - 1);
            if (grass[randint][0] < error)
            {
                grass[randint][0] = 1;
                grassNum++;
                tot_grass++;
            }
        }
    }
}
void catch_sheep(double wolve[N][T], double sheep[N][T], int *curX, int *curY, int tw, int ts) 
// void catch_sheep(double **wolve, double **sheep, int *curX, int *curY, int tw, int ts)
{
    int count = 0;
    int xlist[8] = {0}, ylist[8] = {0};
    int count_list[9] = {0};
    int idx = (*curX) + (*curY) * half;
    gen_surroundxy(*curX, *curY, xlist, ylist);
    
    for (int i = 0; i < 8; i++)
    {
        bool check = check_sheep_exist(sheep, xlist[i], ylist[i], ts);
        if (check == true){
            count_list[i] = 1;
            count++;
        }
    }
    if (sheep[idx][ts] > error)
    {
        count_list[8] = 1;
        count++;
    }

    if (count > 0)
    {
        // int randint = (float)rand() / (float)(RAND_MAX) * (count - 1);
        int randint = 8;
        for (int i = 8; i >= 0;i--){
            if (count_list[i] > 0)
            {
                randint = i;
                break;
            }
        }
        sheep[xlist[randint] + half * ylist[randint]][ts] = 0; //kill
        tot_sheep--;
        wolve[idx][tw] += wolveGainFromFood;
    }
}
void eatgrass(double sheep[N][T], int grass[N][T], int i, int t)
{
    // t+1 sheep eat t grass
    if (grass[i][t]>0){
        sheep[i][t + 1] += sheepGainFromFood;
        grass[i][t] = 0;
    }
}

void ask_patch(int grass[N][T], int i, int t)
{
    if (grass[i][t]==1){
        grass[i][t + 1] = 1;
    }
    else{
        float randint = (float)rand() / (float)(RAND_MAX)*100;
        if (randint < (float)grassRegrowth)
            grass[i][t+1] = 1;
        else
            grass[i][t + 1] = 0;
    }
}

void ask_sheep(double sheep[N][T], int grass[N][T], int i, int t)
{
    // void ask_sheep(double** sheep,double** grass, int i, int t){
    if (sheep[i][t] < error)
        return;
    int curY = i / half;
    int curX = i - curY * half;
    int newX, newY;
    move(sheep, curX, curY, t, t+1, &newX, &newY);
    if (Grass == 1) 
    {
        int newIdx = newX + newY * half;
        sheep[newIdx][t+1]--;
        eatgrass(sheep, grass, newIdx, t); // t+1 sheep eat t grass
    }
    death(sheep, newX, newY, t+1, 0);
    reproduce(sheep, newX, newY, t+1, 0);
}
void ask_wolve(double wolve[N][T], double sheep[N][T], int i, int t) 
// void ask_wolve(double **wolve, double **sheep, int i, int t)
{
    if (wolve[i][t] < error)
        return;
    int curY = i / half;
    int curX = i - curY * half;
    int newX, newY;
    move(wolve, curX, curY, t, t+1, &newX, &newY);
    wolve[newX + newY*half][t + 1] -= 1;
    catch_sheep(wolve, sheep, &newX, &newY, t+1, t+1);
    death(wolve, newX, newY, t+1, 1);
    reproduce(wolve, newX, newY, t+1, 0);
}

int main(void)
{
    // double **sheep = malloc(N * sizeof(double *));
    // double **wolve = malloc(N * sizeof(double *));
    // for (int i = 0; i < N; i++){
    //     sheep[i] = malloc(T * sizeof(double));
    //     wolve[i] = malloc(T * sizeof(double));
    // }
    double sheep[N][T] = {0};
    double wolve[N][T] = {0};
    int grass[N][T] = {0};

    init_sheep_wolve(sheep, wolve, grass);
    // for (int t = 0; t < T; t++)
    // {
    //     for (int i = 0; i < N; i++)
    //     {
    //         sheep[i][t] = 1;
    //         wolve[i][t] = 1;
    //     }
    // }
    
    for (int t = 0; t < T-1; t++)
    {
        //     // printf("sheep = %d\n", tot_sheep);
        //     // printf("wolve = %d\n", tot_wolve);
        for (int i = 0; i < N; i++)
        {
            // printf("t,i %d, %d, sheep value %f\n", t, i, sheep[i][t]);
            ask_sheep(sheep, grass, i, t);
            ask_wolve(wolve, sheep, i, t);
            if (Grass == 1)
                ask_patch(grass, i, t);
        }
    }

    mat2hdf5(sheep, FILE);
    mat2hdf5(wolve, FILEw);
    mat2hdf5int(grass, FILEg);

    // for (int i = 0; i < N; i++)
    // {
    //     free(sheep[i]);
    //     free(wolve[i]);
    // }
    // free(sheep);
    // free(wolve);
    return 0;
}
