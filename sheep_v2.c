#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#define FILE "sheep.h5"
#define FILEw "wolve.h5"
#define FILEg "grass.h5"
#define DATASET "DS1"
#define N 250000
#define T 100
#define initSheepNum 50
#define sheepGainFromFood 4
#define sheepReproduce 4 //%
#define initWolveNum 50
#define wolveGainFromFood 5
#define wolveReproduce 5 //%
#define Grass 1
#define initGrass 100
#define grassRegrowth 10 //time
#define max(a, b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define min(a, b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })
#define error 0.00001
#define debug 0
#define dymatrix

int half = sqrt(N);
int tot_sheep=0;
int tot_wolve = 0;
int tot_grass = 0;

#ifdef dymatrix
int mat2hdf5(double** wdata, char *filename) 
#else
int mat2hdf5(double wdata[N][T], char *filename)
#endif
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
#ifdef dymatrix
int mat2hdf5int(int** wdata, char *filename) 
#else
int mat2hdf5int(int wdata[N][T], char *filename)
#endif
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
#ifdef dymatrix
bool check_sheep_exist(double* animal, int x, int y) 
#else
bool check_sheep_exist(double animal[N], int x, int y) 
#endif
{
    if ((x >= 0) && (x < half) && (y >= 0) && (y < half) && animal[x + y * half] > error)
        return true;
    return false;
}
#ifdef dymatrix
void move(double* animal, int curX, int curY, int *newX, int *newY) 
#else
void move(double animal[N], int curX, int curY, int *newX, int *newY) 
#endif
{
    int flag = false;
    int tmpX = curX, tmpY = curY;
    int count = 0;
    int xlist[8] = {0}, ylist[8] = {0};
    
    // assume they can move randomly within step
    int step = min(half, 50);
    step = step < 4 ? 4 : step;
    int halfstep = step / 2;
    int randx = (float)rand() / (float)(RAND_MAX) * (float)step;
    int randy = (float)rand() / (float)(RAND_MAX) * (float)step;
    tmpX += randx - halfstep;
    tmpX = tmpX < 0 ? (tmpX + half) : tmpX;
    tmpX = tmpX >= half ? (tmpX - half) : tmpX;
    tmpY += randy - halfstep;
    tmpY = tmpY < 0 ? (tmpY + half) : tmpY;
    tmpY = tmpY >= half ? (tmpY - half) : tmpY;

    if (animal[tmpX + tmpY * half] <= error)
    {
        *newX = tmpX;
        *newY = tmpY;
        animal[tmpX + tmpY * half] = animal[curX + curY * half];
        animal[curX + curY * half] = 0;
    }
    else{
        *newX = curX;
        *newY = curY;
    }
    
}
#ifdef dymatrix
void death(double* animal, int curX, int curY, int flag)
#else
void death(double animal[N], int curX, int curY, int flag)
#endif
{
    if (animal[curX + curY * half] <= error)
    {
        animal[curX + curY * half] = 0;
        if (flag == 0){
            tot_sheep -= 1;
            if (debug == 1)
                printf("dead one sheep\n");
        }
        else
        {
            tot_wolve -= 1;
            if (debug == 1)
                printf("dead one wolve\n");
        }
    }
}
#ifdef dymatrix
void create_animal(double* animal, int curX, int curY, int flag)
#else
void create_animal(double animal[N], int curX, int curY, int flag)
#endif
{
    int gainFood = flag == 0 ? sheepGainFromFood : wolveGainFromFood;
    int Energy = max(1.0, 2 * gainFood * (float)rand() / (float)(RAND_MAX));
    int count = 0, cnt = 0;
    int idx = curX + curY *half;
    while ((count < 1) && cnt < N)
    {
        int randint = (float)rand() / (float)(RAND_MAX) * (N - 1);
        if (animal[randint] <= error)
        {
            animal[randint] = max(1.0, Energy);
            count++;
        }
        cnt++;
    }
    if (count>0){
        if (flag == 0){
            if (debug==1)
                printf("create one sheep\n");
            tot_sheep += 1;
        }
        else{
            if (debug == 1)
                printf("create one wolve\n");
            tot_wolve += 1;
        }
    }
    // else{
    //     printf("no place for reproduce flag %d\n",flag);
    // }
}
#ifdef dymatrix
void reproduce(double* animal, int curX, int curY, int flag) 
#else
void reproduce(double animal[N], int curX, int curY, int flag) 
#endif
{
    float randint = (float)rand() / (float)(RAND_MAX)*100;
    int cond_reproduce_rate = (flag == 0) ? sheepReproduce : wolveReproduce;
    if (randint < (float)cond_reproduce_rate){
        animal[curX + curY * half] /= 2;
        create_animal(animal, curX, curY,flag);
    }
}
#ifdef dymatrix
void init_sheep_wolve(double* sheep, double* wolve, int* grass)
#else
void init_sheep_wolve(double sheep[N], double wolve[N], int grass[N])
#endif
{
    int sheepNum = 0;
    while (sheepNum < initSheepNum){
        int randint = (float)rand() / (float)(RAND_MAX) * (N-1);
        if (sheep[randint] < error){
            sheep[randint] = max(2.0, 2 * sheepGainFromFood * (float)rand() / (float)(RAND_MAX));
            // printf("randint %d, sheepNum %d, sheep %f, out %f \n", randint, sheepNum, sheep[randint][0], 2 * sheepGainFromFood * (float)rand() / (float)(RAND_MAX));
            sheepNum++;
            tot_sheep++;
        }
    }
    int wolveNum = 0;
    while (wolveNum < initWolveNum)
    {
        int randint = (float)rand() / (float)(RAND_MAX) * (N-1);
        if (wolve[randint] < error)
        {
            wolve[randint] = max(2.0, 2 * wolveGainFromFood * (float)rand() / (float)(RAND_MAX));
            wolveNum++;
            tot_wolve++;
        }
    }
    if (Grass == 0){
        for (int i =0; i<N; i++)
            grass[i] = 1;
    }else{
        for (int i = 0; i < N; i++)
            grass[i] = (-1) * grassRegrowth;
        int grassNum = 0;
        while (grassNum < initGrass)
        {
            int randint = (float)rand() / (float)(RAND_MAX) * (N - 1);
            if (grass[randint] < error)
            {
                grass[randint] = 1;
                grassNum++;
                tot_grass++;
            }
        }
    }
}
#ifdef dymatrix
void catch_sheep(double* wolve, double* sheep, int *curX, int *curY) 
#else
void catch_sheep(double wolve[N], double sheep[N], int *curX, int *curY) 
#endif
{
    int count = 0;
    int xlist[9] = {0}, ylist[9] = {0};
    int count_list[9] = {0};
    int idx = (*curX) + (*curY) * half;

    // if (sheep[idx] > error)
    // {
    //     sheep[idx] = 0; //kill
    //     tot_sheep--;
    //     wolve[idx] += wolveGainFromFood;
    // }
    
    gen_surroundxy(*curX, *curY, xlist, ylist);
    for (int i = 0; i < 8; i++)
    {
        bool check = check_sheep_exist(sheep, xlist[i], ylist[i]);
        if (check == true)
        {
            count_list[i] = 1;
            count++;
        }
    }
    if (sheep[idx] > error)
    {
        xlist[8] = (*curX);
        ylist[8] = (*curY);
        count_list[8] = 1;
        count++;
    }
    if (count > 0)
    {
        // int randint = (float)rand() / (float)(RAND_MAX) * (count - 1);
        int randint = 8;
        for (int i = 8; i >= 0; i--)
        {
            if (count_list[i] > 0)
            {
                randint = i;
                break;
            }
        }
        sheep[xlist[randint] + half * ylist[randint]] = 0; //kill
        if (debug == 1)
            printf("kill one sheep\n");
        tot_sheep--;
        wolve[idx] += wolveGainFromFood;
    }
}
#ifdef dymatrix
void eatgrass(double* sheep, int* grass, int i)
#else
void eatgrass(double sheep[N], int grass[N], int i)
#endif
{
    // t+1 sheep eat t grass
    if (grass[i]>0){
        sheep[i] += sheepGainFromFood;
        grass[i] = (-1) * grassRegrowth;
        tot_grass--;
    }
}
#ifdef dymatrix
void ask_patch(int* grass, int i)
#else
void ask_patch(int grass[N], int i)
#endif
{
    if (grass[i]==0)
        tot_grass++;
    if (grass[i] < 1)
        grass[i]++;
}
#ifdef dymatrix
void ask_sheep(double *sheep, int *grass, int i)
#else
void ask_sheep(double sheep[N], int grass[N], int i)
#endif
{
    if (sheep[i] < error)
        return;
    int curY = i / half;
    int curX = i - curY * half;
    int newX, newY;
    move(sheep, curX, curY, &newX, &newY);
    if (Grass == 1) 
    {
        int newIdx = newX + newY * half;
        sheep[newIdx]--;
        eatgrass(sheep, grass, newIdx); // t+1 sheep eat t grass
    }
    death(sheep, newX, newY, 0);
    reproduce(sheep, newX, newY, 0);
}
#ifdef dymatrix
void ask_wolve(double *wolve, double *sheep, int i)
#else
void ask_wolve(double wolve[N], double sheep[N], int i) 
#endif
{
    if (wolve[i] < error)
        return;
    int curY = i / half;
    int curX = i - curY * half;
    int newX, newY;
    move(wolve, curX, curY, &newX, &newY);
    wolve[newX + newY*half] -= 1;
    catch_sheep(wolve, sheep, &newX, &newY);
    death(wolve, newX, newY, 1);
    reproduce(wolve, newX, newY, 1);
}
#ifdef dymatrix
void save2mat(double **matTime, double *mat, int t, bool *fmat)
#else
void save2mat(double matTime[N][T], double mat[N], int t, bool* fmat)
#endif
{
    for (int i = 0; i < N; i++){
        matTime[i][t] = mat[i];
        if (mat[i] > error)
            fmat[i] = true;
        else
            fmat[i] = false;
    }
}
#ifdef dymatrix
void save2matInt(int** matTime, int* mat, int t)
#else
void save2matInt(int matTime[N][T], int mat[N], int t)
#endif
{
    for (int i = 0; i < N; i++)
    {
        matTime[i][t] = mat[i];
    }
}

int main(void)
{
    clock_t start_t, end_t, total_t;
    start_t = clock();
    assert(initSheepNum <= N && "sheep init number should be smaller than N");
    assert(initWolveNum <= N && "wolve init number should be smaller than N");
    assert(initGrass <= N && "grass init number should be smaller than N");

#ifdef dymatrix
    double **sheep = malloc(N * sizeof(double *));
    double **wolve = malloc(N * sizeof(double *));
    int **grass = malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++){
        sheep[i] = malloc(T * sizeof(double));
        wolve[i] = malloc(T * sizeof(double));
        grass[i] = malloc(T * sizeof(int));
    }
    double* curSheep = malloc(N * sizeof(double));
    double* curWolve = malloc(N * sizeof(double));
    int* curGrass = malloc(N * sizeof(int));

#else
    double sheep[N][T] = {0};
    double wolve[N][T] = {0};
    int grass[N][T] = {0};

    double curSheep[N] = {0};
    double curWolve[N] = {0};
    int curGrass[N] = {0};
#endif
    bool *fsheep = malloc(N * sizeof(bool));
    bool *fwolve = malloc(N * sizeof(bool));

    init_sheep_wolve(curSheep, curWolve, curGrass);
    // testing for all 1 case
    // for (int i = 0; i < N; i++)
    // {
    //     curSheep[i] = 1;
    //     curWolve[i] = 1;
    //     curGrass[i] = 1;
    // }
    save2mat(sheep, curSheep, 0, fsheep);
    save2mat(wolve, curWolve, 0, fwolve);
    save2matInt(grass, curGrass, 0);

    for (int t = 1; t < T; t++)
    {
        int curSheep_tot = tot_sheep;
        int curWolve_tot = tot_wolve;
        int curGrass_tot = tot_grass;

        for (int i = 0; i < N; i++)
        {
            if (fsheep[i] == true)
                ask_sheep(curSheep, curGrass, i);
            if (fwolve[i] == true)
                ask_wolve(curWolve, curSheep, i);
            if (Grass == 1)
                ask_patch(curGrass, i);
        }
        int acc_s = 0;
        int acc_w = 0;
        int acc_g = 0;
        for (int i = 0; i < N; i++)
        {
            if (curWolve[i] > 0){
                acc_w++;
                // printf("wolf %f, ", curWolve[i]);
            }
                
            if (curSheep[i] > 0)
                acc_s++;
            if (curGrass[i] > 0)
                acc_g++;
        }
        // printf("\n");
        if (debug == 1){
            printf("t %d\n", t);
            printf("s %d,%d,%d\n", curSheep_tot, tot_sheep, acc_s);
            printf("w %d,%d,%d\n", curWolve_tot, tot_wolve, acc_w);
            printf("g %d,%d,%d\n", curGrass_tot, tot_grass, acc_g);
            printf("\n");
        }
        // else
        //     printf("t %d, s %d, w %d, g %d\n", t, tot_sheep, tot_wolve, tot_grass);

        save2mat(sheep, curSheep, t, fsheep);
        save2mat(wolve, curWolve, t, fwolve);
        save2matInt(grass, curGrass, t);
    }

    mat2hdf5(sheep, FILE);
    mat2hdf5(wolve, FILEw);
    mat2hdf5int(grass, FILEg);

    end_t = clock();
    total_t = (double)(end_t - start_t) / CLOCKS_PER_SEC;
    printf("End of program total_t = %f\n", total_t);

#ifdef dymatrix
    for (int i = 0; i < N; i++)
    {
        free(sheep[i]);
        free(wolve[i]);
        free(grass[i]);
    }
    free(sheep);
    free(grass);
    free(wolve);
#endif
    free(fsheep);
    free(fwolve);

    return 0;
}
