#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#define FILE "animal.h5"
#define N 2500
#define T 1000
#define initSheepNum 100
#define sheepGainFromFood 4
#define sheepReproduce 4 //%
#define initWolveNum 50
#define wolveGainFromFood 20
#define wolveReproduce 5 //%
#define Grass 1
#define initGrass 1200
#define grassRegrowth 30 //time
#define max(a, b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define min(a, b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })
#define error 0.00001
#define debug 0

int half = sqrt(N);
int tot_sheep=0;
int tot_wolve = 0;
int tot_grass = 0;


int mat2hdf5(double *sdata, double *wdata, int *gdata, int *count, int *setting, char *filename)
{
    hid_t file, space, space_c, space_set, dsets, dsetw, dsetg, dsetc, dsetset; /* Handles */
    herr_t status;
    hsize_t dimset[1] = {2};
    hsize_t dims[1] = {N * T};
    hsize_t dims_count[1] = {3 * T};

    char DATASETs[20] = "Ds_sheep";
    char DATASETw[20] = "Ds_wolve";
    char DATASETg[20] = "Ds_grass";
    char DATASETc[20] = "Ds_count";
    char DATASETset[20] = "Ds_set";

    file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    space = H5Screate_simple(1, dims, NULL);
    space_c = H5Screate_simple(1, dims_count, NULL);
    space_set = H5Screate_simple(1, dimset, NULL);
    dsets = H5Dcreate(file, DATASETs, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dsetw = H5Dcreate(file, DATASETw, H5T_IEEE_F64LE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dsetg = H5Dcreate(file, DATASETg, H5T_STD_I64BE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dsetc = H5Dcreate(file, DATASETc, H5T_STD_I64BE, space_c, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dsetset = H5Dcreate(file, DATASETset, H5T_STD_I64BE, space_set, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dsets, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, sdata);
    status = H5Dwrite(dsetw, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, wdata);
    status = H5Dwrite(dsetg, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, gdata);
    status = H5Dwrite(dsetc, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, count);
    status = H5Dwrite(dsetset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, setting);

    status = H5Dclose(dsets);
    status = H5Dclose(dsetw);
    status = H5Dclose(dsetg);
    status = H5Dclose(dsetset);
    status = H5Sclose(space);
    status = H5Sclose(space_c);
    status = H5Sclose(space_set);
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
bool check_sheep_exist(double* animal, int x, int y) 
{
    if ((x >= 0) && (x < half) && (y >= 0) && (y < half) && animal[x + y * half] > error)
        return true;
    return false;
}
void valify_move(int* tmpX, int* tmpY)
{
    *tmpX = *tmpX < 0 ? (*tmpX + half) : *tmpX;
    *tmpX = *tmpX >= half ? (*tmpX - half) : *tmpX;
    *tmpY = *tmpY < 0 ? (*tmpY + half) : *tmpY;
    *tmpY = *tmpY >= half ? (*tmpY - half) : *tmpY;
}
void move(double *animal, int *animaldir, int curX, int curY, int *newX, int *newY)
{
    int xlist[8] = {0}, ylist[8] = {0};
    
    // assume they can move randomly in 8 direction
    gen_surroundxy(curX, curY, xlist, ylist);
    for (int i = 0; i < 8; i++)
        valify_move(&xlist[i], &ylist[i]);
    int count = 0, cnt = 0;
    while ((count < 1) && (cnt < 9))
    {
        // int randdir = (float)rand() / (float)(RAND_MAX)*7.99;
        int randdir = (float)rand() / (float)(RAND_MAX)*2.99; //Brownian movement
        // int randdir = (float)rand() / (float)(RAND_MAX)*2.99; //Brownian movement
        randdir += animaldir[curX + curY * half] - 1;
        randdir = randdir < 0 ? (randdir + 8) : randdir;
        randdir = randdir > 7 ? (randdir - 8) : randdir;
        int randidx = xlist[randdir] + ylist[randdir] * half;
        // printf("%d, %d, %d\n", animaldir[curX + curY * half], animaldir[randidx], randdir);

        if (animal[randidx] <= error)
        {
            int curidx = curX + curY * half;
            *newX = xlist[randdir];
            *newY = ylist[randdir];
            animal[randidx] = animal[curidx];
            animaldir[randidx] = animaldir[curidx];
            animal[curidx] = 0;
            animaldir[curidx] = -1;
            count++;
        }
        cnt++;
    }
    if (count < 1){
        printf("all occupied\n");
        *newX = curX;
        *newY = curY;
    }
        
}
void move50(double *animal, int *animaldir, int curX, int curY, int *newX, int *newY)
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
    tmpY += randy - halfstep;
    valify_move(&tmpX, &tmpY);

    if (animal[tmpX + tmpY * half] <= error)
    {
        *newX = tmpX;
        *newY = tmpY;
        animal[tmpX + tmpY * half] = animal[curX + curY * half];
        animal[curX + curY * half] = 0;
    }
    else
    {
        *newX = curX;
        *newY = curY;
    }
}
void death(double *animal, int *animaldir, int curX, int curY, int flag)
{
    if (animal[curX + curY * half] <= error)
    {
        animal[curX + curY * half] = 0;
        animaldir[curX + curY * half] = -1;
        if (flag == 0)
        {
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
void create_animal(double *animal, int* animaldir, int curX, int curY, int flag)
{
    int gainFood = flag == 0 ? sheepGainFromFood : wolveGainFromFood;
    int Energy = max(1.0, 2 * gainFood * (float)rand() / (float)(RAND_MAX));
    int count = 0, cnt = 0;
    int idx = curX + curY *half;
    int xlist[8] = {0}, ylist[8] = {0};

    // assume they can move randomly in 8 direction
    gen_surroundxy(curX, curY, xlist, ylist);
    for (int i = 0; i < 8; i++){
        valify_move(&xlist[i], &ylist[i]);
        int randint = xlist[i] + ylist[i] * half;
        if (animal[randint] <= error)
        {
            animal[randint] = max(1.0, Energy);
            animaldir[randint] = (float)rand() / (float)(RAND_MAX)*7.99;
            count++;
            break;
        }
    }
    while ((count < 1) && cnt < N)
    {
        int randint = (float)rand() / (float)(RAND_MAX) * (N - 1);
        if (animal[randint] <= error)
        {
            animal[randint] = max(1.0, Energy);
            animaldir[randint] =  (float)rand() / (float)(RAND_MAX) * 7.99;
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
    else{
        printf("no place for reproduce flag %d\n",flag);
    }
}
void reproduce(double* animal, int* animaldir, int curX, int curY, int flag) 
{
    float randint = (float)rand() / (float)(RAND_MAX)*100;
    int cond_reproduce_rate = (flag == 0) ? sheepReproduce : wolveReproduce;
    if (randint < (float)cond_reproduce_rate){
        animal[curX + curY * half] /= 2;
        create_animal(animal, animaldir, curX, curY, flag);
    }
}
void init_sheep_wolve(double* sheep, int* sheepdir, double* wolve, int* wolfdir, int* grass)
{
    int sheepNum = 0;
    while (sheepNum < initSheepNum){
        int randint = (float)rand() / (float)(RAND_MAX) * (N-1);
        if (sheep[randint] < error){
            sheep[randint] = max(2.0, 2 * sheepGainFromFood * (float)rand() / (float)(RAND_MAX));
            // printf("randint %d, sheepNum %d, sheep %f, out %f \n", randint, sheepNum, sheep[randint][0], 2 * sheepGainFromFood * (float)rand() / (float)(RAND_MAX));
            sheepdir[randint] = 7.99 * (float)rand() / (float)(RAND_MAX);
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
            wolve[randint] = max(1.0, 2 * wolveGainFromFood * (float)rand() / (float)(RAND_MAX));
            wolfdir[randint] = 7.99 * (float)rand() / (float)(RAND_MAX);
            wolveNum++;
            tot_wolve++;
        }
    }
    if (Grass == 0){
        for (int i =0; i<N; i++)
            grass[i] = 1;
        tot_grass = N;
    }else{
        for (int i = 0; i < N; i++)
            grass[i] = (-1) * max(1.0, grassRegrowth * (float)rand() / (float)(RAND_MAX) * (N - 1));
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
void catch_sheep(double* wolve, double* sheep, int* sheepdir, int *curX, int *curY) 
{
    int count = 0;
    int xlist[9] = {0}, ylist[9] = {0};
    int count_list[9] = {0};
    int idx = (*curX) + (*curY) * half;

    if (sheep[idx] > error)
    {
        sheep[idx] = 0; //kill
        sheepdir[idx] = -1;
        tot_sheep--;
        wolve[idx] += wolveGainFromFood;
        if (debug == 1)
            printf("kill one sheep\n");
    }
    
    // gen_surroundxy(*curX, *curY, xlist, ylist);
    // for (int i = 0; i < 8; i++)
    // {
    //     bool check = check_sheep_exist(sheep, xlist[i], ylist[i]);
    //     if (check == true)
    //     {
    //         count_list[i] = 1;
    //         count++;
    //     }
    // }
    // if (sheep[idx] > error)
    // {
    //     xlist[8] = (*curX);
    //     ylist[8] = (*curY);
    //     count_list[8] = 1;
    //     count++;
    // }
    // if (count > 0)
    // {
    //     // int randint = (float)rand() / (float)(RAND_MAX) * (count - 1);
    //     int randint = 8;
    //     for (int i = 8; i >= 0; i--)
    //     {
    //         if (count_list[i] > 0)
    //         {
    //             randint = i;
    //             break;
    //         }
    //     }
    //     sheep[xlist[randint] + half * ylist[randint]] = 0; //kill
    //     if (debug == 1)
    //         printf("kill one sheep\n");
    //     tot_sheep--;
    //     wolve[idx] += wolveGainFromFood;
    // }
}
void eatgrass(double* sheep, int* grass, int i)
{
    // t+1 sheep eat t grass
    if (grass[i]>0){
        sheep[i] += sheepGainFromFood;
        grass[i] = (-1) * grassRegrowth;
        tot_grass--;
    }
}
void ask_patch(int* grass, int i)
{
    if (grass[i]==0)
        tot_grass++;
    if (grass[i] < 1)
        grass[i]++;
}
void ask_sheep(double *sheep, int* sheepdir, int *grass, int i)
{
    if (sheep[i] < error)
        return;
    int curY = i / half;
    int curX = i - curY * half;
    int newX, newY;
    move(sheep, sheepdir, curX, curY, &newX, &newY);
    if (Grass == 1) 
    {
        int newIdx = newX + newY * half;
        sheep[newIdx]--;
        eatgrass(sheep, grass, newIdx); // t+1 sheep eat t grass
    }
    death(sheep, sheepdir, newX, newY, 0);
    reproduce(sheep, sheepdir, newX, newY, 0);
}
void ask_wolve(double *wolve, int *wolfdir, double *sheep, int *sheepdir, int i)
{
    if (wolve[i] < error)
        return;
    int curY = i / half;
    int curX = i - curY * half;
    int newX, newY;
    move(wolve, wolfdir, curX, curY, &newX, &newY);
    wolve[newX + newY*half] -= 1;
    catch_sheep(wolve, sheep, sheepdir, & newX, &newY);
    death(wolve, wolfdir, newX, newY, 1);
    reproduce(wolve, wolfdir, newX, newY, 1);
}
void save2mat(double *matTime, double *mat, int t, bool *fmat)
{
    int base = N * t;
    for (int i = 0; i < N; i++){
        matTime[i + base] = mat[i];
        if (mat[i] > error)
            fmat[i] = true;
        else
            fmat[i] = false;
    }
}
void save2matInt(int* matTime, int* mat, int t)
{
    int base = N * t;
    for (int i = 0; i < N; i++)
        matTime[i + base] = mat[i];
}

int main(void)
{
    srand(time(NULL));
    clock_t start_t, end_t, total_t;
    start_t = clock();
    assert(half * half == N && "N should be able to square root");
    assert(initSheepNum <= N && "sheep init number should be smaller than N");
    assert(initWolveNum <= N && "wolve init number should be smaller than N");
    assert(initGrass <= N && "grass init number should be smaller than N");

    double *sheep = malloc(N * T * sizeof(double));
    double *wolve = malloc(N * T * sizeof(double));
    int *grass = malloc(N * T * sizeof(int));
    double *curSheep = malloc(N * sizeof(double));
    double *curWolve = malloc(N * sizeof(double));
    int *curGrass = malloc(N * sizeof(int));
    int *curSheepDir = malloc(N * sizeof(int));
    int *curWolveDir = malloc(N * sizeof(int));
    int *animalNum = malloc(3 * T * sizeof(int));
    int *setting = malloc(2 * sizeof(int));
    for (int i = 0 ; i < 3*T; i++)
        animalNum[i] = 0;
    setting[0] = N;
    setting[1] = T;

    for (int i = 0 ; i < N; i++){
        curSheepDir[i] = -1;
        curWolveDir[i] = -1;
    }

    bool *fsheep = malloc(N * sizeof(bool));
    bool *fwolve = malloc(N * sizeof(bool));

    init_sheep_wolve(curSheep, curSheepDir, curWolve, curWolveDir, curGrass);
    animalNum[0] = tot_sheep;
    animalNum[1] = tot_wolve;
    animalNum[2] = tot_grass;
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
            if (fsheep[i] == true){
                ask_sheep(curSheep, curSheepDir, curGrass, i);
            }
            if (fwolve[i] == true)
                ask_wolve(curWolve, curWolveDir, curSheep, curSheepDir, i);
            if (Grass == 1)
                ask_patch(curGrass, i);
        }
        int acc_s = 0;
        int acc_w = 0;
        int acc_g = 0;
        for (int i = 0; i < N; i++)
        {
            if (curWolve[i] > 0)
            {
                acc_w++;
                // printf("wolf %f, ", curWolve[i]);
            }

            if (curSheep[i] > 0)
                acc_s++;
            if (curGrass[i] > 0)
                acc_g++;
        }
        animalNum[0 + t * 3] = acc_s;
        animalNum[1 + t * 3] = acc_w;
        animalNum[2 + t * 3] = acc_g;
        if (debug == 1)
        {

            // printf("\n");
            printf("t %d\n", t);
            printf("s %d,%d,%d\n", curSheep_tot, tot_sheep, acc_s);
            printf("w %d,%d,%d\n", curWolve_tot, tot_wolve, acc_w);
            printf("g %d,%d,%d\n", curGrass_tot, tot_grass, acc_g);
            printf("\n");
        }
        else
            printf("t %d, s %d, w %d, g %d\n", t, tot_sheep, tot_wolve, tot_grass);
        save2mat(sheep, curSheep, t, fsheep);
        save2mat(wolve, curWolve, t, fwolve);
        save2matInt(grass, curGrass, t);
    }

    mat2hdf5(sheep, wolve, grass, animalNum,setting, FILE);
    end_t = clock();
    total_t = (double)(end_t - start_t) / (double)CLOCKS_PER_SEC;
    printf("End of program total_t = %u, %u, %u, %u\n", total_t, end_t, start_t, CLOCKS_PER_SEC);

    free(sheep);
    free(grass);
    free(wolve);
    free(curSheep);
    free(curWolve);
    free(curGrass);
    free(fsheep);
    free(fwolve);
    free(animalNum);
    free(setting);

    return 0;
    }
