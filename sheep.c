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

const int half = N/2;
int tot_sheep = initSheepNum;
int tot_wolve = initWolveNum;
void gen_surroundxy(int curX, int curY, int xlist[8], int ylist[8]){
    int xl[8] = {-1, -1 ,0 ,1, 1,1,0,-1};
    int yl[8] = {0  ,-1,-1,-1, 0,1,1, 1};
    for (int i = 0; i < 8; i++){
        xlist[i] = curX + xl[i];
        ylist[i] = curY + yl[i];
    }
}

bool check_move_valid(double animal[N][T], int curX, int curY, int x, int y, int t)
{
    if ((x >= 0) && (x < half) && (y >= 0) && (y < half) && animal[x + y * half][t]==0)
        animal[x + y * half][t] = animal[curX + curY * half][t];
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

void move(double animal[N][T], int curX, int curY, int t, int* newX, int* newY)
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
        flag = check_move_valid(animal, curX, curY, tmpX, tmpY, t);
        if (flag == true){
            *newX = tmpX;
            *newY = tmpY;
            break;
        }
    }
    if (flag == false)
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
    int tmpEnergy = animal[curIdx][t];
    int newX, newY;

    animal[curIdx][t] = Energy;
    move(animal, curX, curY, t, &newX, &newY);
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
        printf("%d\n",randint);
        if (sheep[randint][0] != 0)
            sheep[randint][0] = 2 * sheepGainFromFood * (float)rand() / (float)(RAND_MAX);
        sheepNum++;
    }
    int wolveNum = 0;
    while (wolveNum < initWolveNum)
    {
        int randint = (float)rand() / (float)(RAND_MAX) * (N - 1);
        if (wolve[randint][0] != 0)
            wolve[randint][0] = 2 * wolveGainFromFood * (float)rand() / (float)(RAND_MAX);
        wolveNum++;
    }
}
void catch_sheep(double wolve[N][T], double sheep[N][T], int *curX, int *curY, int t)
{
    int count = 0;
    int xlist[8] = {0}, ylist[8] = {0};
    int count_list[8] = {0};
    int idx = (*curX) + (*curY) * half;
    gen_surroundxy(*curX, *curY, xlist, ylist);
    for (int i = 0; i < 8; i++)
    {
        bool check = check_sheep_exist(sheep, xlist[i], ylist[i], t);
        if (check == true){
            count_list[i] = 1;
            count+=1;
        }
    }
    if (count > 0){
        int randint = (float)rand() / (float)(RAND_MAX) * (count - 1);
        sheep[xlist[randint] + half * ylist[randint]][t] = 0;
        wolve[idx][t] += wolveGainFromFood;
    }
}

void ask_sheep(double sheep[N][T], int i, int t){
    int curY = i/half;
    int curX = i - curY * half;
    int newX, newY;
    move(sheep, curX, curY, t, &newX, &newY);
    // if Grass == 1{
    //     sheep[newX + newY * half][t] = sheep[newX + newY * half][t]-1
    // }
    death(sheep, newX, newY, t, 0);
    reproduce(sheep, newX, newY, t, 0);
}
void ask_wolve(double wolve[N][T], double sheep[N][T], int i, int t)
{
    int curY = i / half;
    int curX = i - curY * half;
    int newX, newY;
    move(wolve, curX, curY, t, &newX, &newY);
    wolve[i][t] -= 1;
    catch_sheep(wolve, sheep, &newX, &newY, t);
    death(wolve, newX, newY, t, 1);
    reproduce(wolve, newX, newY, t, 0);
}

int main(void)
{
    double sheep[N][T] = {0}, wolve[N][T] = {0};
    init_sheep_wolve(sheep, wolve);
    for (int t = 0; t < T; t++)
    {
        printf("%d", tot_sheep);
        printf("%d", tot_wolve);
        for (int i = 0; i < N; i++)
        {
                ask_sheep(sheep, i, t);
                ask_wolve(wolve, sheep, i, t);
                // if (Grass == 1)
                //     ask_patch(grass,i,t)
                
        }
    }

    // mat2hdf5(time_series);
    return 0;
}
