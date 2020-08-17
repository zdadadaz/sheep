#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>

#define FILE "sheep.h5"
#define DATASET "DS1"
#define N 100
#define T 100
#define initSheepNum 10
#define sheepGainFromFood 4
#define sheepReproduce 4 //%
#define initWolveNum 5
#define wolveGainFromFood 20
#define sheepReproduce 5 //%
#define initGrass 1


void init_sheep_wolve(double** sheep, double** wolve)
{
    int sheepNum = 0;
    while (sheepNum < initSheepNum){
        int randint = (float)rand() / (float)(RAND_MAX)* N;
        if (sheep[randint][0] != 0)
            sheep[randint][0] = 2 * sheepGainFromFood * (float)rand() / (float)(RAND_MAX);
        sheepNum++;
    }
    int wolveNum = 0;
    while (wolveNum < initWolveNum)
    {
        int randint = (float)rand() / (float)(RAND_MAX)*N;
        if (wolve[randint][0] != 0)
            wolve[randint][0] = 2 * wolveGainFromFood * (float)rand() / (float)(RAND_MAX);
        wolveNum++;
    }
}

int main(void)
{
    double sheep[N][T] = {0}, wolve[N][T] = {0};
    float random_float;
    
    for (int i = 0; i < N; i++)
    {
        for (int t = 0; t < T; t++)
        {
            if (t == 0)
            {
                time_series[i][t] = 1; // initial radiation
            }
            else
            {
                random_float = (float)rand() / (float)(RAND_MAX); //random float from 0 to 1
                if (random_float < DECAY_CHANCE)
                { // less than DECAY_CHANCE, particle decay
                    time_series[i][t] = decay_rate * time_series[i][t - 1];
                }
                else
                { // no decay
                    time_series[i][t] = time_series[i][t - 1];
                }
            }
        }
    }

    // mat2hdf5(time_series);
    return 0;
}
