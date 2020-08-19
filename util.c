#include "util.h"

bool check_move_valid(int x, int y)
{
    int half = N / 2;
    if ((x >= 0) && (x < half) && (y >= 0) && (y < half))
        return true;
    return false;
}

void move(double animal[N][T], int curX, int curY, int t)
{
    int half = N / 2;
    int randint = (float)rand() / (float)(RAND_MAX)*7;
    if (randint == 0)
    {
        if (check_move_valid(curX - 1, curY))
        {
            animal[curX - 1 + curY * half][t] = animal[curX + curY * half][t] - 1;
            animal[curX + curY * half][t] = 0;
        }
    }
    else if (randint == 1)
    {
        if (check_move_valid(curX - 1, curY))
        {
            animal[curX - 1 + curY * half][t] = animal[curX + curY * half][t] - 1;
            animal[curX + curY * half][t] = 0;
        }
    }
}