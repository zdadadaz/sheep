#include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#define FILE "animal.h5"
#define N 25
#define T 1000
#define initSheepNum 25
#define sheepGainFromFood 4
#define sheepReproduce 4 //%
#define initWolveNum 25
#define wolveGainFromFood 20
#define wolveReproduce 5 //%
#define Grass 1
#define initGrass 25
#define grassRegrowth 30 //time
#define max(a, b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define min(a, b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })
#define error 0.00001
#define debug 0
#include <list>
using namespace std;
int half = sqrt(N);
int tot_sheep=0;
int tot_wolve = 0;
int tot_grass = 0;

class Grassclass{
	int number;
	int sx,sy;
	public:
		Grassclass(int x, int y, int count){sx = x; sy = y; number=count;};
		void countdown();

};
void Grassclass::countdown(){
    if (number==0)
        tot_grass++;
    if (number < 1)
        number++;
}

class RandomWalk{
    int sx,sy,sd;
    public:
        RandomWalk(int x, int y, int d);
        void valify_move(int *tmpX, int* tmpY);
        void gen_surroundxy(int curX, int curY, int xlist[8], int ylist[8]);
        void move();
        int x(){return sx;};
        int y(){return sy;};
        int d(){return sd;};
};
RandomWalk::RandomWalk(int inx, int iny, int ind){sx = inx; sy=iny; sd = ind;}
void RandomWalk::gen_surroundxy(int curX, int curY, int xlist[8], int ylist[8]){
    int xl[8] = {-1, -1 ,0 ,1, 1,1,0,-1};
    int yl[8] = {0  ,-1,-1,-1, 0,1,1, 1};
    for (int i = 0; i < 8; i++){
        xlist[i] = curX + xl[i];
        ylist[i] = curY + yl[i];
    }
}
void RandomWalk::valify_move(int* tmpX, int* tmpY){
    *tmpX = *tmpX < 0 ? (*tmpX + half) : *tmpX;
    *tmpX = *tmpX >= half ? (*tmpX - half) : *tmpX;
    *tmpY = *tmpY < 0 ? (*tmpY + half) : *tmpY;
    *tmpY = *tmpY >= half ? (*tmpY - half) : *tmpY;
}
void RandomWalk::move(){
    int curX = sx;
    int curY = sy;

    int xlist[8] = {0}, ylist[8] = {0};
    // assume they can move randomly in 8 direction
    gen_surroundxy(curX, curY, xlist, ylist);
    for (int i = 0; i < 8; i++)
        valify_move(&xlist[i], &ylist[i]);
    int randdir = (float)rand() / (float)(RAND_MAX)*2.99; //Brownian movement
    randdir += sd - 1;
    randdir = randdir < 0 ? (randdir + 8) : randdir;
    randdir = randdir > 7 ? (randdir - 8) : randdir;
    sx = xlist[randdir];
    sy = ylist[randdir];
    sd = randdir;
}

class Animal: public RandomWalk{
	int sflag;
	float energy;
	float gainFood;
	public:
		Animal(int x, int y, int d, int flag): RandomWalk(x,y,d){
		sflag = flag;
     	gainFood = flag == 0 ? sheepGainFromFood : wolveGainFromFood;
     	double energy_rand = max(1.0, (double)2 * gainFood * (double)rand() / (double)(RAND_MAX));
     	energy = max(1.0, energy_rand);
		};
		//void reduceEnergy(){energy -= 1;};
		//void eat();
		//void die(list<Animal>& animallist);
		//void reproduce();
};
//class linkedlist{
//	private:
//		Animal *head, *tail;
//	public:
//		
//
//};

int main(void)
{
    srand(time(NULL));
    assert(N>=0 && half * half == N && "N should be able to square root and greater equals to zero");
    assert(initSheepNum>=0 && initSheepNum <= N && "sheep init number should be smaller than N and greater equals to zero");
    assert(initWolveNum>=0 && initWolveNum <= N && "wolve init number should be smaller than N and greater equals to zero");
    assert(initGrass>=0 && initGrass <= N && "grass init number should be smaller than N and greater equals to zero");

    //double *sheep = new double(N * T * sizeof(double));
	std::list<Animal> sheeplist;
	std::list<Animal> wolflist;
	//Animal sheep (2,3,0,0);
	//printf("%d, %d, %d\n",sheep.x(),sheep.y(),sheep.d());
	//sheep.move();
	//printf("%d, %d, %d\n",sheep.x(),sheep.y(),sheep.d());

	//init
	//
	//ask sheep
	//ask wolf
	//ask grass

	
return 0;
}
