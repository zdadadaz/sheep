// #include "hdf5.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <math.h>
// #define FILE "animal.h5"
#define N 25
#define T 100
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
#include <vector>
#include <unordered_set>
#include "hdf5.h"

using namespace std;
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

void gen_surroundxy(int curX, int curY, int xlist[8], int ylist[8])
{
    int xl[8] = {-1, -1, 0, 1, 1, 1, 0, -1};
    int yl[8] = {0, -1, -1, -1, 0, 1, 1, 1};
    for (int i = 0; i < 8; i++)
    {
        xlist[i] = curX + xl[i];
        ylist[i] = curY + yl[i];
    }
}
void valify_move(int *tmpX, int *tmpY)
{
    *tmpX = *tmpX < 0 ? (*tmpX + half) : *tmpX;
    *tmpX = *tmpX >= half ? (*tmpX - half) : *tmpX;
    *tmpY = *tmpY < 0 ? (*tmpY + half) : *tmpY;
    *tmpY = *tmpY >= half ? (*tmpY - half) : *tmpY;
}
class Grassclass{
	int number;
	int sx,sy;
	public:
		Grassclass(int x, int y, int count){sx = x; sy = y; number=count;};
		void countdown();
        int x() { return sx; };
        int y() { return sy; };
        int gNum() { return number; };
        void sNum(int count) { number = count;};
        void reset(){number = (-1) * grassRegrowth;};
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
        void move();
        int x(){return sx;};
        int y(){return sy;};
        int d(){return sd;};
};
RandomWalk::RandomWalk(int inx, int iny, int ind){sx = inx; sy=iny; sd = ind;}

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
		void reduceEnergy(){energy -= 1;};
        void addEnergy() { energy += gainFood; }
        void divideEnergy() { energy /= 2;}
        float gEnergy() { return energy; };
        int gFlag() {return sflag;};
        //void die(list<Animal>& animallist);
		//void reproduce();
};
void init_sheep_wolve(list<Animal> &sheeps, int flag)
{
    int sheepNum = 0;
    int direction,x,y;
    int initnumber;
    int* tmpTotl;

    if (flag == 0)
    {
        tmpTotl = &tot_sheep;
        initnumber = initSheepNum;
    }
    else
    {
        tmpTotl = &tot_wolve;
        initnumber = initWolveNum;
    }
    
    while (sheepNum < initnumber)
    {
        int randint = (float)rand() / (float)(RAND_MAX) * (N - 1);
        direction = 7.99 * (float)rand() / (float)(RAND_MAX);
        y = randint / half;
        x = randint - y * half;
        Animal animal(x, y, direction, flag);
        sheeps.push_back(animal);
        sheepNum++;
        (*tmpTotl)++;
    }
    
}
void init_grass(vector<Grassclass> &grasses)
{
    int grassNum = 0;
    int count = 0;
    int x, y;
    for (int i = 0; i < N; i++){
        y = i / half;
        x = i - y * half;
        count = (-1) * max(1.0, grassRegrowth * (float)rand() / (float)(RAND_MAX));
        Grassclass grass(x, y, count);
        grasses.push_back(grass);
    }
    std::unordered_set<int> grassSet;
    while (grassNum < initGrass)
    {
        int randint = (float)rand() / (float)(RAND_MAX) * (N - 1);
        if (grassSet.find(randint) != grassSet.end())
            continue;
        grassSet.insert(randint);
        grasses[randint].sNum(1);
        grassNum++;
        tot_grass++;
    }

    // for (auto i = g1.begin(); i != g1.end(); ++i)
}
void eatGrass(Animal &sheep, vector<Grassclass> &grasses)
{
    if (grasses[sheep.x() + half * sheep.y()].gNum() == 1)
    {
        sheep.addEnergy();
        grasses[sheep.x() + half * sheep.y()].reset();
        tot_grass--;
    }
}
void eatSheep(Animal &wolf, list<Animal> &sheeplist)
{
    for (std::list<Animal>::iterator it = sheeplist.begin(); it != sheeplist.end(); ++it)
    {
        if ((wolf.x() == (*it).x()) && (wolf.y() == (*it).y()) && ((*it).gEnergy()))
        {
            wolf.addEnergy();
            sheeplist.erase(it);
            tot_sheep--;
        }
    }
}
void create_animal(list<Animal> &animals, Animal &animal)
{
    int xlist[8] = {0}, ylist[8] = {0};
    // assume they can move randomly in 8 direction
    gen_surroundxy(animal.x(), animal.y(), xlist, ylist);
    int i = (float)rand() / (float)(RAND_MAX)*7.99; //choose one direction to move
    valify_move(&xlist[i], &ylist[i]);
    int d = (float) rand() / (float)(RAND_MAX)*7.99; // initialize direction
    Animal newAnimal(xlist[i], ylist[i], d, animal.gFlag());
    animals.push_back(newAnimal);
}
void reproduce(list<Animal> &animals, Animal &animal)
{
    float randint = (float)rand() / (float)(RAND_MAX)*100;
    int cond_reproduce_rate = (animal.gFlag() == 0) ? sheepReproduce : wolveReproduce;
    if (randint < (float)cond_reproduce_rate)
    {
        animal.divideEnergy();
        create_animal(animals, animal);
    }
}
void death(list<Animal> &animals, std::list<Animal>::iterator it)
{
    if ((*it).gEnergy() <= error)
    {
        animals.erase(it);
        if ((*it).gFlag() == 0)
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
void ask_sheep(list<Animal>& sheeplist, vector<Grassclass>& grasslist)
{
    for (std::list<Animal>::iterator it = sheeplist.begin(); it != sheeplist.end(); ++it)
    {
        (*it).move();
        if (Grass == 1)
        {
            (*it).reduceEnergy();
            eatGrass((*it), grasslist);
        }
        if ((*it).gEnergy() > error)
        {
            reproduce(sheeplist, (*it));
        }
        death(sheeplist, it);
    }
}
void ask_wolf(list<Animal> &wolflist, list<Animal> &sheeplist)
{
    for (std::list<Animal>::iterator it = wolflist.begin(); it != wolflist.end(); ++it)
    {
        (*it).move();
        if (Grass == 1)
        {
            (*it).reduceEnergy();
            eatSheep((*it), sheeplist);
        }
        if ((*it).gEnergy() > 0)
        {
            reproduce(wolflist, (*it));
        }
        death(wolflist, it);
    }
}
void ask_patch(vector<Grassclass> &grasslist)
{
    for (Grassclass grass : grasslist){
        if (grass.gNum() == 0)
            tot_grass++;
        if (grass.gNum() < 1)
            grass.sNum(grass.gNum()+1);
    }
}
void save2mat(double *matTime, list<Animal> *mat, int t)
{
    int base = N * t;
    for (Animal animal : mat){
        matTime[animal.x() + half * animal.y() + base] = animal.gEnerge();
    }
}
void save2matInt(int *matTime, vector<Grassclass> &grasslist, int t)
{
    int base = N * t;
    for (int i = 0; i < N; i++)
        matTime[i + base] = grasslist[i].gNum();
}

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
    std::vector<Grassclass> grasslist;
    double *sheep = new double (N * T * sizeof(double));
    double *wolve = new double (N * T * sizeof(double));
    int *grass = new int (N * T * sizeof(int));
    int *animalNum = new int (3 * T * sizeof(int));
    for (int i = 0; i < 3 * T; i++)
        animalNum[i] = 0;
    int *setting = new int (2 * sizeof(int));
    animalNum[0] = tot_sheep;
    animalNum[1] = tot_wolve;
    animalNum[2] = tot_grass;
    setting[0] = N;
    setting[1] = T;
    for (int t = 1; t < T; t++){
        //init
        init_sheep_wolve(sheeplist, 0);
        init_sheep_wolve(wolflist, 1);
        init_grass(grasslist);
        // printf("%d, %d, %d\n", tot_sheep, tot_wolve, tot_grass);
        //ask sheep
        ask_sheep(sheeplist, grasslist);
        //ask wolf
        ask_wolf(wolflist, sheeplist);
        //ask grass
        ask_patch(grasslist);
        save2mat(sheep, sheeplist, t);
        save2mat(wolve, wolflist, t);
        save2matInt(grass, grasslist, t);
        animalNum[0 + t * 3] = tot_sheep;
        animalNum[1 + t * 3] = tot_wolve;
        animalNum[2 + t * 3] = tot_grass;
    }
    mat2hdf5(sheep, wolve, grass, animalNum, setting, FILE);
    free(sheep);
    free(grass);
    free(wolve);
    free(animalNum);
    free(setting);
    return 0;
}
