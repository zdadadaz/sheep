#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <math.h>
// #define N 10000
// #define T 1000
// #define initSheepNum 400
#define sheepGainFromFood 4
#define sheepReproduce 4 //%
// #define initWolveNum 200
#define wolveGainFromFood 20
#define wolveReproduce 5 //%
#define Grass 1
// #define initGrass 5000
#define grassRegrowth 30 //time
#include <algorithm>
#define merror 0.00001
#define reproduceThreshold 2.0
#define debug 0
// #define visualization 1
// #define outputPopulation 1
#include <vector>
#include <unordered_set>
#include <omp.h>
#include <fstream>
#include <iterator>
#include <chrono>
#include <random>

#ifdef visualization
#include "hdf5.h"
#endif
static const char filename[] = "animal.h5";

using namespace std;
int half;
// int half = sqrt(N);
int tot_sheep=0;
int tot_wolve = 0;
int tot_grass = 0;
int N,T,initSheepNum,initWolveNum,initGrass;

std::chrono::microseconds TimeTotal;
std::chrono::microseconds TimeParallel;
std::chrono::microseconds TimeEat;
std::chrono::microseconds TimeInit;
std::chrono::microseconds TimeRenewStat;

float randFloat(float low, float high)
{
    thread_local static std::random_device rd;
    thread_local static std::mt19937 rng(rd());
    thread_local std::uniform_real_distribution<float> urd;
    return urd(rng, decltype(urd)::param_type{low, high});
}

#ifdef visualization
int mat2hdf5(double *sdata, double *wdata, int *gdata, int *count, int *setting, const char *filename)
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
#endif

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
        Grassclass(){sx = 0; sy = 0; number=0;}
        Grassclass(int x, int y, int count){sx = x; sy = y; number=count;}
        int countdown();
        int x() { return sx; };
        int y() { return sy; };
        int gNum() { return number; };
        void sNum(int count) { number = count;};
        void reset(){number = (-1) * grassRegrowth;};
        void assign(int x, int y, int count){sx = x; sy = y; number=count;}
};
int Grassclass::countdown(){
    bool out = 0;
    if (number==0){
        // tot_grass++;
        out = 1;
    }
    if (number < 1)
        number++;
    return out;
}

class RandomWalk{
    int sx,sy,sd;
    public:
        RandomWalk(int x, int y, int d);
        void move();
        int x(){return sx;};
        int y(){return sy;};
        int d(){return sd;};
        void assign(int x, int y, int d){sx = x; sy=y; sd= d;};
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
    // int randdir = (float)rand()  *2.99; //Brownian movement
    int randdir = (float)randFloat(0., 1.) *2.99; //Brownian movement

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
     	double energy_rand = std::max(1.0, (double)2 * gainFood * (double)randFloat(0., 1.) );
     	energy = std::max(1.0, energy_rand);
		};
        Animal(int x, int y, int d, int flag, float inenergy): RandomWalk(x,y,d){
		sflag = flag;
     	gainFood = flag == 0 ? sheepGainFromFood : wolveGainFromFood;
     	energy = inenergy;
		};
		void reduceEnergy(){energy -= 1;};
        void addEnergy() { energy += gainFood; }
        void divideEnergy() { energy /= 2;}
        float gEnergy() { return energy; };
        void sEnergy(float en) {energy = en;}
        int gFlag() {return sflag;};
};
void init_sheep_wolve(vector<Animal> &sheeps, int flag)
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
        int randint = (float)randFloat(0., 1.) * (N - 1);
        direction = 7.99 * (float)randFloat(0., 1.) ;
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
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < N; i++){
            y = i / half;
            x = i - y * half;
            count = (-1) * std::max(1.0, (double)grassRegrowth * (float)randFloat(0., 1.)  );
            grasses[i].assign(x,y,count);
            // Grassclass grass(x, y, count);
            // grasses.push_back(grass);
        }
    }
    std::unordered_set<int> grassSet;
    while (grassNum < initGrass)
    {   
        int randint = (float)randFloat(0., 1.) * (N - 1);
        if (grassSet.find(randint) != grassSet.end())
            continue;
        grassSet.insert(randint);
        grasses[randint].sNum(1);
        grassNum++;
        tot_grass++;
    }

}
int eatGrass(Animal &sheep, vector<Grassclass> &grasses)
{
    int out= 0;
    // omp_lock_t writelock;
    // omp_init_lock(&writelock);
    // omp_set_lock(&writelock);
    if (grasses[sheep.x() + half * sheep.y()].gNum() == 1)
    {
        sheep.addEnergy();
        grasses[sheep.x() + half * sheep.y()].reset();
        out = 1;
    }
    // omp_unset_lock(&writelock);
    return out;
}
int death(vector<Animal> &animals, Animal &animal)
{
    int out = 0;
    if (animal.gEnergy() <= merror)
    {
        if (animal.gFlag() == 0)
        {
            // tot_sheep -= 1;
            if (debug == 1)
                printf("dead one sheep\n");
        }
        else
        {
            // tot_wolve -= 1;
            if (debug == 1)
                printf("dead one wolve\n");
        }
        out = 1;
        animal.sEnergy(0);
    }
    return out;
}
int eatSheep(Animal &wolf, vector<Animal> &sheeplist)
{
    int out = 0;
    for (std::vector<Animal>::iterator it = sheeplist.begin(); it != sheeplist.end(); ++it)
    {
        if ((wolf.x() == (*it).x()) && (wolf.y() == (*it).y()) && ((*it).gEnergy()> merror))
        {
            wolf.addEnergy();
            (*it).sEnergy(-1);
            out = death(sheeplist, (*it));
            // it = sheeplist.erase(it);
            // tot_sheep--;
            break;
        }
    }
    return out;
}
int eatSheep_parallel(Animal &wolf, Animal &sheep, vector<Animal> &sheeplist){
    int out = 0;
    if ((wolf.x() == sheep.x()) && (wolf.y() == sheep.y()) && (sheep.gEnergy()> merror))
    {
        wolf.addEnergy();
        sheep.sEnergy(-1);
        out = death(sheeplist, sheep);
    }
    return out;

}
int eatSheep_overlap(Animal &wolf, vector<Animal> &sheeplist, vector<int>& sheepPos ){
    int vec_size = sheepPos.size();
    int out = 0;
    for (int i = 1; i < vec_size; i++){
        if (sheeplist[i].gEnergy()> merror){
            wolf.addEnergy();
            sheeplist[i].sEnergy(-1);
            out = death(sheeplist, sheeplist[i]);
            break;
        }
    }
    return out;
}

int create_animal(vector<Animal> &animals, Animal &animal)
{
    int xlist[8] = {0}, ylist[8] = {0};
    int out = 0;
    // assume they can move randomly in 8 direction
    gen_surroundxy(animal.x(), animal.y(), xlist, ylist);
    int i = (float)randFloat(0., 1.)  *7.99; //choose one direction to move
    valify_move(&xlist[i], &ylist[i]);
    int d = (float) randFloat(0., 1.)  *7.99; // initialize direction
    Animal newAnimal(xlist[i], ylist[i], d, animal.gFlag());
    animals.push_back(newAnimal);
    if (animal.gFlag()==0){
        // tot_sheep++;
        out = 1;
    }else{
        // tot_wolve++;
        out = 1;
    }
    return out;
}
int reproduce(vector<Animal> &animals, Animal &animal)
{
    float randint = (float)randFloat(0., 1.)  *100;
    int cond_reproduce_rate = (animal.gFlag() == 0) ? sheepReproduce : wolveReproduce;

    int out = 0;
    if ((randint < (float)cond_reproduce_rate) && ((int)animals.size() < N))
    {
        animal.divideEnergy();
        out = create_animal(animals, animal);
    }
    return out;
}

void get_state(vector<Animal> &sheeps, vector<Animal> &wolves, vector<Grassclass> &grasses)
{
    int tot_sheep_acc = sheeps.size();
    int tot_wolf_acc = wolves.size();
    int tot_grass_acc = 0;
    int i, vec_size = grasses.size();
#pragma omp parallel
    {
        int local_tot = 0;
        #pragma omp for
        // for (Grassclass a : grasses)
        for (i = 0; i < vec_size; i++)
        {
            if (grasses[i].gNum() > 0)
            {
                local_tot+=1;
            }
        }
        #pragma omp critical
        tot_grass_acc += local_tot;
    }
    tot_grass = tot_grass_acc;
    tot_sheep = tot_sheep_acc;
    tot_wolve = tot_wolf_acc;
    // printf("s %d, w %d, g %d\n", tot_sheep, tot_wolve, tot_grass);
    // printf("s (%d, %d), w (%d, %d), g (%d, %d)\n", tot_sheep, tot_sheep_acc, tot_wolve, tot_wolf_acc, tot_grass, tot_grass_acc);
}
void ask_sheep(vector<Animal>& sheeplist, vector<Grassclass>& grasslist)
{
    int vec_size = sheeplist.size();
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < vec_size; i++)
        {
            sheeplist[i].move();
            if (Grass == 1)
            {
                sheeplist[i].reduceEnergy();
            }
        }
    }
    
    int local_tot = 0;
    for (int i = 0; i < vec_size; i++)
    {
        if (Grass == 1)
        {
            
            int out = eatGrass(sheeplist[i], grasslist);
            local_tot -= out;
        }
    }
    tot_grass += local_tot;

    #pragma omp parallel
    {
        std::vector<Animal> sheeplist_tmp;
        int local_tot = 0;
        #pragma omp for
        for (int i = 0; i < vec_size; i++)
        {
            if (sheeplist[i].gEnergy() > reproduceThreshold)
            {
                int out = reproduce(sheeplist_tmp, sheeplist[i]);
                local_tot += out;
            }
        }
        #pragma omp critical
        {
            sheeplist.insert(sheeplist.end(), sheeplist_tmp.begin(), sheeplist_tmp.end());
            tot_sheep += local_tot;
        }
    }
    #pragma omp parallel
    {
        int local_tot = 0;
        #pragma omp for
        for (int i = 0; i < vec_size; i++)
        {
            int out = death(sheeplist, sheeplist[i]);
            local_tot -= out;
        }
        #pragma omp critical
        tot_sheep += local_tot;
    }
}

void ask_wolf(vector<Animal> &wolflist, vector<Animal> &sheeplist)
{
    int vec_size = wolflist.size();
    #pragma omp parallel
    {
    #pragma omp for
        for (int i = 0; i < vec_size; i++)
        {
            wolflist[i].move();
            wolflist[i].reduceEnergy();
        }
    }
    for (int i = 0; i < vec_size; i++)
    {
        int out = eatSheep(wolflist[i], sheeplist);
        tot_sheep -= out;
    }
    #pragma omp parallel
    {
        std::vector<Animal> wolflist_tmp;
        int local_tot = 0;
        #pragma omp for
        for (int i = 0; i < vec_size; i++)
        {
            if (wolflist[i].gEnergy() > reproduceThreshold)
            {
                int out = reproduce(wolflist_tmp, wolflist[i]);
                local_tot += out;
            }   
        }
        #pragma omp critical
        {
            wolflist.insert(wolflist.end(), wolflist_tmp.begin(), wolflist_tmp.end());
            tot_wolve += local_tot;
        }
    }
    #pragma omp parallel
    {
        int local_tot = 0;
        #pragma omp for
        for (int i = 0; i < vec_size; i++)
        {
            int out = death(wolflist, wolflist[i]);
            local_tot -= out;
        }
        #pragma omp critical
        tot_wolve += local_tot;
    }
}
// move, reduceEnergy, death, reproduce
void ask_animal(int agent, vector<Animal>& animallist, vector<Animal>& newanimallist, int start, int end)
{
    int vec_size = end;
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = start; i < vec_size; i++)
        {
            animallist[i].move();
            if ((agent!=1) || ( (agent==1) && (Grass == 1)))
            {
                animallist[i].reduceEnergy();
            }
        }
    }
    #pragma omp parallel
    {
        int local_tot = 0;
        #pragma omp for
        for (int i = start; i < vec_size; i++)
        {
            int out = death(animallist, animallist[i]);
            local_tot -= out;
        }
    }

    #pragma omp parallel
    {
        std::vector<Animal> animallist_tmp;
        int local_tot = 0;
        #pragma omp for
        for (int i = start; i < vec_size; i++)
        {
            if (animallist[i].gEnergy() > reproduceThreshold) // merror) //
            {
                int out = reproduce(animallist_tmp, animallist[i]);
                local_tot += out;
            }
        }
        #pragma omp critical
        {
            newanimallist.insert(newanimallist.end(), animallist_tmp.begin(), animallist_tmp.end());
            // tot_sheep += local_tot;
        }
    }
    
}
void ask_patch(vector<Grassclass> &grasslist, int start, int end)
{
     int vec_size = end;
    #pragma omp parallel
    {
        int local_tot = 0;
        #pragma omp for
         for (int i = start; i < vec_size; i++)
         {
             int out = grasslist[i].countdown();
             local_tot += out;
         }
         #pragma omp critical
         tot_grass += local_tot;
     }
}
void save2mat(double *matTime, vector<Animal> &mat, int t)
{
    // int ID = omp_get_thread_num();
    // printf("mat %d\n", ID);

    int base = N * t;
    int vec_size = mat.size();
    #pragma omp parallel
    {
        #pragma omp for
        // for (std::vector<Animal>::iterator it = mat.begin(); it != mat.end(); ++it)
        for (int i = 0; i < vec_size; i++)
        {
            matTime[mat[i].x() + half * mat[i].y() + base] = mat[i].gEnergy();
            // matTime[(*it).x() + half * (*it).y() + base] = (*it).gEnergy();
        }
    }
}
void save2matInt(int *matTime, vector<Grassclass> &grasslist, int t)
{
    int base = N * t;
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = 0; i < N; i++)
        {
            matTime[i + base] = grasslist[i].gNum();
        }
    }
}
void renew_vector(vector<Animal> &mat,vector<Animal>& newMat, int start, int end){
    int vec_size = end;
    #pragma omp parallel
    {
        vector<Animal> local;
        #pragma omp for
        for (int i = start; i < vec_size; i++){
            if (mat[i].gEnergy() > merror){
                local.push_back(mat[i]);
            }
        }
        #pragma omp critical
        newMat.insert(newMat.end(), local.begin(), local.end());
    }
}

void animalvec2mat(vector<Animal> &mat, float *energy, int *xydf, int start, int end){
    int vec_size = end;
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = start; i < vec_size; i++)
        {
            energy[i] = mat[i].gEnergy();
            xydf[i*4]= mat[i].x();
            xydf[i*4+1]= mat[i].y();
            xydf[i*4+2]= mat[i].d();
            xydf[i*4+3] = mat[i].gFlag();
        }
    }
}
void mat2animalvec(vector<Animal> &mat, float* energy, int* xydf, int energySize){
    int vec_size = energySize;
    #pragma omp parallel
    {
        vector<Animal> local;
        #pragma omp for
        for (int i = 0; i < vec_size; i++)
        {
            Animal animal(xydf[i*4],xydf[i*4+1],xydf[i*4+2],xydf[i*4+3], energy[i]);
            local.push_back(animal);
        }
        #pragma omp critical
        mat.insert(mat.end(), local.begin(), local.end());
    }
}

void grassvec2mat(vector<Grassclass> &mat, int* xyn, int start, int end){
    int vec_size = end;
    #pragma omp parallel
    {
        #pragma omp for
        for (int i = start; i < vec_size; i++)
        {
            xyn[i*3]= mat[i].x();
            xyn[i*3+1]= mat[i].y();
            xyn[i*3+2]= mat[i].gNum();
        }
    }
}
void mat2grassvec(vector<Grassclass> &mat, int* xyn, int start, int end){
    int vec_size = end;
    if (mat.size()==0){
        #pragma omp parallel
        {
            vector<Grassclass> local;
            #pragma omp for
            for (int i = start; i < vec_size; i++)
            {
                Grassclass grass(xyn[i*3],xyn[i*3+1],xyn[i*3+2]);
                local.push_back(grass);
            }
            #pragma omp critical
            mat.insert(mat.end(), local.begin(), local.end());
        }
    }else{
        #pragma omp parallel
        {
            #pragma omp for
            for (int i = start; i < vec_size; i++)
            {
                mat[i].assign(xyn[i*3],xyn[i*3+1],xyn[i*3+2]);
            }
        }
    }
    
}

void animal_eat(vector<Animal> &sheeplist, vector<Animal> &wolflist, vector<Grassclass> &grasslist){
    int vec_size;
    vec_size = sheeplist.size();
    for (int i = 0; i < vec_size; i++)
    {
        if (Grass == 1){
            if (sheeplist[i].gEnergy()>merror)
                eatGrass(sheeplist[i], grasslist);
        }
        
    }
    vec_size = wolflist.size();
    for (int i = 0; i < vec_size; i++)
    {
        if (wolflist[i].gEnergy()>merror)
            eatSheep(wolflist[i], sheeplist);
    }
    // vec_size = sheeplist.size();
    // vector<vector<int>> sheep2d(N, vector<int> (1, 0));
    // vector<vector<int>> wolf2d(N, vector<int> (1, 0));
    // #pragma omp parallel
    // {
    //     #pragma omp for
    //     for (int i = 0; i < vec_size; i++)
    //     {
    //         if ((Grass == 1) && (sheeplist[i].gEnergy()>merror))
    //         {
    //             eatGrass(sheeplist[i], grasslist);
    //         }
    //     }
    // }
    // // speed bottle neck
    // // store the index of sheep in the vector for look up while wolf eat sheep
    // for (int i = 0; i < vec_size; i++){
    //     if (sheeplist[i].gEnergy()>merror){
    //         sheep2d[sheeplist[i].x() + half * sheeplist[i].y()].push_back(i);
    //     }
    // }
    // vec_size = wolflist.size();
    // for (int i = 0; i < vec_size; i++)
    //     wolf2d[wolflist[i].x() + half * wolflist[i].y()].push_back(i);

    // // parallel for non overlap agent
    // #pragma omp parallel
    // {
    //     #pragma omp for
    //     for (int i = 0; i < vec_size; i++)
    //     {
    //         int pos = wolflist[i].x() + wolflist[i].y() * half;
    //         if ((wolflist[i].gEnergy() > merror) 
    //         && (wolf2d[pos].size()==2) 
    //         && (sheep2d[pos].size()==2)){
    //             eatSheep_parallel( wolflist[i], sheeplist[sheep2d[pos][1]], sheeplist);
    //         }
    //     }
    // }
    // // serial for overlapping agents
    // for (int i = 0; i < vec_size; i++)
    // {
    //     int pos = wolflist[i].x() + wolflist[i].y() * half;
    //     if ((wolflist[i].gEnergy() > merror)
    //     && (sheep2d[pos].size() >= 2) 
    //     && ((wolf2d[pos].size() != 2) ||
    //         (sheep2d[pos].size() != 2) ) ) {
    //         eatSheep_overlap( wolflist[i], sheeplist, sheep2d[pos]);
    //     }
    // }


}


void act_master(vector<Animal> &sheeplist, vector<Animal> &wolflist, vector<Grassclass> &grasslist, int time){


    auto StartParallel = std::chrono::high_resolution_clock::now();
    // serial
    // //ask grass
    if (Grass != 0)
        ask_patch(grasslist, 0, N);

    // ask animal for move, reduce energy, death, reproduce
    vector<Animal> newsheeplist;
    vector<Animal> newwolflist;
    int size_s =  sheeplist.size();
    int size_w =  wolflist.size();
    ask_animal(1, sheeplist, sheeplist,0,size_s);
    ask_animal(2, wolflist, wolflist,0,size_w);

    auto EndParallel = std::chrono::high_resolution_clock::now();
    
    TimeParallel += std::chrono::duration_cast<std::chrono::microseconds>(EndParallel-StartParallel);

    auto StartEat = std::chrono::high_resolution_clock::now();
    // animal eat grass or eat sheep
    animal_eat(sheeplist, wolflist, grasslist);
    auto EndEat = std::chrono::high_resolution_clock::now();
    TimeEat += std::chrono::duration_cast<std::chrono::microseconds>(EndEat-StartEat);

    StartParallel = std::chrono::high_resolution_clock::now();
    // pick alive animal in the vector
    vector<Animal> new_sheeplist, new_wolflist;
    renew_vector(sheeplist, new_sheeplist, 0, sheeplist.size());
    renew_vector(wolflist, new_wolflist, 0, wolflist.size());
    sheeplist = new_sheeplist;
    wolflist = new_wolflist;
    // count green grass
    get_state(sheeplist, wolflist, grasslist);

    EndParallel = std::chrono::high_resolution_clock::now();
    TimeRenewStat += std::chrono::duration_cast<std::chrono::microseconds>(EndParallel-StartParallel);
    TimeParallel += std::chrono::duration_cast<std::chrono::microseconds>(EndParallel-StartParallel);


}
int main(int argc, char** argv)
{
    srand(time(NULL));
    if (argc < 6)
    {
        std::cerr << "expected: N, T, initSheepNum, initWolveNum, initGrass \n";
        return 1;
    }
    N = atoi(argv[1]);
    T = atoi(argv[2]);
    initSheepNum = atoi(argv[3]);
    initWolveNum = atoi(argv[4]);
    initGrass = atoi(argv[5]);
    half = sqrt(N);
    printf("N %d, T %d, iS %d, iW %d, iG %d\n", N,T,initSheepNum,initWolveNum,initGrass);

    assert(N>=0 && half * half == N && "N should be able to square root and greater equals to zero");
    assert(initSheepNum>=0 && initSheepNum <= N && "sheep init number should be smaller than N and greater equals to zero");
    assert(initWolveNum>=0 && initWolveNum <= N && "wolve init number should be smaller than N and greater equals to zero");
    assert(initGrass>=0 && initGrass <= N && "grass init number should be smaller than N and greater equals to zero");


    TimeTotal = std::chrono::microseconds::zero();
    TimeParallel = std::chrono::microseconds::zero();
    TimeEat = std::chrono::microseconds::zero();
    TimeRenewStat = std::chrono::microseconds::zero();

	std::vector<Animal> sheeplist;
    std::vector<Animal> wolflist;
    std::vector<Grassclass> grasslist(N);

    auto Startall = std::chrono::high_resolution_clock::now();


#ifdef visualization
    double *sheep = new double [N * T * sizeof(double)];
    double *wolve = new double [N * T * sizeof(double)];
    int *grass = new int [N * T * sizeof(int)];
#endif
    int *animalNum = new int [3 * T * sizeof(int)];
    for (int i = 0; i < 3 * T; i++)
        animalNum[i] = 0;
    vector<int> animalNumVec(3 * T,0);
    int *setting = new int [2 * sizeof(int)];
    setting[0] = N;
    setting[1] = T;

    auto StartMultiply = std::chrono::high_resolution_clock::now();
    //init
	init_sheep_wolve(sheeplist, 0);
    init_sheep_wolve(wolflist, 1);
    if (Grass != 0)
    {
        init_grass(grasslist);
    }
    else{
        tot_grass = N;
    }
    auto EndMultiply = std::chrono::high_resolution_clock::now();
    TimeInit = std::chrono::duration_cast<std::chrono::microseconds>(EndMultiply-StartMultiply);
    TimeParallel += std::chrono::duration_cast<std::chrono::microseconds>(EndMultiply-StartMultiply);

    for (int t = 0; t < T; t++)
    {
    #ifdef visualization
            save2mat(sheep, sheeplist, t);
            save2mat(wolve, wolflist, t);
            save2matInt(grass, grasslist, t);
    #endif
            animalNum[0 + t * 3] = tot_sheep;
            animalNum[1 + t * 3] = tot_wolve;
            animalNum[2 + t * 3] = tot_grass;
            animalNumVec[0 + t * 3] = tot_sheep;
            animalNumVec[1 + t * 3] = tot_wolve;
            animalNumVec[2 + t * 3] = tot_grass;
            act_master(sheeplist, wolflist, grasslist, t);
            // printf("t %d, s %d, w %d, g %d\n",t, tot_sheep, tot_wolve, tot_grass);
        #ifdef outputPopulation
        std::ofstream f("population_dynamic.txt");
        int count = 0;
        f << T << ", " << N << ", " << 0 << "\n";
        for(vector<int>::const_iterator i = animalNumVec.begin(); i != animalNumVec.end(); ++i) {
            if (count % 3 == 2){
                f << *i <<  '\n';
            }else{
                f << *i <<  ", ";
            }
            count++;
        }
        #endif
    }
#ifdef visualization
	mat2hdf5(sheep, wolve, grass, animalNum, setting, filename);
#endif
    
    auto Endall = std::chrono::high_resolution_clock::now();
    TimeTotal = std::chrono::duration_cast<std::chrono::microseconds>(Endall-Startall);

    std::cout << "Total time:           " << std::setw(12) << TimeTotal.count() << " us\n";
    std::cout << "Total init:           " << std::setw(12) << TimeInit.count() << " us\n";
    std::cout << "Act time:             " << std::setw(12) << TimeParallel.count() << " us\n";
    std::cout << "Eat time:             " << std::setw(12) << TimeEat.count() << " us\n";
    std::cout << "Renew & State time:   " << std::setw(12) << TimeRenewStat.count() << " us\n";

#ifdef visualization
    delete [] sheep;
    delete [] grass;
    delete [] wolve;
#endif
    delete [] animalNum;
    delete [] setting;

    return 0;
}

