#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#define N 10000
#define T 100
#define initSheepNum 400
#define sheepGainFromFood 4
#define sheepReproduce 4 //%
#define initWolveNum 200
#define wolveGainFromFood 20
#define wolveReproduce 5 //%
#define Grass 1
#define initGrass 5000
#define grassRegrowth 30 //time
#include <algorithm>
#define merror 0.00001
#define reproduceThreshold 2.0
#define debug 0
// #define visualization 1
#define outputPopulation 1
#include <vector>
#include <unordered_set>
#include <mpi.h>
#include <omp.h>
#include <fstream>
#include <iterator>

#ifdef visualization
#include "hdf5.h"
#endif
static const char filename[] = "animal.h5";

using namespace std;
int half = sqrt(N);
int tot_sheep=0;
int tot_wolve = 0;
int tot_grass = 0;

#include <random>
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
        // tot_grass--;
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
        // animal.sEnergy(animals[animals.size() - 1].gEnergy());
        // animal.assign(animals[animals.size() - 1].x(), animals[animals.size() - 1].y(), animals[animals.size() - 1].d());
        // animals.pop_back();
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
        int idx = sheepPos[i];
        if (sheeplist[idx].gEnergy()> merror){
            wolf.addEnergy();
            sheeplist[idx].sEnergy(-1);
            out = death(sheeplist, sheeplist[idx]);
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
    if ((randint < (float)cond_reproduce_rate) && (animals.size() < N))
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
void gen_agent_rankid(int agent, int dim, int* idlist){
    for (int j= 0; j < dim; j++){
            idlist[j] = agent + j*3;
    }
}
void gen_agent_processor(int* dim){
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size>2){
        for (int i =0;i<world_size; i++){
            dim[i%3] += 1;
        }
    }
}
void send_animal(vector<Animal>& animallist){
    int size = animallist.size();
    float* energy = new float[size * sizeof(float)];
    int* xydf = new int[size * 4 * sizeof(int)];
    animalvec2mat(animallist, energy, xydf, 0, size);
    MPI_Send(energy, size, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(xydf, size*4, MPI_INT, 0, 0, MPI_COMM_WORLD);
    delete [] energy;
    delete [] xydf;

}
void receive_animal(vector<Animal>& animallist, int agent){
    int iniNum = agent == 1? initSheepNum:initWolveNum;
    float* energy = new float [iniNum * sizeof(float)];
    int* xydf = new int [iniNum * 4 * sizeof(int)];
    MPI_Recv(energy, iniNum, MPI_FLOAT, agent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(xydf, iniNum * 4, MPI_INT, agent, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    mat2animalvec(animallist, energy, xydf, iniNum);
    delete [] energy;
    delete [] xydf;
}
void initialize_parallel(vector<Grassclass> &grasslist, vector<Animal> &sheeplist, vector<Animal> &wolflist){
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    if (world_size < 3){
        if(world_rank == 0){
            init_sheep_wolve(sheeplist, 0);
            init_sheep_wolve(wolflist, 1);
            if (Grass != 0)
            {
                init_grass(grasslist);
            }
            else{
                tot_grass = N;
            }
        }
    }else{
        if (world_rank == 0){
            if (Grass != 0)
            {
                init_grass(grasslist);
                tot_grass = initGrass;
            }
            receive_animal(sheeplist, 1);
            receive_animal(wolflist, 2);
            tot_sheep = initSheepNum;
            tot_wolve = initWolveNum;

        }
        if (world_rank == 1){
            init_sheep_wolve(sheeplist, 0);
            send_animal(sheeplist);
        }
        if (world_rank == 2){
            init_sheep_wolve(wolflist, 1);        
            send_animal(wolflist);    
        }
    }
}
template <class S>
void my_scatter_send(S* data, int* count, MPI_Datatype datatype, int root, int* idlist, int numSub, MPI_Comm communicator){
    int world_rank;
    MPI_Comm_rank(communicator, &world_rank);
    int world_size;
    MPI_Comm_size(communicator, &world_size);

    if (world_rank == root) {
        // If we are the root process, scatter our data to everyone
        for (int i = 0; i < numSub; i++) {
            if (idlist[i] != world_rank) {
                MPI_Send(data+i*(*count), (*count), datatype, idlist[i], 0, communicator);
            }
        }
    } 
}
template <class S>
void my_scatter_rec(S** data, int* count, MPI_Datatype datatype, int root, int* idlist, int numSub, MPI_Comm communicator){
    // If we are a receiver process, receive the data from the root
    int number_amount;
    MPI_Status status;
    MPI_Probe(root, 0, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, datatype, &number_amount);
    *data = new S [sizeof(S) * number_amount];
    *count = number_amount;
    MPI_Recv(*data, number_amount, datatype, root, 0, communicator, MPI_STATUS_IGNORE);
}

template <class S>
void my_gather(S* data, int* count, MPI_Datatype datatype, int root, int* idlist, int numSub, MPI_Comm communicator){
    int world_rank;
    MPI_Comm_rank(communicator, &world_rank);
    int world_size;
    MPI_Comm_size(communicator, &world_size);
    MPI_Status status;
    int number_amount;
    int acc_count = (*count);

    if (world_rank == root) {
        // If we are the root process, receive our data from everyone
        int i;
        for (i = 0; i < numSub; i++) {
            // root always either 0, 1 or 2
            if (idlist[i] != world_rank) {
                MPI_Probe(idlist[i], 0, MPI_COMM_WORLD, &status);
                MPI_Get_count(&status, datatype, &number_amount);
                MPI_Recv(data + acc_count, number_amount, datatype, idlist[i], 0, communicator, MPI_STATUS_IGNORE);
                acc_count += number_amount;
            }
        }
        (*count) = acc_count;
    } else {
        // If we are a receiver process, send the data to the root
        MPI_Send(data, *count, datatype, root, 0, communicator);
    }

}

void act_grass_master(vector<Grassclass>& grasslist, int* dim, int* idlist, MPI_Comm communicator){
    int size = grasslist.size();
    int num_elements_per_proc = size/dim[0];
    int num_elements_per_proc3 = num_elements_per_proc * 3;
    int* xyn = new int[N * 3 * sizeof(int)];
    
    // send data from proc 0 to rest proc included root proc
    grassvec2mat(grasslist, xyn, 0, size);
    my_scatter_send(xyn, &num_elements_per_proc3, MPI_INT, 0, idlist, dim[0], communicator);

    // run grass patch for root's part
    ask_patch(grasslist, 0, num_elements_per_proc);
    // corner case taken cared by root proc for the agent in the far back
    if (num_elements_per_proc*dim[0] != size){
        ask_patch(grasslist, num_elements_per_proc*dim[0], size);
    }
    // receive data from slaves
    my_gather(xyn, &num_elements_per_proc3, MPI_INT, 0, idlist, dim[0], MPI_COMM_WORLD);
    mat2grassvec(grasslist, xyn, num_elements_per_proc, num_elements_per_proc*dim[0]);
    delete [] xyn;
}
void act_grass_slaves(int* dim, int* idlist){
    // receive from root
    int* sub_xyn;
    int num_elements_per_proc;

    // receive data from master 0 
    my_scatter_rec(&sub_xyn, &num_elements_per_proc, MPI_INT, 0, idlist, dim[0], MPI_COMM_WORLD);
    // create sub grasslist in slave processor
    vector<Grassclass> sub_grasslist;
    num_elements_per_proc /= 3;
    mat2grassvec(sub_grasslist, sub_xyn, 0, num_elements_per_proc);
    // run ask patch
    ask_patch(sub_grasslist, 0, num_elements_per_proc);
    // slave proc (sub_xyn) to vector and send back to root
    grassvec2mat(sub_grasslist, sub_xyn, 0, num_elements_per_proc);
    num_elements_per_proc *= 3;
    my_gather(sub_xyn, &num_elements_per_proc, MPI_INT, 0, idlist, dim[0], MPI_COMM_WORLD);
    delete [] sub_xyn;
}
void act_animal_slaves(int agent, int* dim, int* idlist){
    // receive from root
    float* sub_energy;
    int* sub_xydf;
    int num_elements_per_proc;
    int num_elements_per_proc4;
    
    //receive from master 0
    my_scatter_rec(&sub_energy, &num_elements_per_proc, MPI_FLOAT, 0, idlist, dim[agent], MPI_COMM_WORLD);
    my_scatter_rec(&sub_xydf, &num_elements_per_proc4, MPI_INT, 0, idlist, dim[agent], MPI_COMM_WORLD);
    // create sub grasslist in slave processor
    vector<Animal> sub_animal;
    vector<Animal> sub_newanimal;
    mat2animalvec(sub_animal, sub_energy, sub_xydf, num_elements_per_proc);
    // run animal move, reduction, death, reproduce for slaves
    ask_animal(agent, sub_animal, sub_newanimal, 0, num_elements_per_proc);
    // clean death animal
    vector<Animal> sub_animal_clean;
    // renew_vector(sub_animal, sub_animal_clean, 0, num_elements_per_proc);
    // // append new born animal to original animal
    // sub_animal_clean.insert(sub_animal_clean.end(), sub_newanimal.begin(), sub_newanimal.end());

    // no clean death animal
    sub_animal.insert(sub_animal.end(), sub_newanimal.begin(), sub_newanimal.end());
    sub_animal_clean = sub_animal;

    // slave proc (sub_xydf, sub_energy) to vector and send back to master 0
    num_elements_per_proc = sub_animal_clean.size(); //if clean
    num_elements_per_proc4 = num_elements_per_proc * 4;
    animalvec2mat(sub_animal_clean, sub_energy, sub_xydf, 0, num_elements_per_proc);

    // send size first
    int send_num = 1;
    int* send_size = new int [sizeof(int)*send_num];
    send_size[0] = num_elements_per_proc;
    my_gather(send_size, &send_num, MPI_INT, 0, idlist, dim[agent], MPI_COMM_WORLD);
    // send data later
    my_gather(sub_energy, &num_elements_per_proc, MPI_FLOAT, 0, idlist, dim[agent], MPI_COMM_WORLD);
    my_gather(sub_xydf, &num_elements_per_proc4, MPI_INT, 0, idlist, dim[agent], MPI_COMM_WORLD);
    delete [] sub_energy;
    delete [] sub_xydf;

}
void act_animal_corner(vector<Animal>& animal, vector<Animal>& animal_corner, int agent, int* dim){
    int size = animal.size();
    int num_elements_per_proc = size/dim[agent];
    vector<Animal> newanimallist;
        
    // corner case taken cared by root proc for the agent in the far back
    if (num_elements_per_proc * dim[agent] != size){
        ask_animal(agent, animal, newanimallist, num_elements_per_proc*dim[agent], size);
    }
    // clean death sheep
    // renew_vector(animal, animal_corner, num_elements_per_proc*dim[agent], size);

    // append new born animal to original animal
    animal_corner.insert(animal_corner.end(), newanimallist.begin(), newanimallist.end());

}
void animal_scatter(vector<Animal>& animal,int agent, int* dim, int* idlist, MPI_Comm communicator){
    int size = animal.size();
    int num_elements_per_proc = size/dim[agent];
    int num_elements_per_proc4 = num_elements_per_proc * 4;
    float* energy = new float[size * sizeof(float)];
    int* xydf = new int[size * 4 * sizeof(int)];
    // printf("dim %d, idlist %d,%d\n",dim[agent], idlist[0],idlist[1]);
    
    // send data from master 0 to rest proc, sheep (1,4,7..)/ wolf(2,5,8,...)
    animalvec2mat(animal, energy, xydf, 0, size);
    my_scatter_send(energy, &num_elements_per_proc, MPI_FLOAT, 0, idlist, dim[agent], communicator);
    my_scatter_send(xydf, &num_elements_per_proc4, MPI_INT, 0, idlist, dim[agent], communicator);

    delete [] energy;
    delete [] xydf;
}

void animal_gather(vector<Animal>& animal,int agent, int* dim, int* idlist, MPI_Comm communicator){
    // new born
    // the maximum new born animals two times more than original numbers.
    int newBornNumber = 0;
    int* new_size = new int [dim[agent] * sizeof(int)];
    int new_dim = 0;
    int new_num_elements_per_proc = 0;
    int new_num_elements_per_proc4 = 0;
    vector<Animal> outAnimal;

    // receive size of data from each processor
    my_gather(new_size, &new_dim, MPI_INT, 0, idlist, dim[agent], MPI_COMM_WORLD);
    
    for (int i = 0; i< new_dim; i++){
        newBornNumber += new_size[i];
    }
    // printf("newBornNumber %d\n",newBornNumber);
    float* new_energy = new float[newBornNumber * sizeof(float)];
    int* new_xydf = new int[newBornNumber * 4 * sizeof(int)];
    
    // receive data from each processor with dynamic size
    my_gather(new_energy, &new_num_elements_per_proc, MPI_FLOAT, 0, idlist, dim[agent], MPI_COMM_WORLD);
    my_gather(new_xydf, &new_num_elements_per_proc4, MPI_INT, 0, idlist, dim[agent], MPI_COMM_WORLD);

    // printf("total animal %d\n", new_num_elements_per_proc);
    mat2animalvec(outAnimal, new_energy, new_xydf, new_num_elements_per_proc);
    animal = outAnimal;
    delete [] new_energy;   
    delete [] new_xydf;
}
void animal_eat(vector<Animal> &sheeplist, vector<Animal> &wolflist, vector<Grassclass> &grasslist){
    int vec_size;
    vec_size = sheeplist.size();
    // vector<vector<int>> sheep2d(N, vector<int> (1, 0));
    // vector<vector<int>> wolf2d(N, vector<int> (1, 0));
    
    // parallel eatgrass while one sheep at one position
    // #pragma omp parallel
    // {
    //     #pragma omp for
    //     for (int i = 0; i < vec_size; i++)
    //     {
    //         //int pos = sheeplist[i].x() + sheeplist[i].y() * half;
    //         if ((Grass == 1) )// && (sheeplist[i].gEnergy()>merror)) // && (sheep2d[pos].size()==2))
    //         {
    //             eatGrass(sheeplist[i], grasslist);
    //         }
    //     }
    // }
    // // serial eatgrass while multiple sheeps at one position
    // // for (int i = 0; i < vec_size; i++)
    // // {
    // //     int pos = sheeplist[i].x() + sheeplist[i].y() * half;
    // //     if ((Grass == 1) && (sheeplist[i].gEnergy()>merror) && (sheep2d[pos].size()!=2))
    // //     {
    // //         eatGrass(sheeplist[i], grasslist);
    // //     }
    // // }

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
    //         if ( //(wolflist[i].gEnergy() > merror) 
    //         /*&&*/ (wolf2d[pos].size()==2) 
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

}

void sendTime(int t){
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    for (int j = 1; j < world_size; j++){
        MPI_Send(&t, 1, MPI_INT, j, 44, MPI_COMM_WORLD);
    }
}

void act_master(vector<Animal> &sheeplist, vector<Animal> &wolflist, vector<Grassclass> &grasslist, int time){
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int* dim = new int [3 * sizeof(int)]();
    gen_agent_processor(dim);
    int* idlist = new int [sizeof(int) * dim[0]];
    int* idlist_s = new int [sizeof(int) * dim[1]];
    int* idlist_w = new int [sizeof(int) * dim[2]];

    if (world_size < 3){
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

        // animal eat grass or eat sheep
        animal_eat(sheeplist, wolflist, grasslist);

        // pick alive animal in the vector
        vector<Animal> new_sheeplist, new_wolflist;
        renew_vector(sheeplist, new_sheeplist, 0, sheeplist.size());
        renew_vector(wolflist, new_wolflist, 0, wolflist.size());
        sheeplist = new_sheeplist;
        wolflist = new_wolflist;
        // count green grass
        get_state(sheeplist, wolflist, grasslist);
    }
    else{
        // parallel
        // grass master, scatter grass to slave grass (3,6)
        if (Grass != 0)
        {
            gen_agent_rankid(0, dim[0], idlist);
            act_grass_master(grasslist, dim, idlist, MPI_COMM_WORLD);
        }
        // scatter to slaves sheep (1,4,7)/ wolf (2,5,8)
        gen_agent_rankid(1, dim[1], idlist_s);
        animal_scatter(sheeplist, 1, dim, idlist_s, MPI_COMM_WORLD);
        gen_agent_rankid(2, dim[2], idlist_w);
        animal_scatter(wolflist, 2, dim, idlist_w, MPI_COMM_WORLD);
        // do corner case
        vector<Animal> sheeplist_corner;
        vector<Animal> wolflist_corner;
        // act_animal_corner( sheeplist, sheeplist_corner, 1, dim);
        // act_animal_corner( wolflist, wolflist_corner, 1, dim);

        // gather from sheep, wolf 
        animal_gather(sheeplist, 1, dim, idlist_s, MPI_COMM_WORLD);
        animal_gather(wolflist, 2, dim, idlist_w, MPI_COMM_WORLD);
        // append corner case to main vector
        sheeplist.insert(sheeplist.end(), sheeplist_corner.begin(), sheeplist_corner.end());
        wolflist.insert(wolflist.end(), wolflist_corner.begin(), wolflist_corner.end());
        
        // sheep eat grass, wolf eat sheep
        animal_eat(sheeplist, wolflist, grasslist);

        // pick alive animal in the vector
        vector<Animal> new_sheeplist;
        vector<Animal> new_wolflist;
        renew_vector(sheeplist, new_sheeplist, 0, sheeplist.size());
        renew_vector(wolflist, new_wolflist, 0, wolflist.size());

        sheeplist = new_sheeplist;
        wolflist = new_wolflist;
        // count green grass
        get_state(sheeplist, wolflist, grasslist);
    }

    delete [] dim;
    delete [] idlist;
    delete [] idlist_s;
    delete [] idlist_w;
}
bool main_slave(int world_rank){
    MPI_Status status;
    int* dim = new int [3 * sizeof(int)]();
    int agent = world_rank % 3;
    // break while t is the last frame
    int time;
    MPI_Recv(&time, 1, MPI_INT, 0, 44, MPI_COMM_WORLD, &status);
    
    gen_agent_processor(dim);
    // wolf slaves
    int* idlist = new int [sizeof(int) * dim[agent]];
    gen_agent_rankid(agent, dim[agent], idlist);
    if (agent == 0){
        act_grass_slaves(dim, idlist);
    }
    else{
        act_animal_slaves(agent, dim, idlist);
    }

    if (time == T-1){
        MPI_Finalize();
        return false;
    }
    return true;
}
int main(void)
{
    srand(time(NULL));
    assert(N>=0 && half * half == N && "N should be able to square root and greater equals to zero");
    assert(initSheepNum>=0 && initSheepNum <= N && "sheep init number should be smaller than N and greater equals to zero");
    assert(initWolveNum>=0 && initWolveNum <= N && "wolve init number should be smaller than N and greater equals to zero");
    assert(initGrass>=0 && initGrass <= N && "grass init number should be smaller than N and greater equals to zero");

    MPI_Init(NULL, NULL);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // printf("Hello from process %d out of %d\n", world_rank, world_size);
	std::vector<Animal> sheeplist;
    std::vector<Animal> wolflist;
    std::vector<Grassclass> grasslist(N);

    

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
    //init
    // world_size < 3 use serial
    // world_size >= 3 use one processor for one agent
    initialize_parallel(grasslist, sheeplist, wolflist);

    if (world_rank == 0){
        for (int t = 0; t < T; t++)
        {
    #ifdef visualization
            save2mat(sheep, sheeplist, t);
            save2mat(wolve, wolflist, t);
            save2matInt(grass, grasslist, t);
    #endif
            sendTime(t);
            animalNum[0 + t * 3] = tot_sheep;
            animalNum[1 + t * 3] = tot_wolve;
            animalNum[2 + t * 3] = tot_grass;
            animalNumVec[0 + t * 3] = tot_sheep;
            animalNumVec[1 + t * 3] = tot_wolve;
            animalNumVec[2 + t * 3] = tot_grass;
            act_master(sheeplist, wolflist, grasslist, t);
            // printf("t %d, s %d, w %d, g %d\n",t, tot_sheep, tot_wolve, tot_grass);
        }
        MPI_Finalize();
        
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
    else{
        bool flag = true;
        while(flag == true){
            flag = main_slave(world_rank);
        }
    }
#ifdef visualization
	mat2hdf5(sheep, wolve, grass, animalNum, setting, filename);
#endif
    


#ifdef visualization
    delete [] sheep;
    delete [] grass;
    delete [] wolve;
#endif
    delete [] animalNum;
    delete [] setting;

    return 0;
}

