#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#define N 1000000
#define T 10
#define initSheepNum 400000
#define sheepGainFromFood 4
#define sheepReproduce 4 //%
#define initWolveNum 200000
#define wolveGainFromFood 20
#define wolveReproduce 5 //%
#define Grass 1
#define initGrass 500000
#define grassRegrowth 30 //time
#include <algorithm>
#define merror 0.00001
#define reproduceThreshold 2.0
#define debug 0
// #define visualization 1
#include <vector>
#include <unordered_set>
#include <mpi.h>
#include <omp.h>


// #include "hdf5.h"
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
    // omp_lock_t writelock;
    // omp_init_lock(&writelock);
    // omp_set_lock(&writelock);
    int out= 0;
    if (grasses[sheep.x() + half * sheep.y()].gNum() == 1)
    {
        sheep.addEnergy();
        grasses[sheep.x() + half * sheep.y()].reset();
        out = 1;
        // tot_grass--;
    }
    return out;
    // omp_unset_lock(&writelock);
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
    printf("s (%d, %d), w (%d, %d), g (%d, %d)\n", tot_sheep, tot_sheep_acc, tot_wolve, tot_wolf_acc, tot_grass, tot_grass_acc);
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
void renew_vector(vector<Animal> &mat){
    std::vector<Animal> newMat;
    for (Animal i: mat){
        if (i.gEnergy() > merror){
            newMat.push_back(i);
        }
    }
    mat = newMat;
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
void gen_agent_rankid(int dim, int* idlist){
    for (int j= 0; j < dim; j++){
            idlist[j] = j*3;
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
            }
        }
        if (world_rank == 1){
            init_sheep_wolve(sheeplist, 0);
        }
        if (world_rank == 2){
            init_sheep_wolve(wolflist, 1);            
        }
    }
}
template <class S>
void my_scatter_send(S* data, int* count, MPI_Datatype datatype, int* idlist, int numSub, MPI_Comm communicator){
    int world_rank;
    MPI_Comm_rank(communicator, &world_rank);
    int world_size;
    MPI_Comm_size(communicator, &world_size);
    int root = 0;

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
void my_scatter_rec(S** data, int* count, MPI_Datatype datatype, int* idlist, int numSub, MPI_Comm communicator){
    int root = 0;
    // If we are a receiver process, receive the data from the root
    int number_amount;
    MPI_Status status;
    MPI_Probe(0, 0, MPI_COMM_WORLD, &status);
    MPI_Get_count(&status, MPI_INT, &number_amount);
    *data = new S [sizeof(S) * number_amount];
    *count = number_amount;
    MPI_Recv(*data, number_amount, datatype, root, 0, communicator, MPI_STATUS_IGNORE);
}

template <class S>
void my_gather(S* data, int* count, MPI_Datatype datatype, int* idlist, int numSub, MPI_Comm communicator){
    int world_rank;
    MPI_Comm_rank(communicator, &world_rank);
    int world_size;
    MPI_Comm_size(communicator, &world_size);
    int root = 0;

    if (world_rank == root) {
        // If we are the root process, receive our data from everyone
        int i;
        for (i = 0; i < numSub; i++) {
            if (idlist[i] != world_rank) {
                MPI_Recv(data+i*(*count), (*count), datatype, idlist[i], 0, communicator, MPI_STATUS_IGNORE);
            }
        }
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
    my_scatter_send(xyn, &num_elements_per_proc3, MPI_INT, idlist, dim[0], communicator);

    // run grass patch for root's part
    ask_patch(grasslist, 0, num_elements_per_proc);
    // corner case taken cared by root proc for the agent in the far back
    if (num_elements_per_proc*dim[0] != size){
        ask_patch(grasslist, num_elements_per_proc*dim[0], size);
    }
    // receive data from slaves
    my_gather(xyn, &num_elements_per_proc3, MPI_INT, idlist, dim[0], MPI_COMM_WORLD);
    mat2grassvec(grasslist, xyn, num_elements_per_proc, num_elements_per_proc*dim[0]);
    delete [] xyn;
}
void act_grass_slaves(int* dim, int* idlist){
    // receive from root
    int* sub_xyn;
    int num_elements_per_proc;
    
    my_scatter_rec(&sub_xyn, &num_elements_per_proc, MPI_INT, idlist, dim[0], MPI_COMM_WORLD);
    // create sub grasslist in slave processor
    vector<Grassclass> sub_grasslist;
    num_elements_per_proc /= 3;
    // printf("%d, %d\n",num_elements_per_proc, (int)N);
    mat2grassvec(sub_grasslist, sub_xyn, 0, num_elements_per_proc);
    // run ask patch
    ask_patch(sub_grasslist, 0, num_elements_per_proc);
    // slave proc (sub_xyn) to vector and send back to root
    grassvec2mat(sub_grasslist, sub_xyn, 0, num_elements_per_proc);
    num_elements_per_proc *= 3;
    my_gather(sub_xyn, &num_elements_per_proc, MPI_INT, idlist, dim[0], MPI_COMM_WORLD);
    delete [] sub_xyn;
}

void act_master(vector<Animal> &sheeplist, vector<Animal> &wolflist, vector<Grassclass> &grasslist){
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int* dim = new int [3 * sizeof(int)]();
    gen_agent_processor(dim);
    int* idlist = new int [sizeof(int) * dim[0]];

    if (world_size < 3){
        // serial
        // //ask grass
        if (Grass != 0)
            ask_patch(grasslist, 0, N);
        //  ask sheep
        ask_sheep(sheeplist, grasslist);
        // //ask wolf
        ask_wolf(wolflist, sheeplist);
        get_state(sheeplist, wolflist, grasslist);
        renew_vector(sheeplist);
        renew_vector(wolflist);
    }
    else{
        // parallel
        // grass master
        if (Grass != 0)
        {
            gen_agent_rankid(dim[0], idlist);
            act_grass_master(grasslist, dim, idlist, MPI_COMM_WORLD);
        }
        // receive from sheep, wolf 

        // agent eat 

        // state calculate

        // send back to agent master and slaves
    }

    delete [] dim;
    delete [] idlist;
}

void stopCondition(int t){
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    for (int j = 1; j < world_size; j++){
        MPI_Send(&t, 1, MPI_INT, j, 44, MPI_COMM_WORLD);
    }
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
            // #pragma omp parallel num_threads(3)
            // {
            //     int i = omp_get_thread_num();
            //     if (i ==0){
                    save2mat(sheep, sheeplist, t);
                // }
                // if (i ==1){
                    save2mat(wolve, wolflist, t);
                // }
                // if (i==2){
                    save2matInt(grass, grasslist, t);
            //     }
            // }
    #endif
            animalNum[0 + t * 3] = tot_sheep;
            animalNum[1 + t * 3] = tot_wolve;
            animalNum[2 + t * 3] = tot_grass;
            act_master(sheeplist, wolflist, grasslist);
            printf("t %d\n",t);
            stopCondition(t);
            
        }
        MPI_Finalize();
    }else if(world_rank == 1){
        while(1){
            MPI_Status status;
            int* dim = new int [3 * sizeof(int)]();
            gen_agent_processor(dim);
            // sheep master

            // break while t is the last frame
            int tmp;
            MPI_Recv(&tmp, 1, MPI_INT, 0, 44, MPI_COMM_WORLD, &status);
            if (tmp == T-1){
                MPI_Finalize();
                break;
            }
        }
    }else if(world_rank == 2){
        while(1){
            MPI_Status status;
            int* dim = new int [3 * sizeof(int)]();
            gen_agent_processor(dim);
            // wolf master

            // break while t is the last frame
            int tmp;
            MPI_Recv(&tmp, 1, MPI_INT, 0, 44, MPI_COMM_WORLD, &status);
            if (tmp == T-1){
                MPI_Finalize();
                break;
            }
        }
    }
    else{
        if (world_rank % 3 == 0){
            while (true){
                MPI_Status status;
                int* dim = new int [3 * sizeof(int)]();
                gen_agent_processor(dim);

                // grass slaves
                int* idlist = new int [sizeof(int) * dim[0]];
                gen_agent_rankid(dim[0], idlist);
                act_grass_slaves(dim, idlist);

                // break while t is the last frame
                int tmp;
                MPI_Recv(&tmp, 1, MPI_INT, 0, 44, MPI_COMM_WORLD, &status);
                if (tmp == T-1){
                    MPI_Finalize();
                    break;
                }
            }
        }
        if (world_rank % 3 == 1){
            while (true){
                MPI_Status status;
                int* dim = new int [3 * sizeof(int)]();
                gen_agent_processor(dim);
                // sheep slaves
                int* idlist = new int [sizeof(int) * dim[1]];
                gen_agent_rankid(dim[1], idlist);

                // break while t is the last frame
                int tmp;
                MPI_Recv(&tmp, 1, MPI_INT, 0, 44, MPI_COMM_WORLD, &status);
                if (tmp == T-1){
                    MPI_Finalize();
                    break;
                }
            }
        }
        if (world_rank % 3 ==2){
            while (true){
                MPI_Status status;
                int* dim = new int [3 * sizeof(int)]();
                gen_agent_processor(dim);
                // wolf slaves
                int* idlist = new int [sizeof(int) * dim[2]];
                gen_agent_rankid(dim[2], idlist);


                // break while t is the last frame
                int tmp;
                MPI_Recv(&tmp, 1, MPI_INT, 0, 44, MPI_COMM_WORLD, &status);
                if (tmp == T-1){
                    MPI_Finalize();
                    break;
                }

            }
        }

    }
// 	// mat2hdf5(sheep, wolve, grass, animalNum, setting, filename);

#ifdef visualization
    delete [] sheep;
    delete [] grass;
    delete [] wolve;
#endif
    delete [] animalNum;
    delete [] setting;

    // MPI_Finalize();
    return 0;
}

