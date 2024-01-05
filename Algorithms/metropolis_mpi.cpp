#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include<fstream>
#include<ctime>
#include <chrono>
#include <random>
#include <mpi.h>



#define L 100
#define N (L*L)
#define A 50//Block side lenght
#define J 1.00
#define IT 2*1e4//number of iterations
#define CHUNKSIZE 1e6


void print_lattice(std::vector<int>& lattice) {

    int i;
    for (i = 0; i < N; i++) {
        if(i%L == 0) std::cout<<std::endl;
        if (lattice[i] == -1) {
            std::cout << "o" << " ";
        } else {
            std::cout << "x" << " ";


        }
    }
    std::cout<<std::endl;
}


//if needed function to correctly evaluate energy
float evaluate(std::vector<int>& lattice) {
    int sum=0;
    for (int i=0;i<N;i++){
        //energy contibute
        if (i >= L) { //NO FIRST ROW
            sum += lattice[i - L]*lattice[i]*2;
        }
        if (i%L != 0) { //NO FIRST COLUMN
            sum += lattice[i-1]*lattice[i]*2;
        }

        if(i >= L*(L - 1)) { //LAST ROW
            sum += lattice[i-L*(L-1)]*lattice[i]*2; //times 2 two count also contribute where i = 0
        }
        if((i + 1) % L  == 0) { //LAST COLUMN
            sum += lattice[i-(L-1)]*lattice[i]*2;
        }
    }
    return -J*sum;
}

int atomicflip(std::vector <int> & lattice, std::vector<float>& prob, int site) {
    int sum = 0;
    int spin;
    if (site < L) {
        spin = lattice[site+L*(L-1)];
        sum += lattice[site+L*(L-1)];
    }
    else {
        spin = lattice[site-L];
        sum += lattice[site-L];
    }
    if (site % L == 0) {
        spin = lattice[site+(L-1)];
        sum += lattice[site + (L - 1)];
    }
    else {
        spin = lattice[site-1];
        sum += lattice[site - 1];
    }

    if (site >= L*(L - 1)) {
        spin = lattice[site-L*(L-1)];
        sum += lattice[site - L*(L-1)];
    }
    else {
        spin = lattice[site+L];
        sum += lattice[site + L];
    }
    if ((site+1) % L == 0) {
        spin = lattice[site-(L-1)];
        sum += lattice[site - (L-1)];
    }
    else {
        spin = lattice[site+1];
        sum += lattice[site + 1];
    }
    int delta = 2*sum*lattice[site];
    if (delta <= 0) {
#pragma omp atomic write
        lattice[site] = -lattice[site];
    }

    else if (delta == 4) {
        float rnd = (rand() % 10000)/1e4;
        if (rnd < prob[0] ){
#pragma omp atomic write
            lattice[site] = -lattice[site];
        }
        else{
            return 0;
        }
    }
    else if (delta==8){
        float rnd = (rand() % 10000)/1e4;
        if (rnd < prob[1]) {
#pragma omp atomic write
            lattice[site] = -lattice[site];
        }
        else{
            return 0;
        }
    }
    return 2*lattice[site];

}

int flip(std::vector <int> & lattice, std::vector<float>& prob, int site) {

    int sum = 0;

    if (site < L) {
        sum += lattice[site+L*(L-1)];
    }
    else {
        sum += lattice[site-L];
    }
    if (site % L == 0) {
        sum += lattice[site + (L - 1)];
    }
    else {
        sum += lattice[site - 1];
    }

    if (site >= L*(L - 1)) {
        sum += lattice[site - L*(L-1)];
    }
    else {
        sum += lattice[site + L];
    }
    if ((site+1) % L == 0) {
        sum += lattice[site - (L-1)];
    }
    else {
        sum += lattice[site + 1];
    }
    int delta = 2*sum*lattice[site];
    if (delta <= 0) {
        lattice[site] = -lattice[site];
    }

    else if (delta == 4) {
        float rnd = (rand() % 10000)/1e4;
        if (rnd < prob[0] ){
            lattice[site] = -lattice[site];
        }
        else{
            return 0;
        }
    }
    else if (delta==8){
        float rnd = (rand() % 10000)/1e4;
        if (rnd < prob[1]) {
            lattice[site] = -lattice[site];
        }
        else{
            return 0;
        }
    }
    return 2*lattice[site];

}

void initialize_lattice(std::vector<int>& lattice, float& energy, int& M) {
    int k;
    int size = N;
    int sum = 0;
    for (int i = 0; i < size; i++) {
        k = rand() % 2;
        if (k == 0) {
            lattice[i] = -1;
            M -= 1;
        }
        else {
            lattice[i] = 1;
            M += 1;
        }
        //energy contibute
        if (i >= L) { //NO FIRST ROW
            sum += lattice[i - L]*lattice[i]*2;
        }
        if (i%L != 0) { //NO FIRST COLUMN
            sum += lattice[i-1]*lattice[i]*2;
        }

        if(i >= L*(L - 1)) { //LAST ROW
            sum += lattice[i-L*(L-1)]*lattice[i]*2; //times 2 two count also contribute where i = 0
        }
        if((i + 1) % L  == 0) { //LAST COLUMN
            sum += lattice[i-(L-1)]*lattice[i]*2;
        }

    }
    energy+= -J*sum;
}






int write_file(std::vector<float>& energy_vec, std::vector<float>& m, std::vector<int>& t){
    std::ofstream myfile ("data.txt");
    if (myfile.is_open())
    {
        for(int i = 0; i < t.size(); i ++){
            myfile << energy_vec[i] << " "<<abs(m[i])<<' '<<t[i]<<'\n';
        }
        myfile.close();
    }
    else std::cout << "Unable to open file";
    return 0;
}



float simulate(float T, std::vector<int>& lattice, std::vector<int> num_vect, const int NUMBLOCKS, std::vector<bool> boundary, int rank, int size) {
    std::vector<float> prob(2);
    prob[0] = exp(-4 * J / T);
    prob[1] = exp(-8 * J / T);

    // Calculate the number of iterations per process
    int local_iterations = IT / size;
    int start_iter = rank * local_iterations;
    int end_iter = (rank + 1) * local_iterations;

    int M_local = 0;  // Local magnetization

    for (int iter = start_iter; iter < end_iter; ++iter) {
        int n = num_vect[iter];
        if (n % A == 0 || (n + 1) % A == 0 || n % (A * L) == 0 || (n + L) % (A * L) == 0) {
            M_local += atomicflip(lattice, prob, n);
        } else {
            M_local += flip(lattice, prob, n);
        }
    }

    // Aggregate local magnetizations to get the total magnetization
    int globalM;
    MPI_Allreduce(&M_local, &globalM, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return globalM;
}




void create_rand_vect(std::vector<int> &rand_vect,const int NUMBLOCKS, std::vector<int> & tStart,std::vector<bool>& boundary) {
    int i;

    for (unsigned long int i = 0; i < IT/CHUNKSIZE; i++) {
        //first term is first site in our block
        // shift = column*NUMTHREAD + row*NUMBLOCKLINE*NUMTHREAD due to current htread work
        for (int j = 0;j<CHUNKSIZE;j++) {
            int r = rand() % A, c = rand() % A;
            rand_vect.push_back(tStart[i%NUMBLOCKS] + r * L + c);
            if(r==0 || r==(A-1)||c==0||c==(A-1)){
                boundary.push_back(true);
            }
            else{
                boundary.push_back(false);
            }
        }
    }
}


//return the number of block in a line = num  blocks in a col
const int setBlockSize(int dimSideBlock,std::vector<int>& tStart) {
    if (L % dimSideBlock == 0) {
        const int numCellRow = L / dimSideBlock; // = numero di celle per colonna
        const int NUMBLOCKS = numCellRow*numCellRow;
        for(int i=0;i<  NUMBLOCKS;i++){
            tStart.push_back(floor(i/numCellRow)*A*L + i%numCellRow*A); //thread at which each block start
        }
        return NUMBLOCKS;
    }
    else {
        std::cout << "La dimensione del reticolo non è multiplo della dimensione del blocco";
    }
    return 0;
}


int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Initialize variables
    std::vector<int> lattice(N);
    std::vector<int> tStart;
    std::vector<bool> boundary;
    int NUMBLOCKS = setBlockSize(A, tStart);
    float energy = 0;
    int M = 0;

    // Root process initializes the lattice
    if (rank == 0) {
        unsigned seed = time(0);
        srand(seed);
        initialize_lattice(lattice, energy, M);
        print_lattice(lattice);
    }

    // Broadcast the initialized lattice to all processes
    MPI_Bcast(lattice.data(), lattice.size(), MPI_INT, 0, MPI_COMM_WORLD);

    float T = 0.1;
    std::vector<float> results(1);
    results[0] = 1;

    std::vector<int> n_vect;
    create_rand_vect(n_vect, NUMBLOCKS, tStart, boundary);

    auto start = std::chrono::high_resolution_clock::now();
    T=0.1;

    while (T < 0.2) {
        float result = simulate(T, lattice, n_vect, NUMBLOCKS, boundary, rank, size);
        if (rank == 0) {  // Only root process handles the result
            results.push_back(result);
            print_lattice(lattice);
            T += 0.1;
            std::cout << abs(results.back()) << std::endl;
        }
        MPI_Bcast(&T, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);  // Synchronize T value among all processes
    }

    auto end = std::chrono::high_resolution_clock::now();
    if (rank == 0) {  // Only root process prints the elapsed time
        std::chrono::duration<double> elapsed = end - start;
        std::cout << elapsed.count() << std::endl;
    }

    MPI_Finalize();
    return 0;
}