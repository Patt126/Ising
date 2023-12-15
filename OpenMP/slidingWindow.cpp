#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include<fstream>
#include<ctime>
#include <chrono>
#include <omp.h>
#include <random>



#define L 200
#define N (L*L)
#define A 50//Block side lenght
#define J 1.00
#define IT 6*1e8//number of iterations
#define CHUNKSIZE 1e5
#define TSWITCH 1e7 // TSWITCH>CHUNKSIZE*NUMBLOCKS else no sense


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




void translateMatrix( std::vector<int>& inputMatrix,std::vector<int>& helpMatrix) {
//#pragma omp parallel for num_threads(16) schedule(static, (int)CHUNKSIZE)

        for (int i = 0; i < N; ++i) {
            // Check if the new indices are within bounds
            if ((i + L) < N) { //we are not moving last row
                if ((i + 1) % L != 0) //We are not in last column
                    helpMatrix[i + L + 1] = inputMatrix[i];
                else { //if we are in last column but not last row put element at the beginning of new row
                    helpMatrix[i + 1] = inputMatrix[i];
                }
            } else if ((i + 1) % L != 0) //We are not in last column of last row
            {
                // If out of bounds, put the element at the beginning
                helpMatrix[i + L + 1 - N] = inputMatrix[i]; //+1 to translate col -N to move to first row
            } else {//edge
                helpMatrix[0] = inputMatrix[i];
            }
        }

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



float simulate(float T,std::vector <int> & lattice, std::vector<int> num_vect, const int NUMBLOCKS,std::vector<bool> boundary) {

    using namespace std;
    vector<float> prob(2);
    vector<int> lattice1(N);
    prob[0] = exp(-4 * J / T);
    prob[1] = exp(-8 * J / T);

    int M_local,M=0;

    vector<int> t_axis(1);
    t_axis[0] = 0;
    //const int chunckSize = floor(IT/(NUMBLOCKS));

#pragma omp parallel  num_threads(NUMBLOCKS)
    {
        for (unsigned long int s = 0; s < (IT / (TSWITCH));s++)  {
            #pragma omp for schedule(static, (int)CHUNKSIZE)
            for (unsigned long int i = 0; i < TSWITCH; i++) {
                //first term is first site in our block
                // shift = column*NUMTHREAD + row*NUMBLOCKLINE*NUMTHREAD due to current htread work
                //int r = num_vect[i], c = block_vect[i];
                //int flipingSite = tStart[omp_get_thread_num()] + r * L  + c;
                //if (boundary[i]) {
                int n = num_vect[i + s*TSWITCH];
                if (n % A != 0 && (n + 1) % A != 0 && n % (A * L) != 0 && (n + L) % (A * L) != 0) {
                    flip(lattice, prob, n);
                }
            }
            #pragma omp single
            {
            translateMatrix(lattice,lattice1);
            lattice = lattice1;
        }
        }

    }
    //cout<<"COLOR: "<<color<<" ITERATION: "<<iteration;




    //float m = (float)M/N;
    return M;
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


int main() {
    using namespace std;
    unsigned seed = time(0);
    srand(seed);
    vector<int> lattice(N);
    vector<int> tStart;
    vector<bool> boundary;
    const int NUMBLOCKS = setBlockSize(A,tStart);
    cout << NUMBLOCKS << endl;
    float energy = 0;
    int M = 0;
    initialize_lattice(lattice, energy, M);
    print_lattice(lattice);
    float T = 0.1;
    vector<float> results(1);
    results[0] = 1;

    vector<int> n_vect;
    create_rand_vect(n_vect,NUMBLOCKS,tStart,boundary);



    auto start = std::chrono::high_resolution_clock::now();  // Start timing before simulation
    while (T <= 0.2) {

        results.push_back(simulate(T, lattice, n_vect, NUMBLOCKS,boundary));
        T += 0.1;

        print_lattice(lattice);
        cout << abs(results.back()) << endl;
    }
    auto end = std::chrono::high_resolution_clock::now();  // End timing after simulation
    std::chrono::duration<double> elapsed = end - start;  // Calculate elapsed time
    cout << elapsed.count() << endl;


    return 0;

}


