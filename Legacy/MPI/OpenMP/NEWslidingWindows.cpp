#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include<fstream>
#include<ctime>
#include <chrono>
#include <omp.h>
#include <random>



constexpr int L = 100  ;
constexpr int N = (L*L);
constexpr int THREADPERSIDE = 4; //is nothing but blocks per side, 
constexpr int NUMTHREAD = THREADPERSIDE*THREADPERSIDE; //is nothing but number of block,s
constexpr int J = 1.00;




void print_lattice(std::array<int,N>& lattice) {

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
float evaluate(std::array<int,N>& lattice) {
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

void flip(std::array <int,N> & lattice, std::array<float,2>& prob, int site, int& M_loc, int& E_loc ) {

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
            return ;
        }
    }
    else if (delta==8){
        float rnd = (rand() % 10000)/1e4;
        if (rnd < prob[1]) {
            lattice[site] = -lattice[site];
        }
        else{
            return ;
        }
    }
    M_loc += 2*lattice[site];

}

void initialize_lattice(std::array<int,N>& lattice, float& energy, int& M) {
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



void simulate(std::array<float,2> prob, std::array <int,N> & lattice, std::array<int,N/NUMTHREAD> randomVect, int& M_loc ,int& E_loc , int offset) {

    
        for (unsigned long int i = 0; i < (N/ (NUMTHREAD));i++) {
                //first term is first site in our block
                // shift = column*NUMTHREAD + row*NUMBLOCKLINE*NUMTHREAD due to current htread work
                //int r = num_vect[i], c = block_vect[i];
                //int flipingSite = tStart[omp_get_thread_num()] + r * L  + c;
                //if (boundary[i]) {
                int n = randomVect[i];
                if (n != -1) {
                    flip(lattice, prob, n + offset, M_loc, E_loc);
                }
            }
        

    
}



void create_rand_vect(std::array<int,N/NUMTHREAD> &rand_vect, std::array<int,NUMTHREAD> & tStart, const int A) {
    //  N/NUMBLOCKS is the number of iteration for each task
    for (int j = 0;j<N/NUMTHREAD;j++) {
                int r = rand() % A, c = rand() % A;
                if(r%A != 0 && c%A!=0 ){ //not last or first row or column
                    rand_vect[j] = (r * L + c);
                }
                else{
                    rand_vect[j] = -1; //will mean do nothing
                }
            }
    }



//return the number of block in a line = num  blocks in a col
const int setBlockSize(std::array<int,NUMTHREAD>& tStart) {
    
    if (L % THREADPERSIDE == 0) {
        const int A = L / THREADPERSIDE; // = larghezze di un blocco
    for(int i=0;i<  NUMTHREAD;i++){
            tStart[i] = (floor(i/THREADPERSIDE)*A*L + i%THREADPERSIDE*A); //thread at which each block start
        }
        return A;
    }
    else {
        std::cout << "It's not possible to fill a line of lenght: "<<L<<" with: "<<THREADPERSIDE<<" blocks"<<std::endl;
    }
    return 0;
}

//Evaluate exact equilibrium value of per site magnetization from Onsager's formula for 2D case
auto Mexact =[](int T){return N*std::pow((1.0 - std::pow(std::sinh(2 * J / T), -4)), 1.0 / 8.0);} ;


int main() {
    unsigned seed = time(0);
    srand(seed);
    std::array<int,N> lattice;
    std::array<int,NUMTHREAD> tStart; //starting point of each block
    const int A = setBlockSize(tStart);

    float energy = 0;
    int M = 0;
    double mExact = 1;
    initialize_lattice(lattice, energy, M);
    print_lattice(lattice);
    float T = 0.1;
  
    int M_loc = 0;
    int E_loc = 0;

    const double tollerance = std::pow(10,-5);
    int step = 0;
    auto start = std::chrono::high_resolution_clock::now();  // Start timing before simulation
    int offset = 0;
    std::array<float,2> prob;
    std::array<int,N/NUMTHREAD> random;
    
    #pragma omp parallel num_threads(NUMTHREAD) shared(M) 
    {
        #pragma omp single nowait
        {
            while (T <= 0.2) {
                
                prob[0] = exp(-4 * J / T);
                prob[1] = exp(-8 * J / T);
                step = 0;
                mExact = Mexact(T);
                std::cout<<mExact<<std::endl;
                
                  while(M-mExact>tollerance){

                    create_rand_vect(random,tStart,A);

                    for(int taskNum = 0; taskNum < NUMTHREAD; taskNum++ ){
                        #pragma omp task  private(M_loc, E_loc)
                        {
                            M_loc = 0;
                            simulate(prob,lattice, random, M_loc, E_loc, tStart[taskNum]);
                            #pragma omp atomic update
                            M += M_loc;
                        }
                    }
                    #pragma omp taskwait
                    step++;
                  }  
                  T+=0.1;
                }
        }        
    }
    auto end = std::chrono::high_resolution_clock::now();  // End timing after simulation
    std::chrono::duration<double> elapsed = end - start;  // Calculate elapsed time
    std::cout << elapsed.count() << std::endl;
    print_lattice(lattice);



    return 0;

}

