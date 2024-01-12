#include "SwendsenWangParallel.h"
#include <cmath>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <chrono>
#include <unordered_set>
#include <queue>
#include <map>
#include <numeric>
#include <omp.h>
#include <random>
// Global vector of thread-local RNGs


int main(){
    int L ;
    int L_MAX = 4096;   


    for(L = 64; L <= L_MAX; L *= 2){
        int IT = 30;
         std::cout<<"Simulation start for L = "<<L<<std::endl;
        SwendsenWangParallel simulation(1.0f, L, 0.1f, 1.0f, 0.1f, IT);
        simulation.simulate_phase_transition();
        std::cout<<"Simulation done for L = "<<L<<std::endl;
        simulation.store_results_to_file();
    }
    
    }
