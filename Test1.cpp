#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <iostream>


#define L 100
#define N (L*L)
#define T 2.4
#define beta (1/T)
#define J 1.00


int flip(std::vector < std::vector<int> >& matrix, float DE1, float DE2) {

	int r, c;
	r = rand() % L;
	c = rand() % L;
	int sum = 0;

	if (r == 0) {
		sum += matrix[L - 1][c];
	}
	else {
		sum += matrix[r - 1][c];
	}
	if (c == 0) {
		sum += matrix[r][L - 1];
	}
	else {
		sum += matrix[r][c - 1];
	}

	if (r == L - 1) {
		sum += matrix[0][c];
	}
	else {
		sum += matrix[r + 1][c];
	}
	if (c == L - 1) {
		sum += matrix[r][0];
	}
	else {
		sum += matrix[r][c + 1];
	}
    int delta = 2*J*sum*matrix[r][c];
	if (delta <= 0) {
        //std::cout<<"neg ";
		matrix[r][c] = -matrix[r][c];
	}
	else if (delta == 4*J) {
        float rnd = (rand() % 10000)/1e4;
		if (rnd < exp(-DE1 * beta)) {
			matrix[r][c] = -matrix[r][c];
		}
        else{
            return 0;
        }
	}
	else if (delta==8*J){
        float rnd = (rand() % 10000)/1e4;
		if (rnd < exp(-DE2 * beta)) {
			matrix[r][c] = -matrix[r][c];
		}
        else{
            return 0;
        }
	}

	return delta;
}

float initialize_lattice(std::vector < std::vector<int> >& matrix) {
	int k,m=0;
    int energy = 0;
	int size = L;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			k = rand() % 2;
			if (k == 0) {
				matrix[i][j] = -1;
			}
			else {
				matrix[i][j] = 1;
			}
            //energy contibute
            if (i != 0) {
                energy += matrix[i - 1][j]*matrix[i][j]*2;
            }
            if (j != 0) {
                energy += matrix[i][j-1]*matrix[i][j]*2;
            }

            if(i == L - 1) {
                energy += matrix[0][j]*matrix[i][j]*2; //times 2 two count also contribute where i = 0
            }
            if(j == L - 1) {
                energy += matrix[i][0]*matrix[i][j]*2;
            }
            }
		}
    std::cout<<std::endl<<energy<<' ';
    return -J*energy;
	}



void print_lattice(std::vector < std::vector<int> >& matrix) {

    int i, j;
    for (i = 0; i < L; i++) {
        std::cout << std::endl;
        for (j = 0; j < L; j++) {
            if (matrix[i][j] == -1) {
                std::cout << "o" << " ";
            } else {
                std::cout << "x" << " ";

            }
        }
    }
}






int main() {

	using namespace std;
	vector<vector<int> > lattice(L, vector<int>(L));
	float energy = initialize_lattice(lattice);
	float prob1, prob2;
	float DE1 = 4 * J;
	float DE2 = 8 * J;
	vector<float> energy_vec(1);
    energy_vec[0] =  energy;
	vector<int> t_axis(1);
    t_axis[0] = 0;

    cout << energy << endl;
	for (int i = 0; i < 6*1e7; i++){
		energy += flip(lattice, DE1, DE2);
        if (i % N == 0) {
            // aggiornare l'energia ogni N step (definizione di tempo nel metropolis come da libro)
           cout << energy << endl;
            energy_vec.push_back(energy_vec.back() + energy);
            t_axis.push_back(i/N);

		}
	}





	return 0;
}
