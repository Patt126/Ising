#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <iostream>


#define L 25
#define N (L*L)
#define XNN 1
#define YNN 1
#define beta 0.417
#define J 1.00



void sweep()
{
	double prob[5];
	int s[N];
	int i, k;
	int nn, sum, delta;
	float DE1, DE2;
	DE1 = 4 * J;
	DE2 = 8 * J;
	for (k = 0; k < N; k++) {
		/* Choose a site */
		i = rand() % N;
		printf("%d\n", i);
		/* Calculate the sum of the neighbouring spins */
		if ((nn = i + 1) >= N) nn -= N;

		sum = s[nn];
		if ((nn = i - 1) < 0) nn += N;
		sum += s[nn];
		if ((nn = i + 1) >= N) nn -= N;
		sum += s[nn];
		if ((nn = i - 1) < 0) nn += N;
		sum += s[nn];
		/* Calculate the change in energy */
		delta = sum * s[i];
		/* Decide whether to flip spin */
		if (delta <= 0) {
			s[i] = -s[i];
		}
		else if (rand() < prob[delta]) {
			s[i] = -s[i];
		}
	}
}

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

	if (sum <= 0) {
		matrix[r][c] = -matrix[r][c];
	}
	else if (sum == 2) {
		if (rand() / RAND_MAX < exp(-DE1 * beta)) {
			matrix[r][c] = -matrix[r][c];
		}
	}
	else {
		if (rand() / RAND_MAX < exp(-DE2 * beta)) {
			matrix[r][c] = -matrix[r][c];
		}
	}
	return sum;
}

void initialize_matrix(std::vector < std::vector<int> >& matrix) {
	int k;
	int size = L;
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			k = rand() % 2;
			if (k == 0) {
				matrix[i][j] = -1;
			}
			else {
				matrix[i][j] = -1;
			}
		}
	}

}

void print_lattice(std::vector < std::vector<int> >& matrix) {

	int i, j;
	for (i = 0; i < L; i++) {
		std::cout << std::endl;
		for (j = 0; j < L; j++) {
			if (matrix[i][j] == -1) {
				std::cout << "o" << " ";
			}
			else {
				std::cout << "x" << " ";

			}
		}
	}
	return;
}






int main() {

	using namespace std;
	vector<vector<int>> lattice(L, vector<int>(L));
	initialize_matrix(lattice);
	float prob1, prob2;
	float DE1 = 4 * J;
	float DE2 = 8 * J;
	print_lattice(lattice);
	int i;
	vector<float> energy_vec = { -L * L * J };
	vector<int> x_axis = { 0 };
	float temp =  0;
	for (i = 0; i < 1e6; i++) {
		temp = energy_vec.back();
		energy_vec.push_back(temp + flip(lattice, DE1, DE2));
		x_axis.push_back(i);

		if (i % 1000 == 0) {
			cout << energy_vec.back() << endl;
		}
	}


	cout<< endl;
	print_lattice(lattice);


	return 0;
}