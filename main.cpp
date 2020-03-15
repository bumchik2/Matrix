#include "rational.h"
#include "finite.h"
#include "matrix.h"
#include "square_matrix.h"


#include <iostream>
#include <ctime>

using std::cin;
using std::cout;
using std::endl;
using std::vector;


template<typename T>
std::ostream& operator << (std::ostream& os, std::vector<T> v) {
	for (int i = 0; i < static_cast<int>(v.size()); ++i) {
		os << v[i] << ' ';
	}
	os << endl;
	return os;
}


template<unsigned N, unsigned M, typename Field>
void generate(Matrix<N, M, Field>& m) {
	for (int i = 0; i < static_cast<int>(N); ++i) {
		for (int j = 0; j < static_cast<int>(M); j++) {
			m[i][j] = rand() % 1000;
		}
	}
}


int main() {
	//srand(time(0));
	Matrix<3, 4, Finite<5>> M;
	generate(M);

	SquareMatrix<4, Finite<5>> sqM;
	generate(sqM);

	cout << M << endl <<sqM << endl;
	cout << (M *= sqM) << endl << M << endl;
	cout << M.rank();
	return 0;
}
