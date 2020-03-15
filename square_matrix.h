#pragma once


#include "matrix.h"


#include <vector>
#include <algorithm>
#include <iostream>
#include <stdexcept>


template<unsigned N, typename Field = Rational>
class SquareMatrix : public Matrix<N, N, Field> {

public:

	Field trace () const;
	Field det () const;

	SquareMatrix inverted() const;
	void invert();

private:

	using Matrix<N, N, Field>::_elements;

	SquareMatrix<N-1, Field> _addition (unsigned raw_number, unsigned column_number) const;

	template<unsigned N0, typename Field0>
	friend class SquareMatrix;

};


template<unsigned N, typename Field>
Field SquareMatrix<N, Field>::trace () const {
	Field result;
	for (int i = 0; i < static_cast<int>(N); ++i) {
		result += Matrix<N, N, Field>::_elements[i][i];
	}
	return result;
}


template<unsigned N, typename Field>
Field SquareMatrix<N, Field>::det () const {
	Matrix<N, N, Field> diagonal_matrix = Matrix<N, N>::_toDiagonal();
	Field result = static_cast<Field>(1);
	for (int i = 0; i < static_cast<int>(N); ++i) {
		result *= diagonal_matrix._elements[i][i];
	}
	if (diagonal_matrix._sign_changes == true) {
		result = -result;
	}
	return result;
}


template<unsigned N, typename Field>
SquareMatrix<N-1, Field> SquareMatrix<N, Field>::_addition (unsigned raw_number, unsigned column_number) const {
	SquareMatrix<N-1, Field> result;
	for (int i = 0; i < static_cast<int>(N); ++i) {
		if (i == static_cast<int>(raw_number)) {
			continue;
		} else {
			for (int j = 0; j < static_cast<int>(N); ++j) {
				if (j == static_cast<int>(column_number)) {
					continue;
				} else {
					int delta_raw = static_cast<int>(i > static_cast<int>(raw_number));
					int delta_column = static_cast<int>(j > static_cast<int>(column_number));
					result._elements[i - delta_raw][j - delta_column] = _elements[i][j];
				}
			}
		}
	}
	return result;
}


template<unsigned N, typename Field>
SquareMatrix<N, Field> SquareMatrix<N, Field>::inverted () const {
	SquareMatrix<N, Field> result;
	Field determinant = det();
	if (det() == static_cast<Field>(0)) {
		throw std::runtime_error("uninvertible matrix");
	}
	for (int i = 0; i < static_cast<int>(N); ++i) {
		for (int j = 0; j < static_cast<int>(N); ++j) {
			Field addition_j_i_determinant = (_addition(j, i)).det();
			result._elements[i][j] = addition_j_i_determinant / determinant;
			if ((i + j) % 2 == 1) {
				result._elements[i][j] = - result._elements[i][j];
			}
		}
	}
	return result;
}
