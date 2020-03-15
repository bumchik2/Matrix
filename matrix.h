#pragma once


#include "field_operations.h"


#include <vector>
#include <algorithm>
#include <iostream>


template<typename T>
class CompositeFiniteChecker {
public:
	static const bool is_composite_finite = false;
};

template<>
class CompositeFiniteChecker <HelperClass<false>> {
};

template<unsigned K>
class CompositeFiniteChecker <Finite<K>> {
public:
	static const bool is_composite_finite = CompositeFiniteChecker<HelperClass<Finite<K>::_prime>>::is_composite_finite;
};


unsigned upperTwoPower (unsigned n) {
	unsigned result = 1;
	while (result < n) {
		result *= 2;
	}
	return result;
}


template<typename T>
std::vector<std::vector<T>> rectangle (unsigned n, unsigned m) {
	std::vector<std::vector<T>> result;
	for (int i = 0; i < static_cast<int>(n); ++i) {
		std::vector<T> raw;
		for (int j = 0; j < static_cast<int>(m); ++j) {
			raw.push_back(T());
		}
		result.push_back(raw);
	}
	return result;
}


template<unsigned N, unsigned M, typename Field = Rational>
class Matrix {

public:

	Matrix() : _elements (rectangle<Field>(N, M)) { }

	~Matrix ();

	const std::vector<Field>& operator [] (unsigned index) const;
	std::vector<Field>& operator [] (unsigned index);

	Matrix& operator += (const Matrix& m);
	Matrix& operator -= (const Matrix& m);

	template <unsigned K>
	Matrix<N, K, Field>& operator *= (const Matrix<M, K, Field>& m) {
		int big_square_size = std::max(upperTwoPower(N), std::max(upperTwoPower(M), upperTwoPower(K)));
		std::vector<std::vector<Field>> square1 = _toSquare(big_square_size);
		std::vector<std::vector<Field>> square2 = m._toSquare(big_square_size);
		return _fit(FieldOperations<Field>::multiplicateTwoDegreeSquares(square1, square2));
	}

	template <unsigned K>
	Matrix<N, K, Field> operator * (const Matrix<M, K, Field>& m) const {
		Matrix<N, K, Field> result = *this;
		return result *= m;
	}

	Matrix& operator *= (const Field& f);

	Matrix<M, N, Field> transposed() const;

	std::vector<Field> getRow (unsigned index) const;
	std::vector<Field> getColumn (unsigned index) const;

	unsigned rank () const;

	const static unsigned row_number = N;
	const static unsigned column_number = M;

protected:

	std::vector<std::vector<Field>> _elements;
	static const bool _composite_finite_banned = CompositeFiniteChecker<Field>::is_composite_finite;

	std::vector<std::vector<Field>> _toSquare (unsigned square_size) const;
	Matrix& _fit(const std::vector<std::vector<Field>>& rectangle);

	void _swapRaws(unsigned raw1, unsigned raw2);
	void _swapColumns(unsigned column1, unsigned column2);

	Matrix _toDiagonal() const;
	void _makeZeroColumn(unsigned i, unsigned j);
	void _moveNumbersLeft();

	bool _sign_changes = false;
	void _updateSign();

	template<unsigned N0, unsigned M0, typename Field0>
	friend class Matrix;

};


template<unsigned N, unsigned M, typename Field>
const std::vector<Field>& Matrix<N, M, Field>::operator [] (unsigned index) const {
	return _elements[index];
}


template<unsigned N, unsigned M, typename Field>
std::vector<Field>& Matrix<N, M, Field>::operator [] (unsigned index) {
	return _elements[index];
}


template<unsigned N, unsigned M, typename Field>
std::ostream& operator << (std::ostream& os, Matrix<N, M, Field> m) {
	for (int i = 0; i < static_cast<int>(N); ++i) {
		for (int j = 0; j < static_cast<int>(M); ++j) {
			os << m[i][j] << ' ';
		}
		os << '\n';
	}
	return os;
}


template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>::~Matrix<N, M, Field> () {
	// to count compositeFiniteBanned and not get warnings
	bool foo = _composite_finite_banned;
	foo ^= true;
}


template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator += (const Matrix<N, M, Field>& m) {
	_elements = FieldOperations<Field>::_sum(_elements, m._elements);
	return *this;
}


template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator + (const Matrix<N, M, Field>& m1, const Matrix<N, M, Field>& m2) {
	Matrix<N, M, Field> result = m1;
	return result += m2;
}


template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator -= (const Matrix<N, M, Field>& m) {
	_elements = FieldOperations<Field>::_diff(_elements, m._elements);
	return *this;
}


template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator - (const Matrix<N, M, Field>& m1, const Matrix<N, M, Field>& m2) {
	Matrix<N, M, Field> result = m1;
	return result -= m2;
}


template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::operator *= (const Field& f) {
	for (int i = 0; i < static_cast<int>(N); ++i) {
		for (int j = 0; j < static_cast<int>(M); ++j) {
			_elements[i][j] *= f;
		}
	}
	return *this;
}


template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> operator * (const Matrix<N, M, Field>& m, const Field& f) {
	Matrix<N, M, Field> result = m;
	return result *= f;
}


template<unsigned N, unsigned M, typename Field>
Matrix<M, N, Field> Matrix<N, M, Field>::transposed() const {
	Matrix<M, N, Field> result;
	for (int i = 0; i < static_cast<int>(M); ++i) {
		for (int j = 0; j < static_cast<int>(N); ++j) {
			result._elements[i][j] = _elements[j][i];
		}
	}
	return result;
}


template<unsigned N, unsigned M, typename Field>
std::vector<Field> Matrix<N, M, Field>::getRow(unsigned index) const {
	return _elements.at(index);
}


template<unsigned N, unsigned M, typename Field>
std::vector<Field> Matrix<N, M, Field>::getColumn(unsigned index) const {
	std::vector<Field> result;
	for (int i = 0; i < static_cast<int>(N); ++i) {
		result.push_back(_elements[i][index]);
	}
	return result;
}


template<unsigned N, unsigned M, typename Field>
std::vector<std::vector<Field>> Matrix<N, M, Field>::_toSquare (unsigned square_size) const {
	std::vector<std::vector<Field>> result = rectangle<Field> (square_size, square_size);
	for (int i = 0; i < static_cast<int>(N); ++i) {
		for (int j = 0; j < static_cast<int>(M); ++j) {
			result[i][j] = _elements[i][j];
		}
	}
	return result;
}


template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field>& Matrix<N, M, Field>::_fit (const std::vector<std::vector<Field>>& rectangle) {
	for (int i = 0; i < static_cast<int>(N); ++i) {
		for (int j = 0; j < static_cast<int>(M); ++j) {
			_elements[i][j] = rectangle[i][j];
		}
	}
	return *this;
}


template<unsigned N, unsigned M, typename Field>
void Matrix<N, M, Field>::_makeZeroColumn(unsigned raw_number, unsigned column_number) {
	for (int i = 0; i < static_cast<int>(N); ++i) {
		if (i == static_cast<int>(raw_number)) {
			continue;
		}
		Field multiplier = _elements[i][column_number] / _elements[raw_number][column_number];

		for (int j = 0; j < static_cast<int>(N); ++j) {
			_elements[i][j] -= multiplier * _elements[raw_number][j];
		}
	}
}


template<unsigned N, unsigned M, typename Field>
void Matrix<N, M, Field>::_updateSign () {
	_sign_changes ^= true;
}


template<unsigned N, unsigned M, typename Field>
void Matrix<N, M, Field>::_swapRaws (unsigned raw1, unsigned raw2) {
	std::swap(_elements[raw1], _elements[raw2]);
	_updateSign();
}


template<unsigned N, unsigned M, typename Field>
Matrix<N, M, Field> Matrix<N, M, Field>::_toDiagonal() const {
	Matrix<N, M, Field> result = *this;
	int upper_raw_number = 0;
	for (int column_number = 0; column_number < static_cast<int>(M); ++column_number) {
		for (int raw_number = upper_raw_number; raw_number < static_cast<int>(N); ++raw_number) {
			if (result._elements[raw_number][column_number] != static_cast<Field>(0)) {
				if (upper_raw_number != raw_number) {
					result._swapRaws(upper_raw_number, raw_number);
				}
				result._makeZeroColumn(upper_raw_number, column_number);
				++upper_raw_number;
				break;
			}
		}
	}
	result._moveNumbersLeft();
	return result;
}


template<unsigned N, unsigned M, typename Field>
void Matrix<N, M, Field>::_swapColumns (unsigned column1, unsigned column2) {
	for (int raw_number = 0; raw_number < static_cast<int>(N); ++raw_number) {
		std::swap(_elements[raw_number][column1], _elements[raw_number][column2]);
	}
	_updateSign();
}


template<unsigned N, unsigned M, typename Field>
void Matrix<N, M, Field>::_moveNumbersLeft() {
	int current_left_number = 0;
	for (int column_number = 0; column_number < static_cast<int>(N); ++column_number) {
		if (_elements[current_left_number][column_number] != static_cast<Field>(0)) {
			if (current_left_number != column_number) {
				_swapColumns(current_left_number, column_number);
			}
			++current_left_number;
		}
	}
}


template<unsigned N, unsigned M, typename Field>
unsigned Matrix <N, M, Field>::rank() const {
	Matrix<N, M, Field> diagonal_matrix = _toDiagonal();
	unsigned result = 0;
	for (int i = 0; i < static_cast<int>(std::min(N, M)); ++i) {
		if (diagonal_matrix._elements[i][i] != static_cast<Field>(0)) {
			++result;
		} else {
			break;
		}
	}
	return result;
}
