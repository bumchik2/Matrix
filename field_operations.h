#pragma once


#include <vector>

using std::vector;

template<typename Field>
class FieldOperations {
private:
	static vector<vector<Field>> multiplicateTwoDegreeSquares (
			const vector<vector<Field>>& square1, const vector<vector<Field>>& square2);

	static vector<vector<vector<Field>>> _split (const vector<vector<Field>>& square);

	static vector<vector<Field>> _getSquare (const vector<vector<Field>>& big_square,
				unsigned raw_start, unsigned column_start, unsigned square_size);

	static vector<vector<Field>> _sum (const vector<vector<Field>>& square_1,
				const vector<vector<Field>>& square_2);

	static vector<vector<Field>> _diff (const vector<vector<Field>>& square_1,
				vector<vector<Field>> square_2);

	template <unsigned N0, unsigned M0, typename Field0>
	friend class Matrix;
};


template<typename Field>
vector<vector<Field>> FieldOperations<Field>::multiplicateTwoDegreeSquares (
		const vector<vector<Field>>& square1, const vector<vector<Field>>& square2) {
	if (square1.size() == 1) {
		return {{square1[0][0] * square2[0][0]}};
	} else {
		vector<vector<vector<Field>>> A = _split(square1);
		vector<vector<vector<Field>>> B = _split(square2);
		vector<vector<vector<Field>>> P = {
			multiplicateTwoDegreeSquares(_sum(A[0], A[3]), _sum(B[0], B[3])),
			multiplicateTwoDegreeSquares(_sum(A[2], A[3]), B[0]),
			multiplicateTwoDegreeSquares(A[0], _diff(B[1], B[3])),
			multiplicateTwoDegreeSquares(A[3], _diff(B[2], B[0])),
			multiplicateTwoDegreeSquares(_sum(A[0], A[1]), B[3]),
			multiplicateTwoDegreeSquares(_diff(A[2], A[0]), _sum(B[0], B[1])),
			multiplicateTwoDegreeSquares(_diff(A[1], A[3]), _sum(B[2], B[3]))
		};
		vector<vector<vector<Field>>> C = {
				_sum(P[0], _sum(P[3], _diff(P[6], P[4]))),
				_sum(P[2], P[4]),
				_sum(P[1], P[3]),
				_sum(P[0], _sum(P[2], _diff(P[5], P[1])))
		};

		vector<vector<Field>> result;
		for (int i = 0; i < static_cast<int>(square1.size()); ++i) {
			result.push_back(vector<Field>(square1.size()));
		}
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				for (int raw_number = 0; raw_number < static_cast<int>(square1.size()/ 2); ++raw_number) {
					for (int column_number= 0; column_number < static_cast<int>(square1.size() / 2); ++column_number) {
						result[i * square1.size() / 2 + raw_number][j * square1.size() / 2 + column_number] =
								C[i * 2 + j][raw_number][column_number];
					}
				}
			}
		}
		return result;
	}
}


template<typename Field>
vector<vector<vector<Field>>> FieldOperations<Field>::_split (const vector<vector<Field>>& square) {
	vector<vector<vector<Field>>> result;
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			unsigned square_size = square.size() / 2;
			result.push_back(_getSquare(square, i * square_size, j * square_size, square_size));
		}
	}
	return result;
}


template<typename Field>
vector<vector<Field>> FieldOperations<Field>::_getSquare (const vector<vector<Field>>& big_square,
			unsigned raw_start, unsigned column_start, unsigned square_size) {
	vector<vector<Field>> result;
	for (unsigned i = raw_start; i < raw_start + square_size; ++i) {
		vector<Field> square_raw;
		for (unsigned j = column_start; j < column_start + square_size; ++j) {
			square_raw.push_back(big_square[i][j]);
		}
		result.push_back(square_raw);
	}
	return result;
}


template<typename Field>
vector<vector<Field>> FieldOperations<Field>::_sum (const vector<vector<Field>>& square_1,
			const vector<vector<Field>>& square_2) {
	vector<vector<Field>> result = square_1;
	for (int i = 0; i < static_cast<int>(square_1.size()); ++i) {
		for (int j = 0; j < static_cast<int>(square_1.size()); ++j) {
			result[i][j] += square_2[i][j];
		}
	}
	return result;
}


template<typename Field>
vector<vector<Field>> FieldOperations<Field>::_diff (const vector<vector<Field>>& square_1,
			vector<vector<Field>> square_2) {
	vector<vector<Field>> result = square_1;
	for (int i = 0; i < static_cast<int>(square_1.size()); ++i) {
		for (int j = 0; j < static_cast<int>(square_1.size()); ++j) {
			result[i][j] -= square_2[i][j];
		}
	}
	return result;
}
