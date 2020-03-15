#pragma once


#include <math.h>
#include <iostream>
#include <vector>
#include <stdexcept>


template<bool>
class HelperClass {
};


template <unsigned M, unsigned prefix>
class DividerExists {
public:
	static const bool value = (M % prefix == 0 || DividerExists<M, prefix - 1>::value);
};


template<unsigned M>
class DividerExists <M, 0> {
public:
	static const bool value = true;
};


template<unsigned M>
class DividerExists <M, 1> {
public:
	static const bool value = false;
};

template<unsigned N>
class Finite {

public:

	Finite() : value(0) { }
	Finite(unsigned value) : value(value % N) { }
	Finite(const Finite& f) : value(f.value) { }

	unsigned value;

	static unsigned pow (unsigned number, unsigned degree);

	static unsigned reverse_modulo(unsigned number);

	Finite operator - () const;

	Finite& operator += (const Finite<N>& b);
	Finite& operator -= (const Finite<N>& b);
	Finite& operator *= (const Finite<N>& b);
	Finite& operator /= (const Finite<N>& b);

private:

	static const bool _prime = !DividerExists<N, static_cast<unsigned>(sqrt(N))>::value;

	template <typename T>
	static unsigned _reverse_modulo (T, unsigned number) {
		if (number % N == 0) {
			throw std::invalid_argument("zero has no modulo reverse number");
		}
		return pow (number, N - 2);
	}

	static unsigned _reverse_modulo (HelperClass<false>, unsigned number) {
	   return "a";
	}

	template<unsigned K, unsigned M, typename Field>
	friend class Matrix;

	template <typename Any>
	friend class CompositeFiniteChecker;

};


template<unsigned N>
unsigned Finite<N>::reverse_modulo(unsigned number) {
	return _reverse_modulo(HelperClass<_prime>(), number) % N;
}


template<unsigned N>
std::ostream& operator << (std::ostream& os, const Finite<N>& f) {
	os << f.value;
	return os;
}


template<unsigned N>
Finite<N> Finite<N>::operator - () const {
	return Finite<N>(N - value);
}


template<unsigned N>
Finite<N>& Finite<N>::operator += (const Finite<N>& b) {
	(value += b.value) %= N;
	return *this;
}

template<unsigned N>
Finite<N> operator + (const Finite<N>& a, const Finite<N>& b) {
	Finite<N> result = a;
	return result += b;
}


template<unsigned  N>
Finite<N>& Finite<N>::operator -= (const Finite<N>& b) {
	((value += N) -= b.value) %= N;
	return *this;
}

template<unsigned N>
Finite<N> operator - (const Finite<N>& a, const Finite<N>& b) {
	Finite<N> result = a;
	return (result -= b);
}


template<unsigned  N>
Finite<N>& Finite<N>::operator *= (const Finite<N>& b) {
	(value *= b.value) %= N;
	return *this;
}

template<unsigned N>
Finite<N> operator * (const Finite<N>& a, const Finite<N>& b) {
	Finite<N> result = a;
	return (result *= b);
}


template<unsigned  N>
Finite<N>& Finite<N>::operator /= (const Finite<N>& b) {
	(value *= reverse_modulo(b.value)) %= N;
	return *this;
}

template<unsigned N>
Finite<N> operator / (const Finite<N>& a, const Finite<N>& b) {
	Finite<N> result = a;
	return (result /= b);
}


template<unsigned N>
unsigned Finite<N>::pow (unsigned number, unsigned degree) {
	number %= N;
	if (degree == 0) {
		return number;
	}

	std::vector<unsigned> number_powers_vec;
	unsigned two_power = 1;
	unsigned deg_cnt = 1;
	number_powers_vec.push_back(1);
	number_powers_vec.push_back(number);

	while (two_power * 2 <= degree) {
		two_power *= 2;
		++deg_cnt;
		number_powers_vec.push_back(deg_cnt == 1 ?
				number : (number_powers_vec.back() * number_powers_vec.back()) % N);
	}

	unsigned answer = 1;
	while (two_power > 0) {
		if (degree >= two_power) {
			(answer *= number_powers_vec[deg_cnt]) %= N;
			degree -= two_power;
		}
		two_power /= 2;
		--deg_cnt;
	}
	return answer;
}


template<unsigned N>
bool operator == (const Finite<N>& f1, const Finite<N>& f2) {
	return f1.value == f2.value;
}


template<unsigned N>
bool operator != (const Finite<N>& f1, const Finite<N>& f2) {
	return f1.value != f2.value;
}
