#pragma once


#include <iostream>
#include <vector>
#include <string>


template<typename T>
T max (const T& t1, const T& t2) {
	return t1 < t2 ? t2 : t1;
}

template<typename T>
T min (const T& t1, const T& t2) {
	return t1 < t2 ? t1 : t2;
}

class BigInteger {
private:
	std::vector<long long> digits_;
	static const long long base_ = 1000 * 1000 * 1000;
	static const int degree_ = 9; // base_ = 10 ^ degree_
	int sign_;

	void fix_ ();
	bool lessModule_ (const BigInteger& b) const;
	void countSum_ (const BigInteger& b1, const BigInteger& b2);
	void countDifference_ (const BigInteger& b1, const BigInteger& b2);
public:
	// constructors
	BigInteger();
	BigInteger(int a);
	BigInteger(long long a);
	BigInteger(unsigned long a);
	BigInteger(unsigned long long a);
	BigInteger(const std::string& s);
	BigInteger(const char* c);

	explicit operator int() const;
	explicit operator long long() const;
	explicit operator unsigned long() const;
	explicit operator unsigned long long() const;
	explicit operator bool() const;

	// arithmetic operators
	friend BigInteger operator + (const BigInteger& b1, const BigInteger& b2);
	BigInteger& operator += (const BigInteger& b);
	friend BigInteger operator - (const BigInteger& b1, const BigInteger& b2);
	BigInteger& operator -= (const BigInteger& b);
	friend BigInteger operator * (const BigInteger& b1, const BigInteger& b2);
	BigInteger& operator *= (const BigInteger& b);
	friend BigInteger operator / (const BigInteger& b1, const BigInteger& b2);
	BigInteger& operator /= (const BigInteger& b);
	friend BigInteger operator % (const BigInteger& b1, const BigInteger& b2);
	BigInteger& operator %= (const BigInteger& b);

	BigInteger operator - () const;
	BigInteger& operator ++ ();
	BigInteger& operator -- ();
	BigInteger operator ++ (int);
	BigInteger operator -- (int);

	friend bool operator == (const BigInteger& b1, const BigInteger& b2);
	friend bool operator < (const BigInteger& b1, const BigInteger& b2);
	friend bool operator > (const BigInteger& b1, const BigInteger& b2);
	friend bool operator <= (const BigInteger& b1, const BigInteger& b2);
	friend bool operator >= (const BigInteger& b1, const BigInteger& b2);
	friend bool operator != (const BigInteger& b1, const BigInteger& b2);

	std::string toString() const;

	friend class Rational;
};

BigInteger::BigInteger(){
	sign_ = 0;
}

BigInteger::BigInteger(int a):
	BigInteger(static_cast<long long>(a)) {}

BigInteger::BigInteger(long long a) {
	if (a == 0){
		sign_ = 0;
		return;
	}
	sign_ = a > 0 ? 1 : -1;
	if (a < 0) a = -a;
	while (a > 0){
		digits_.push_back(a % base_);
		a /= base_;
	}
}

BigInteger::BigInteger(unsigned long a):
	BigInteger(static_cast<unsigned long long>(a)) {}

BigInteger::BigInteger(unsigned long long a){
	if (a == 0){
		sign_ = 0;
		return;
	}
	sign_ = 1;
	while (a > 0){
		digits_.push_back(a % base_);
		a /= base_;
	}
}

long long getDigit (const std::string& s, int i1, int i2) {
	long long digit = 0;
	for (int i = i1; i < i2; ++i)
		digit = digit * 10 + (s[i] - '0');
	return digit;
}

BigInteger::BigInteger(const std::string& s){
	if (s.size() == 0 || (s.size() == 1 && s[0] == '-')){
		sign_ = 0;
		return;
	}
	sign_ = 1;
	int begin = 0;
	if (s[0] == '-'){
		begin = 1;
		sign_ = -1;
	}
	int n = static_cast<int>(s.size());
	int i = n - degree_;
	for (; i > begin; i -= degree_){
		digits_.push_back(getDigit(s, i, i + degree_));
	}
	digits_.push_back (getDigit(s, begin, i + degree_));
	fix_();
}

BigInteger::operator int() const {
	return operator long long ();
}

BigInteger::operator long long() const {
	long long result = 0;
	for (unsigned int i = digits_.size(); i != 0; ) {
		--i;
		result = result * base_ + digits_[i];
	}
	return result * sign_;
}

BigInteger::operator unsigned long() const {
	return operator unsigned long long();
}

BigInteger::operator unsigned long long() const {
	unsigned long long result = 0;
	for (unsigned int i = digits_.size(); i != 0; ){
		--i;
		result = result * base_ + digits_[i];
	}
	return result;
}

BigInteger::operator bool() const {
	return sign_ != 0;
}

void BigInteger::fix_(){
	long long delta = 0;
	for (unsigned int i = 0; i < digits_.size(); ++i){
		digits_[i] += delta;
		long long temp = digits_[i] % base_;
		if (temp < 0) {
			temp += base_;
		}
		delta = (digits_[i] - temp) / base_;
		digits_[i] = temp;
		if (i + 1 == digits_.size() && delta != 0){
			digits_.push_back(0);
		}
	}
	while (digits_.size() != 0 && digits_.back() == 0)
		digits_.pop_back();
	if (digits_.size() == 0){
		sign_ = 0;
	}
}

bool BigInteger::lessModule_ (const BigInteger& b) const {
	if (digits_.size() != b.digits_.size())
		return digits_.size() < b.digits_.size();
	if (digits_.size() == 0) return 0;
	for (int i = static_cast<int>(digits_.size()) - 1; i >= 0; --i){
		if (digits_[i] != b.digits_[i])
			return digits_[i] < b.digits_[i];
	}
	return 0;
}

inline long long elementOrZero (const std::vector<long long>& v, int index) {
	if (index >= static_cast<int>(v.size()))
		return 0;
	return v.at(index);
}

void BigInteger::countSum_ (const BigInteger& b1, const BigInteger& b2) {
	for (unsigned int i = 0; i < max (b1.digits_.size(),
			b2.digits_.size()); ++i) {
		digits_.push_back(elementOrZero (b1.digits_, i) +
				elementOrZero (b2.digits_, i));
	}
	fix_();
}

void BigInteger::countDifference_ (const BigInteger& b1, const BigInteger& b2) {
	for (unsigned int i = 0; i < max (b1.digits_.size(),
			b2.digits_.size()); ++i) {
		digits_.push_back(elementOrZero (b1.digits_, i) -
				elementOrZero (b2.digits_, i));
	}
	fix_();
}

BigInteger& BigInteger::operator += (const BigInteger& b){
	return *this = *this + b;
}

BigInteger operator + (const BigInteger& b1, const BigInteger& b2) {
	if (b1.sign_ == 0 || b2.sign_ == 0)
		return b1.sign_ == 0 ? b2 : b1;
	BigInteger result;
	if (b1.sign_ == b2.sign_){
		result.sign_ = b1.sign_;
		result.countSum_ (b1, b2);
	} else {
		if (b1.lessModule_ (b2)){
			result.sign_ = b2.sign_;
			result.countDifference_ (b2, b1);
		} else {
			result.sign_ = b1.sign_;
			result.countDifference_ (b1, b2);
		}
	}
	return result;
}

BigInteger operator - (const BigInteger& b1, const BigInteger& b2) {
	return b1 + (-b2);
}

BigInteger& BigInteger::operator -= (const BigInteger& b){
	return *this = *this + (-b);
}

BigInteger& BigInteger::operator *= (const BigInteger& b) {
	return *this = *this * b;
}

BigInteger operator * (const BigInteger& b1, const BigInteger& b2) {
	if (b1.sign_ == 0 || b2.sign_ == 0)
		return 0;
	BigInteger result;
	result.sign_ = b1.sign_ * b2.sign_;
	for (int i = b1.digits_.size() - 1; i >= 0; --i) {
		for (int j = b2.digits_.size() - 1; j >= 0; --j) {
			while (i + j + 2 > static_cast<int>(result.digits_.size()))
				result.digits_.push_back(0);
			long long newDigit = b1.digits_[i] * b2.digits_[j];
			result.digits_[i+j] += newDigit % b1.base_;
			result.digits_[i+j+1] += newDigit / b1.base_;
		}
	}
	result.fix_();
	return result;
}

BigInteger& BigInteger::operator /= (const BigInteger& b){
	return *this = *this / b;
}

BigInteger operator / (const BigInteger& b1, const BigInteger& b2) {
	BigInteger result;
	result.sign_ = b1.sign_ * b2.sign_;
	for (int i = 0; i < static_cast<int>(b1.digits_.size()); ++i)
		result.digits_.push_back (0);
	for (int i = b1.digits_.size() - 1; i >= 0; --i) {
		long long left = 0;
		long long right = result.base_ - 1;
		while (left != right) {
			long long sr = (left + right + 1) / 2;
			result.digits_[i] = sr;
			if (b1.lessModule_ (result * b2)) {
				right = sr - 1;
			} else {
				left = sr;
			}
		}
		result.digits_[i] = left;
	}
	result.fix_();
	return result;
}

BigInteger operator % (const BigInteger& b1, const BigInteger& b2) {
	return b1 - b1 / b2 * b2;
}

BigInteger& BigInteger::operator %= (const BigInteger& b){
	return *this -= *this / b * b;
}

BigInteger BigInteger::operator - () const {
	BigInteger result = *this;
	result.sign_ = -result.sign_;
	return result;
}

BigInteger& BigInteger::operator ++ (){
	return *this += 1;
}

BigInteger& BigInteger::operator -- (){
	return *this -= 1;
}

BigInteger BigInteger::operator ++ (int){
	BigInteger result = *this;
	*this += 1;
	return result;
}

BigInteger BigInteger::operator -- (int){
	BigInteger result = *this;
	*this -= 1;
	return result;
}

bool operator == (const BigInteger& b1, const BigInteger& b2){
	if (b1.sign_ != b2.sign_ || b1.digits_.size() != b2.digits_.size())
		return 0;
	for (int i = 0; i < static_cast<int>(b1.digits_.size()); ++i){
		if (b1.digits_[i] != b2.digits_[i]) return 0;
	}
	return 1;
}

bool operator < (const BigInteger& b1, const BigInteger& b2){
	if (b1.sign_ != b2.sign_)
		return b1.sign_ < b2.sign_;
	if (b1.sign_ != -1) return b1.lessModule_ (b2);
	return b2.lessModule_(b1);
}

bool operator > (const BigInteger& b1, const BigInteger& b2){
	return b2 < b1;
}

bool operator <= (const BigInteger& b1, const BigInteger& b2){
	return !(b1 > b2);
}

bool operator >= (const BigInteger& b1, const BigInteger& b2){
	return !(b1 < b2);
}

bool operator != (const BigInteger& b1, const BigInteger& b2){
	return !(b1 == b2);
}

std::string formattedDigit(const long long number, const int degree){
	std::string result;
	std::string istring = std::to_string(number);
	for (int j = 0; j < degree - static_cast<int>(istring.size()); ++j) {
		result += '0';
	}
	return result + istring;
}

std::string BigInteger::toString() const {
	if (sign_ == 0) return "0";
	std::string result;
	if (sign_ == -1)
		result += '-';
	for (int i = digits_.size() - 1; i >= 0; --i) {
		if (i + 1 == static_cast<int>(digits_.size()))
			result += std::to_string(digits_[i]);
		else
			result += formattedDigit(digits_[i], degree_);
	}
	return result;
}

std::istream& operator >> (std::istream& is, BigInteger& b){
	std::string s;
	is >> s;
	b = BigInteger(s);
	return is;
}

std::ostream& operator << (std::ostream& os, const BigInteger& b) {
	os << b.toString();
	return os;
}


// Rational part begins

BigInteger gcd (BigInteger b1, BigInteger b2);

class Rational {
public:
	Rational ():
		numerator_(0), denominator_(1) {}
	Rational (const BigInteger& b):
		numerator_ (b), denominator_(1) {}
	Rational (int i):
		numerator_ (i), denominator_(1) {}
	Rational (const std::string& s1, const std::string& s2):
		numerator_(static_cast<BigInteger>(s1)), denominator_(
			static_cast<BigInteger>(s2)) { fix_(); }
	Rational (const BigInteger& b1, const BigInteger& b2):
		numerator_ (b1), denominator_ (b2) { fix_(); }

	Rational& operator += (const Rational& r);
	friend Rational operator + (const Rational& r1, const Rational& r2);
	Rational operator - () const;
	Rational& operator -= (const Rational& r);
	friend Rational operator - (const Rational& r1, const Rational& r2);
	Rational& operator *= (const Rational& r);
	friend Rational operator * (const Rational& r1, const Rational& r2);
	Rational& operator /= (const Rational& r);
	friend Rational operator / (const Rational& r1, const Rational& r2);

	friend bool operator == (const Rational& r1, const Rational& r2);
	friend bool operator < (const Rational& r1, const Rational& r2);
	friend bool operator > (const Rational& r1, const Rational& r2);
	friend bool operator <= (const Rational& r1, const Rational& r2);
	friend bool operator >= (const Rational& r1, const Rational& r2);
	friend bool operator != (const Rational& r1, const Rational& r2);

	explicit operator double ();

	std::string toString () const;
	std::string asDecimal(size_t precision = 0) const;

private:
	static const char separator_ = '.';
	void fix_ ();
	int sign_ () const;

	BigInteger numerator_;
	BigInteger denominator_;

	template<unsigned N, unsigned M, typename Field>
	friend class Matrix;
};

int Rational::sign_ () const {
	return numerator_.sign_ * denominator_.sign_;
}

BigInteger abs (const BigInteger& b) {
	return b < 0 ? -b : b;
}

void Rational::fix_() {
	BigInteger GCD = gcd(abs(numerator_), abs(denominator_));
	numerator_ /= GCD;
	denominator_ /= GCD;
	if (numerator_.sign_ * denominator_.sign_ == -1) {
		numerator_.sign_ = -1;
		denominator_.sign_ = 1;
	} else if (numerator_.sign_ == 0) {
		denominator_.sign_ = 1;
	} else {
		numerator_.sign_ = denominator_.sign_ = 1;
	}
}

Rational& Rational::operator += (const Rational& r) {
	(numerator_ *= r.denominator_) += r.numerator_ * denominator_;
	denominator_ *= r.denominator_;
	fix_();
	return *this;
}

Rational operator + (const Rational& r1, const Rational& r2) {
	Rational result = r1;
	return result += r2;
}

Rational& Rational::operator -= (const Rational& r) {
	(numerator_ *= r.denominator_) -= (r.numerator_ * denominator_);
	denominator_ *= r.denominator_;
	fix_();
	return *this;
}

Rational operator - (const Rational& r1, const Rational& r2) {
	Rational result = r1;
	return result -= r2;
}

Rational Rational::operator - () const {
	Rational result = *this;
	result.numerator_.sign_ = -result.numerator_.sign_;
	return result;
}

Rational& Rational::operator *= (const Rational& r) {
	numerator_ *= r.numerator_;
	denominator_ *= r.denominator_;
	fix_();
	return *this;
}

Rational operator * (const Rational& r1, const Rational& r2) {
	Rational result = r1;
	return result *= r2;
}

Rational& Rational::operator /= (const Rational& r) {
	numerator_ *= r.denominator_;
	denominator_ *= r.numerator_;
	fix_();
	return *this;
}

Rational operator / (const Rational& r1, const Rational& r2) {
	Rational result = r1;
	return result /= r2;
}

bool operator == (const Rational& r1, const Rational& r2) {
	return (r1.numerator_ == r2.numerator_ &&
			r1.denominator_ == r2.denominator_);
}

bool operator < (const Rational& r1, const Rational& r2) {
	return (r1 - r2).sign_() < 0;
}

bool operator > (const Rational& r1, const Rational& r2) {
	return r2 < r1;
}

bool operator <= (const Rational& r1, const Rational& r2) {
	return !(r1 > r2);
}

bool operator >= (const Rational& r1, const Rational& r2) {
	return !(r1 < r2);
}

bool operator != (const Rational& r1, const Rational& r2) {
	return !(r1 == r2);
}

std::string Rational::toString() const {
	std::string result = numerator_.toString();
	if (denominator_ != 1) {
		result += '/';
		result += denominator_.toString();
	}
	return result;
}

std::string Rational::asDecimal(size_t precision) const {
	BigInteger resultb = 1;
	for (size_t i = 0; i < precision; ++i)
		resultb *= 10;
	resultb *= numerator_;
	resultb /= denominator_;
	std::string results = resultb.toString();
	std::string answer;
	int start = 0;
	int n = static_cast<int>(results.size());
	if (results[0] == '-') {
		answer += '-';
		start = 1;
	}
	if (n - start <= static_cast<int>(precision)) {
		answer += '0';
		answer += separator_;
		answer += std::string (precision - n + start, '0');
		answer += std::string (&results[start], results.size()-start);
	} else {
		answer += std::string (&results[start], n-precision-start);
		if (precision != 0)
			answer += separator_;
		answer += std::string (&results[n-precision], precision);
	}
	return answer;
}

double decimalDegree (int number) {
	double answer = 1.0;
	if (number >= 0)
		for (int i = 0; i < number; ++i)
			answer /= 10.0;
	else
		for (int i = 0; i < -number; ++i)
			answer *= 10.0;
	return answer;
}

Rational::operator double () {
	const int decimalPlaces = 30;
	std::string s = this->asDecimal(decimalPlaces);
	double answer = 0.0;
	int counter = decimalPlaces;
	for (unsigned int i = s.size(); i != 0;) {
		--i;
		if (s[i] == separator_ || s[i] == '-')
			continue;
		double k = decimalDegree (counter);
		answer += k * (s[i] - '0');
		--counter;
	}
	if (s[0] == '-')
		answer = -answer;
	return answer;
}

BigInteger gcd (BigInteger b1, BigInteger b2){
	if (b1) return gcd (b2 % b1, b1);
	return b2;
}

std::ostream& operator << (std::ostream& os, const Rational& r) {
	os << r.toString();
	return os;
}

