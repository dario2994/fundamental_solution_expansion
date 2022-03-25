/* Implementation of the rational number for a generic ring. */

#include <map>
#include <iostream>
#include <vector>
#include <cassert>
#include <functional>
#include <numeric>
using namespace std;
typedef long long LL;


template <typename Integer>
struct Rational {
	Integer num;
	Integer den;
	Rational<Integer>(): num(0), den(1) {}
	Rational<Integer>(Integer n): num(n), den(1) {}
	Rational<Integer>(Integer n, Integer d): num(n), den(d) { normalize(); }

	void normalize() {
		assert(den != 0);
		if (den < 0) num = -num, den = -den;
		Integer g = absolute_value(gcd(num, den));
		num /= g, den /= g;
	}
};

template <typename Integer>
ostream& operator <<(ostream& out, const Rational<Integer>& x) {
	out << x.num;
	if (x.den != Integer(1)) out << "/" << x.den;
	return out;
}

template <typename Integer>
bool is_zero(const Rational<Integer>& x) { return x.num == 0; }

template <typename Integer>
bool operator ==(const Rational<Integer>& x, const Rational<Integer>& y) {
	return x.num == y.num and x.den == y.den;
}

template <typename Integer>
bool operator !=(const Rational<Integer>& x, const Rational<Integer>& y) {
	return x.num != y.num or x.den != y.den;
}

template <typename Integer>
bool operator <(const Rational<Integer>& x, const Rational<Integer>& y) {
	return x.num * y.den < x.den * y.num;
}

template <typename Integer>
bool operator >(const Rational<Integer>& x, const Rational<Integer>& y) {
	return x.num * y.den > x.den * y.num;
}

template <typename Integer>
Rational<Integer> abs(const Rational<Integer>& x) { return {absolute_value(x.num), x.den}; }

template <typename Integer>
void operator +=(Rational<Integer>& x, const Rational<Integer>& y) {
	x.num = x.num * y.den + y.num * x.den;
	x.den *= y.den;
	x.normalize();
}

template <typename Integer>
Rational<Integer> operator +(const Rational<Integer>& x, const Rational<Integer>& y) {
	return {x.num * y.den + y.num * x.den, x.den * y.den};
}

template <typename Integer>
void operator -=(Rational<Integer>& x, const Rational<Integer>& y) {
	x.num = x.num * y.den - y.num * x.den;
	x.den *= y.den;
	x.normalize();
}

template <typename Integer>
Rational<Integer> operator -(const Rational<Integer>& x, const Rational<Integer>& y) {
	return {x.num * y.den - y.num * x.den, x.den * y.den};
}

template <typename Integer>
Rational<Integer> operator -(const Rational<Integer>& x) {
	return {-x.num, x.den};
}

template <typename Integer>
void operator *=(Rational<Integer>& x, const Rational<Integer>& y) {
	x.num *= y.num;
	x.den *= y.den;
	x.normalize();
}

template <typename Integer>
Rational<Integer> operator *(const Rational<Integer>& x, const Rational<Integer>& y) {
	return {x.num * y.num, x.den * y.den};
}

template <typename Integer>
void operator /=(Rational<Integer>& x, const Rational<Integer>& y) {
	x.num *= y.den;
	x.den *= y.num;
	x.normalize();
}

template <typename Integer>
Rational<Integer> operator /(const Rational<Integer>& x, const Rational<Integer>& y) {
	return {x.num * y.den, x.den * y.num};
}
