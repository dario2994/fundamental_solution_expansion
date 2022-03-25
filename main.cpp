#include <map>
#include <iostream>
#include <vector>
#include <cassert>
#include <functional>
#include "rational.hpp"
#include "big_int.hpp"
using namespace std;


/////////////////////// GLOBAL CONFIGURATIONS ////////////////////////
// Coefficients type
typedef Rational<Num> F; // Rational numbers, also `typedef double F;` works.
// Dimension of the ambient space
const int N = 3; // Coordinates are 0-based, so {0, 1, 2, ..., N-1}.


// Utility functions
bool is_zero(double z) { return abs(z) < 1e-12; }

template <typename Integer>
Integer absolute_value(Integer x) { return (x > 0) ? x : -x; }

Num gcd(const Num& x, const Num& y) { return Num::gcd(x, y); }

ostream& operator <<(ostream& out, const Num& x) {
	vector<char> str;
	x.print(str);
	out << string(str.begin(), str.end());
	
	return out;
}

// Struct representing a polynomial.
struct Polynomial {
	map<vector<int>, F> coef; // exponent -> coefficient
};

// Struct representing a homogeneous polynomial.
struct HomogeneousPolynomial {
	const static int NULL_POLYNOMIAL = -1234567;
	map<vector<int>, F> coef; // exponent -> coefficient
	int deg() {
		if (coef.empty()) return NULL_POLYNOMIAL;
		int d = 0;
		for (int e: coef.begin()->first) d += e;
		return d;
	}
	HomogeneousPolynomial clean() {
		map<vector<int>, F> new_coef;
		for (auto pp: coef) if (!::is_zero(pp.second)) new_coef.insert(pp);
		coef = new_coef;
		return *this;
	}
	bool is_zero() const { return coef.empty(); }
};
typedef HomogeneousPolynomial HP;

ostream& operator<<(ostream& out, HomogeneousPolynomial p) {
	const static vector<string> var_name = {"x", "y", "z", "a", "b", "c", "d", "e"};
	bool is_null = true;
	bool first_monomial = true;
	for (auto pp: p.coef) {
		F c = pp.second;
		F cpos = abs(pp.second);
		if (c == F(0)) continue;
		is_null = false;

		if (c > F(0) and first_monomial) out << "";
		if (c > F(0) and !first_monomial) out << " + ";
		if (c < F(0) and first_monomial) out << "-";
		if (c < F(0) and !first_monomial) out << " - ";

		if (!is_zero(cpos - F(1))) out << cpos;

		bool no_variables = true;
		auto exp = pp.first;
		for (int i = 0; i < (int)exp.size(); i++) {
			if (exp[i] == 0) continue;
			out << var_name[i];
			no_variables = false;
			if (exp[i] > 1) out << "^" << exp[i];
		}
		if (no_variables and is_zero(cpos - F(1))) out << "1";

		first_monomial = false;
	}

	if (is_null) out << 0;

	return out;
}

void operator *=(HP& p, F lambda) {
	for (auto& pp: p.coef) pp.second *= lambda;
}

HP operator *(HP p, F lambda) {
	p *= lambda;
	return p;
}

void operator /=(HP& p, F lambda) {
	for (auto& pp: p.coef) pp.second /= lambda;
}

HP operator /(HP p, F lambda) {
	p /= lambda;
	return p;
}

void operator +=(HP& p, HP q) {
	for (const auto& pp: q.coef) p.coef[pp.first] += pp.second;
}

void operator -=(HP& p, HP q) {
	for (const auto& pp: q.coef) p.coef[pp.first] -= pp.second;
}

HP operator +(HP p, const HP& q) {
	p += q;
	return p;
}

HP operator -(HP p, const HP& q) {
	p -= q;
	return p;
}

HP Derive(const HP& p, int coord) {
	HP q;
	for (const auto& pp: p.coef) {
		vector<int> exponent = pp.first;
		if (exponent[coord] == 0) continue;
		F c = pp.second;
		c *= F(exponent[coord]);
		exponent[coord]--;
		q.coef[exponent] = c;
	}
	return q;
}

HP MultiplyByX(const HP& p, int coord) {
	HP q;
	for (const auto& pp: p.coef) {
		vector<int> exponent = pp.first;
		exponent[coord]++;
		q.coef[exponent] = pp.second;
	}
	return q;
}

// Returns {a, b} where p * x_coord = b + |x|^2 a.
// The polynomials a, b are guaranteed to be harmonic.
// Assumption: p is harmonic.
pair<HP, HP> HarmonicDecompositionMultiplyByX(HP p, int coord) {
	if (p.deg() == HomogeneousPolynomial::NULL_POLYNOMIAL) return {};

	HP p1 = Derive(p, coord).clean();
	if (!p1.is_zero()) p1 /= F(2 * p.deg() + N - 2);
	HP p2 = MultiplyByX(p, coord);
	for (int i = 0; i < N; i++) 
		p2 -= MultiplyByX(MultiplyByX(p1, i), i);
	return {p1, p2};
}


// Struct representing p|x|^h + q|x|^h\log|x|; where p, q are
// k-homogeneous harmonic polynomials.
struct HomogeneousFunction {
	int k, h; // k >= 0
	HP p, q; // q != 0 only if h = 0, 2, 4, ...
	bool has_log() const { return h >= 0 and h%2 == 0; }
	HomogeneousFunction clean() {
		p.clean();
		q.clean();
		return *this;
	}

	bool is_zero() const { return p.is_zero() and q.is_zero(); }
};
typedef HomogeneousFunction H;

ostream& operator<<(ostream& out, H f) {
	f.clean();

	if (!f.p.is_zero()) {
		bool is_p_monomial = f.p.coef.size() == 1;
		if (is_p_monomial) out << f.p;
		else out << "(" << f.p << ")";
		if (f.h != 0) out << "r^" << f.h;
	}

	if (f.has_log() and !f.q.is_zero()) {
		bool is_q_monomial = f.q.coef.size() == 1;
		if (!f.p.is_zero()) out << " + ";

		if (is_q_monomial) out << f.q;
		else out << "(" << f.q << ")";

		if (f.h != 0) out << "r^" << f.h;
		out << "log(r)";
	}
	return out;
}

void operator*=(H& f, F lambda) {
	f.p *= lambda;
	if (f.has_log()) f.q *= lambda;
}

H operator*(H f, F lambda) {
	f *= lambda;
	return f;
}

void operator+=(H& f, const H& g) {
	assert(f.k == g.k and f.h == g.h);
	f.p += g.p;
	if (f.has_log()) f.q += g.q;
}

H operator+(H f, H g) {
	f += g;
	return f;
}

void operator-=(H& f, const H& g) {
	assert(f.k == g.k and f.h == g.h);
	f.p -= g.p;
	if (f.has_log()) f.q -= g.q;
}

H operator-(H f, H g) {
	f -= g;
	return f;
}

H operator-(H f) {
	return f * F(-1);
}

// Returns a homogeneous function of homogeneity (f.k, f.h-2) whose
// Laplacian gives f.
H InverseLaplacian(H f) {
	int k = f.k, h = f.h;
	F a((h + 2) * (2*k + h + N));
	F b (2 * (k + h) + N + 2);
	if (h < -2 or h%2) f.p *= F(1) / a;
	else if (h == -2) {
		swap(f.p, f.q);
		f.q /= F(2 * k + N - 2);
	} else {
		f.q /= a;
		f.p += f.q * (-b);
		f.p /= a;
	}
	f.h += 2;
	return f;
}

// Struct representing a finite sum of HomogeneousFunctions.
struct AsymptoticExpansion {
	map<pair<int,int>, H> components; // (k, h) -> H^{k,h}-component.
	AsymptoticExpansion() {}
	AsymptoticExpansion(H f) { components[{f.k, f.h}] = f; }
	void add(H z) {
		if (components.count({z.k, z.h})) components[{z.k, z.h}] += z;
		else components[{z.k, z.h}] = z;
	}

	void sub(H z) {
		if (components.count({z.k, z.h})) components[{z.k, z.h}] -= z;
		else components[{z.k, z.h}] = -z;
	}
	AsymptoticExpansion clean() {
		map<pair<int,int>, H> new_components;
		for (auto pp: components) {
			auto id = pp.first;
			H h = pp.second;
			h.clean();
			if (!h.is_zero()) new_components[id] = h;
		}
		components = new_components;
		return *this;
	}

	int min_homogeneity() {
		int min_deg = 1e9;
		clean();
		for (auto pp: components)
			min_deg = min(min_deg, pp.first.first + pp.first.second);
		return min_deg;
	}
	
	void print(int max_homogeneity = 1<<30) {
		vector<H> to_print;
		for (auto pp: components) {
			if (pp.first.first + pp.first.second <= max_homogeneity) 
				to_print.push_back(pp.second);
		}
		sort(to_print.begin(), to_print.end(), [&](H f, H g) {
			return f.k + f.h < g.k + g.h;
		});
		bool is_first = true;
		for (H f: to_print) {
			if (!is_first) cout << " + ";
			cout << f;
			is_first = false;
		}
		cout << endl;
	}
};
typedef AsymptoticExpansion AE;

ostream& operator<<(ostream& out, const AE& f) {
	out << "{";
	bool first_element = true;
	for (auto pp: f.components) {
		if (!first_element) out << ", ";
		first_element = false;
		out << pp.second;
	}
	out << "}";
	return out;
}

void operator*=(AE& f, F lambda) {
	for (auto& comp: f.components) comp.second *= lambda;
}

AE operator*(AE f, F lambda) {
	f *= lambda;
	return f;
}

void operator+=(AE& f, AE g) {
	for (auto comp: g.components) f.add(comp.second);
}

AE operator+(AE f, AE g) {
	f += g;
	return f;
}

void operator-=(AE& f, AE g) {
	for (auto comp: g.components) f.sub(comp.second);
}

AE operator-(AE f, AE g) {
	f -= g;
	return f;
}

AE operator-(const AE& f) { return f * F(-1); }

AE InverseLaplacian(const AE& f) {
	AE g;
	for (const auto& comp: f.components) {
		auto comp2 = InverseLaplacian(comp.second);
		g.add(comp2);
	}
	return g;
}

AE MultiplyByX(H f, int coord) {
	H g1;
	H g2;

	g1.k = f.k-1;
	g1.h = f.h+2;

	g2.k = f.k+1;
	g2.h = f.h;

	pair<HP, HP> p = HarmonicDecompositionMultiplyByX(f.p, coord);
	g1.p = p.first;
	g2.p = p.second;

	if (f.has_log()) {
		pair<HP, HP> q = HarmonicDecompositionMultiplyByX(f.q, coord);
		g1.q = q.first;
		g2.q = q.second;
	}

	AE res;
	if (f.k > 0) res.add(g1);
	res.add(g2);

	return res;
}

AE Derive(H f, int coord) {
	H g1;
	H g2;

	g1.k = f.k-1;
	g1.h = f.h;

	g1.p = Derive(f.p, coord);
	if (f.has_log()) g1.q = Derive(f.q, coord);

	g2.k = f.k;
	g2.h = f.h-2;

	g2.p = f.p * F(f.h);
	if (f.has_log()) {
		g2.p += f.q;
		g2.q = f.q * F(f.h);
	}

	AE res = MultiplyByX(g2, coord);
	if (f.k > 0) res.add(g1);

	return res;
}


AE MultiplyByX(const AE& f, int coord) {
	AE g;
	for (const auto& comp: f.components) g += MultiplyByX(comp.second, coord);
	return g;
}

AE Derive(const AE& f, int coord) {
	AE g;
	for (const auto& comp: f.components) g += Derive(comp.second, coord);
	return g.clean();
}

AE Derive(const AE& f, vector<int> coords) {
	AE g = f;
	for (int coord: coords) g = Derive(g, coord);
	return g;
}

AE MultiplyByPolynomial(const AE& f, Polynomial z) {
	AE res;
	for (auto pp: z.coef) {
		AE curr = f;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < pp.first[i]; j++) curr = MultiplyByX(curr, i);
		}
		res += curr * pp.second;
	}
	return res;
}

AE operator*(const AE& f, Polynomial z) {
	return MultiplyByPolynomial(f, z);
}


// Given an elliptic operator L whose principal symbol at 0 coincides
// with the Laplacian, this function returns the first terms
// of the series expansion of the fundamental solution of the operator
// L.
AE FundamentalSolution(function<AE(AE)> L, int precision = 3) {
	auto Laplacian = [&](AE f) {
		AE res;
		for (int i = 0; i < N; i++) res += Derive(Derive(f, i), i);
		return res;
	};

	// T = (Laplacian - L)(Laplacian^{-1})
	function<AE(AE)> T = [&](AE f) {
		AE g = InverseLaplacian(f);
		AE Tf = Laplacian(g) - L(g);
		return Tf.clean();
	};

	// std_fund_sol is the fundamental solution for the Laplacian
	// operator, up to multiplicative normalization.
	H _std_fund_sol = {.k = 0, .h = 2-N};
	if (N == 2) _std_fund_sol.q.coef = {{vector<int>(N, 0), F(1)}};
	else _std_fund_sol.p.coef = {{vector<int>(N, 0), F(1)}};
	AE std_fund_sol(_std_fund_sol);
	
	// Tdelta = T(delta) = Laplacian(std_fund_sol) - L(std_fund_sol).
	// We cannot call T(delta) directly because the \delta distribution
	// cannot be represented as an AsymptoticExpansion.
	AE Tdelta = -L(std_fund_sol);
	Tdelta.clean();

	// sum_Tdelta = T(delta) + T^2(delta) + T^3(delta) + ...
	AE sum_Tdelta = Tdelta;
	for (int i = 2; i <= precision; i++) {
		Tdelta = T(Tdelta);
		sum_Tdelta += Tdelta;
	}
	sum_Tdelta.clean();

	AE fund_sol = std_fund_sol + InverseLaplacian(sum_Tdelta);
	fund_sol.clean();
	
	// Here we check that L(fund_sol) is, up to high order terms, null.
	// This works because the delta does not appear since it cannot
	// be represented as Asymptotic Expansion.
	AE almost_zero = L(fund_sol).clean();
	assert(almost_zero.min_homogeneity() >= precision + 1 - N);
	
	return fund_sol;
}


int main() {
	/* These lines compute the fundamental solution, in dimension N = 2,
	 * of the operator L which is obtained by the Laplacian through
	 * the transformation Phi(x_0, x_1) = (x_0, x_1 - x_0^2).
	 * The result can be used to deduce that in even dimension the 
	 * construction is not covariant.
	Polynomial pol1;
	pol1.coef = {{{0, 0}, F(1)}};
	Polynomial pol2;
	pol2.coef = {{{1, 0}, F(4)}};
	Polynomial pol3;
	pol3.coef = {{{0, 0}, F(1)}, {{2, 0}, F(4)}};
	Polynomial pol4;
	pol4.coef = {{{0, 0}, F(2)}};
	function<AE(AE)> L = [&](AE f) {
		return Derive(f, {0, 0}) + Derive(f, {1, 1}) * pol3 + Derive(f, {1, 0}) * pol2
				 + Derive(f, 1) * pol4;
	};
	
	cout << FundamentalSolution(L, 5) << endl;
	*/
	
	/* These lines compute exactly what was described above but in the
	 * case N = 3.
	 * The result shows, in this simple setting, the covariance of
	 * the construction.
	 */
	Polynomial pol1;
	pol1.coef = {{{0, 0, 0}, F(1)}};
	Polynomial pol2;
	pol2.coef = {{{1, 0, 0}, F(4)}};
	Polynomial pol3;
	pol3.coef = {{{0, 0, 0}, F(1)}, {{2, 0, 0}, F(4)}};
	Polynomial pol4;
	pol4.coef = {{{0, 0, 0}, F(2)}};
	function<AE(AE)> L = [&](AE f) {
		return Derive(f, {0, 0}) * pol1 
		     + Derive(f, {1, 1}) * pol3
		     + Derive(f, {2, 2}) * pol1
		     + Derive(f, {0, 1}) * pol2
			 + Derive(f, 1) * pol4;
	};
	// Prints the terms of homogeneity <= 0 of the fundamental
	// solution of the operator L.
	FundamentalSolution(L).print(0);
}
