#include <map>
#include <iostream>
#include <vector>
#include <cassert>
#include <functional>
#include "rational.hpp"
#include "big_int.hpp"
// Set here the field of coefficients and the dimension
#include "function_representations.hpp"
using namespace std;




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
	Polynomial pol1; // = 1
	pol1.coef = {{{0, 0, 0}, F(1)}};
	Polynomial pol2; // = 4x
	pol2.coef = {{{1, 0, 0}, F(4)}};
	Polynomial pol3; // = 1 + 4x²
	pol3.coef = {{{0, 0, 0}, F(1)}, {{2, 0, 0}, F(4)}};
	Polynomial pol4; // = 2
	pol4.coef = {{{0, 0, 0}, F(2)}};
	function<AE(AE)> L = [&](AE f) { // d_xx + (1+4x²)d_yy + d_zz + 4x d_xy + 2 d_y
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
