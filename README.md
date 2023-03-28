# Iterative construction of the fundamental solution of a 2nd order elliptic operator.
Companion software for the paper [Expansion of the fundamental solution of a second-order elliptic operator with analytic coefficients](https://arxiv.org/abs/2110.15104).

We implement in `C++` the construction given in the paper.

# Compilation
```g++ -O2 main.cpp -o main```

# Usage
You shall configure the dimension (`N`) at the top of the function_representations.hpp file.

Then, in the `main()` you shall describe the operator you want to consider.
The example currently implemented is for the operator (in dimension `N=3`, with the three variables `x_0`, `x_1`, `x_2`)

![](https://latex.codecogs.com/svg.image?L&space;=&space;\partial_{00}&space;&plus;&space;(1&plus;4x_0^2)\partial_{11}&space;&plus;&space;\partial_{22}&space;&plus;&space;4x_0\partial_{01}&space;&plus;&space;2\partial_1)

<!--- (L = \partial_{00} + (1+4x_0^2)\partial_{11} + \partial_{22} + 4x_0\partial_{01} + 2\partial_1) -->

Then, executing `./main` will print the first few homogenous terms of the expansion of the fundamental solution for `L`.

# Documentation

An operator is a `function<AE(AE)>` taking an `AsymptoticExpansion=AE` and giving back an `AsymptoticExpansion=AE`.
It must be elliptic and its principal symbol must be the same of the Laplacian. See the example in `main.cpp`.
The field of coefficient is denoted by `F`.

The coefficients defining the operator shall be instances of the class `Polynomial`. 
A polynomial is a vector of pairs `{exponent, coefficient}`. 
The exponent is a vector of length `N` (where `N` is the dimension), and the i-th entry represent the exponent of `x_i`.
For example, the polynomial ![](https://latex.codecogs.com/svg.image?1&space;&plus;&space;4x_0^2) is constructed by `{{{0, 0, 0}, F(1)}, {{2, 0, 0}, F(4)}}`.
