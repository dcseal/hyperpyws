Model equations
===============

This directory is where the model equations are defined for use in the
applications.

A hyperbolic conservation law, 

    q_t + ( f(q) )_x = 0,

is defined by the flux function, $f$ together with a method for diagonalizing
it.  See: "Finite Volume Methods for hyperbolic problems", LeVeque, 2002.

In each problem, one is asked to define the following functions:

* f   - the flux function
* J   - the jacobian of the flux function (f'(q))
* eig - a routine for computing the eigenvalues of J
* R   - a routine for constructing a matrix of the right eigenvectors of J
* L   - a routine for constructing the left eigenvectors of J as rows (=inv(J))
* MaxWaveSpeed - a routine for determining the maximum wave speed given a state q.

Examples
--------

[HYPERPYWS](../../README.md) currently has the following model equations
implemented:

### Scalar equations ###

* Advection
* Burger's equation
* Buckley-Leverett

### Systems of equations ###

* Shallow water equations
* Euler's equations

In order to implement a new equation, two steps need to be followed:

1. Define the flux function in this directory.
2. Copy an existing application and run it.

For example, the following line in

    ./apps/euler/shock-tube/shock_tube.py

calls

    from hyperpyws.model_equations.euler import Euler1D

that exists in this part of the code.

