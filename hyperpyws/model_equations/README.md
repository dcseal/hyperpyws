#===============================================================================
= Model equations =
#===============================================================================

This directory is where the model equations are defined for use in the
applications.

A hyperbolic conservation law, $q_t + ( f(q) )_x = 0$,
is defined by the flux function, $f$ together with a method for diagonalizing
it.  See: "Finite Volume Methods for hyperbolic problems", LeVeque, 2002.

In each problem, one is asked to define the following functions:

* f   - the flux function
* J   - the jacobian of the flux function (f'(q))
* eig - a routine for computing the eigenvalues of J
* R   - a routine for constructing a matrix of the right eigenvectors of J
* L   - a routine for constructing the left eigenvectors of J as rows (=inv(J))
* MaxWaveSpeed - a routine for determining the maximum wave speed given a state q.
