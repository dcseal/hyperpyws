#==============================================================================#
# This file is part of HYPERPYWS: Hyperbolic Python WENO Solver
#
#
#   This software is made available for research and instructional use only.
#   You may copy and use this software without charge for these non-commercial
#   purposes, provided that the copyright notice and associated text is
#   reproduced on all copies.  For all other uses (including distribution of
#   modified versions), please contact the author at the address given below.
# 
#   *** This software is made available "as is" without any assurance that it
#   *** will work for your purposes.  The software may in fact have defects, so
#   *** use the software at your own risk.
#
# License: GPL, see COPYING for details
#
# Copyright (C) 2013 
#
#    David Seal,  seal@math.msu.edu,  Michigan State University
#    Yaman Guclu, guclu@math.msu.edu, Michigan State University
#
#===============================================================================

# single module describing iterative solvers

#===============================================================================
# Newton iteration - for a scalar problem.
#
# In the case of a linear problem, the Jacobian, fp(y) should be a matrix that
# gets inverted.
#
#===============================================================================
def NewtonSolveScalar( f, fp, yguess, tol=1e-13, maxiter=1000 ):
    """ Newton iteration to solve for f(y) = 0.
    The user is requested to supply an initial guess, yguess, as well as the
    derivative of the flux function fp.

    This method currently only works on scalar problems.

    Parameters:
    ===========

        f  : callable function.

        fp : analytical derivative of f.

        yguess : initial guesses for iteration

        tol : tolerance used to determine if we found the zero.

        maxiter : maximumun number of iterations allowed

    Returns:
    ========

        y : solution to f(y) = 0.

        n : number of iterations needed to get convergence

    """

    # initial setup
    yn   = yguess
    fnm1 = f(yn)

    # main loop.  keep track of number of iterations
    n = 0
    while( n < maxiter ):

        fn = f(yn)
        yn = yn - fn / fp( yn )
        n += 1
        if( abs(fn) < tol ):
            return yn, n
        
    print("failed to converge;  n = %d" % n)
    raise Exception  # (<<<<< TODO - which exception should we raise? <<<< )

#===============================================================================



#===============================================================================
# Secant method (for scalar problems)
#===============================================================================
def SecantSolveScalar( G, yguess0, yguess1, maxiter=1e5, tol=1e-14):
    """ Secant method for finding the zeros of the function G(y).

    The method depends on the two previous values, so the user needs to supply
    two initial guesses.

    Method:
    
        yn = yn - G(yn) * (yn-ynm1) / ( G(yn) - G(ynm1) ) 

    Parameters:
    ===========

        G : callable function.

        yguess0, yguess1 : initial guesses for iteration

        maxiter : maximumun number of iterations allowed

        tol : tolerance used to determine if we found the zero.

    Returns:
    ========

        y : solution to G(y) = 0.

    """

    # initial setup:
    ynm1 = yguess0;  
    Gnm1 = G(ynm1)
    yn   = yguess1

    n_iters = 0
    while( n_iters < maxiter ):

        Gn   = G(yn)
        if( abs(Gn) < tol ):
            return yn

        # print('Gn = %f, Gnm1 = %f' % (Gn, Gnm1) )

        # update for qn (and swap out the old q):
        tmp = yn - Gn * (yn-ynm1) / ( Gn - Gnm1 )

        # update y (and save old variable):
        ynm1 = yn
        Gnm1 = Gn
        yn   = tmp

#===============================================================================

