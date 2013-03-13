# single free parameter for buckley-leverett (TODO - should pull from problem
# setup ... )
M = 1.0/3.0

#-----------------------------------------------------------------------------
def f (q):
    """ Flux function f(q). """
    
    u0   = q**2
    u1   = (1.0-q)**2
    
   
    return u0 / (u0 + M*u1)
  
#-----------------------------------------------------------------------------
def fp(q):
    """ derivative of flux function for Buckley-Leverett. """

    return (2.*M*q*(1.-q)) / (q**2 + M*(1.-q)**2)**2

#-----------------------------------------------------------------------------
def secant_method( G, yguess0, yguess1, maxiter=1e5, tol=1e-14):
    """ Secant method for finding the zeros of the function G(y).

    The method depends on the two previous values, so the user needs to supply
    two initial guesses.

    Method:
    
        y = y - G(y) * (y-ynm1) / ( G(y) - G(ynm1) ) 

    Parameters:
    ===========

        G : callable function.

        yguess0, yguess1 : initial guesses for iteration

        maxiter : maximumun number of iterations allowed

        tol : tolerance used to determine if we found the zero.

    Returns:

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


#-----------------------------------------------------------------------------
def main():
    """Single function to compute the two shock speeds that are present for
    our numerical example.
   
    This has not been written in a general fashion - certainly we could
    rewrite this so it could be used again in the future.
    """

    # Riemman state that remains static:
    q_fixed  = 1
    f_qfixed = f( q_fixed )
    def G(qs):
        """ zero function that needs to be solved for left-state of our
        example Riemman problem

        R-H conditions:

            f'(qs) = (f(qs) - f(q_fixed)) / (q - q_fixed).

        """
        return ( fp( qs ) - ( f(qs) - f_qfixed ) / (qs-q_fixed) )


    qs_left_problem = secant_method( G, 0.13, 0.130000001 )

    print('left  value qs = %2.15e, speed = %2.15e' % 
        (qs_left_problem, fp(qs_left_problem)  ) )

    q_fixed  = 0
    f_qfixed = f(q_fixed)

    qs_right_problem = secant_method( G, 0.48, 0.480000001 )
    print('right value qs = %2.15e, speed = %2.15e' % 
        (qs_right_problem, fp( qs_right_problem) ) )

    return( qs_left_problem, qs_right_problem )
