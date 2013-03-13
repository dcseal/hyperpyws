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
    """ derivative of flux function. """

    return (2.*M*q*(1.-q)) / (q**2 + M*(1.-q)**2)**2

def G(q, fqr=0, qr=0):
    """ Function used for finding the zero for R-H conditions, which are the
    solution to the equation:

        f'(q) = (f(q) - f(qr)) / (q - qr).

    Note that the known values are qr and fqr = f(qr), and we are looking for
    a solution q to that system.

    This is used in the single case of a Riemann problem
    that produces a Rarefaction + Shock.

    When G(q) = 0, we get a solution to the above equations.
    """
    return ( fp( q ) - ( f(q) - fqr ) / (q-qr) )

def Gp(q, fqr=0, qr=0 ):
    """ Derivative of G from above.  Used for Newton iteration."""

    raise NotImplemented('you suck and need to implement a derivative!')

def secant_right( qguess0=0.4999, qguess1=0.49999, maxiter=1e5, tol=1e-14 ):
    """ Secant method for Newton iteration. 
    
    q = q - G(q) * (q-qnm1) / ( G(q) - G(qnm1) ) 
    """

    # initial guesses:
    qnm1 = qguess0;  Gnm1 = G(qnm1)
    qn   = qguess1

    n_iters = 0
    while( n_iters < maxiter ):

        Gn   = G(qn)
        if( abs(Gn) < tol ):
            return qn

        # print('Gn = %f, Gnm1 = %f' % (Gn, Gnm1) )

        # update for qn (and swap out the old q):
        tmp = qn - Gn * (qn-qnm1) / ( Gn - Gnm1 )

        # update q (and save old variable):
        qnm1 = qn
        Gnm1 = Gn
        qn   = tmp

def secant_left( qguess0=0.13, qguess1=0.130000001, maxiter=1e5, tol=1e-14 ):
    """ Secant method for Newton iteration. 
    
    q = q - G(q) * (q-qnm1) / ( G(q) - G(qnm1) ) 

    For our problem, qs = 0.13397459621556132.
    """

    # initial guesses:
    qnm1 = qguess0;  Gnm1 = G(qnm1)
    qn   = qguess1

    qr  = 1
    fqr = f(qr)

    n_iters = 0
    while( n_iters < maxiter ):

#       print('qn = %f' % qn )
        Gn   = G(qn, fqr, qr)
        if( abs(Gn) < tol ):
            return qn

        # update for qn (and swap out the old q):
        tmp = qn - Gn * (qn-qnm1) / ( Gn - Gnm1 )

        # update q (and save old variable):
        qnm1 = qn
        Gnm1 = Gn
        qn   = tmp

def get_shock_speeds():

    qs_left  = secant_left()
    qs_right = secant_right()

    print('left  value qs = %2.15e, speed = %2.15e' % (qs_left,  fp(qs_left)  ) )
    print('right value qs = %2.15e, speed = %2.15e' % (qs_right, fp(qs_right) ) )

    return( qs_left, qs_right )
