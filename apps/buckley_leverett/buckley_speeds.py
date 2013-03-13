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

def secant( qguess0, qguess1, maxiter=1e5, tol=1e-14 ):
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
        
