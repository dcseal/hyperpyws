#===============================================================================
#
# Module written to solve a single Riemann problem.  We are closely following
# the text of Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics"
#
#===============================================================================

from math import sqrt

# Left and right state variable.
# These are the non-conserved quantities, defined by
# W_K = [rho, u, press], K=l,r.

# ( TODO: <<<< make these user supplied parameters to a function <<<)
W_l = [10.0, 0, 100.0]
W_r = [1.0 , 0,   1.0]  
#W_l = [3.0, 0, 3.0]
#W_r = [1.0, 0, 1.0]  

rho_l = W_l[0]; u_l = W_l[1]; p_l = W_l[2]
rho_r = W_r[0]; u_r = W_r[1]; p_r = W_r[2]

# easy access shortcut variables:
gamma = 1.4;
gm1   = gamma-1.;
gp1   = gamma+1.;

# sound speeds for left and right sides:
c_l = sqrt( abs( gamma*p_l / rho_l ) )
c_r = sqrt( abs( gamma*p_r / rho_r ) )

#===============================================================================
def f_K( ps, W ):
    """ These are identical for left/right values.
    
    See proposition 4.1 in Toro's book.
    """

    # the left (or right) constant variables:
    rho = W[0]
    u   = W[1]
    p   = W[2]

    b = 0 # <<< ideal gas, b = 0.  co-volume gas, b \neq 0.

    # Data dependent constants:
    c = sqrt( abs( gamma*p / rho ) )      # sound speed
    A = 2.0 * (1.0 - b*rho) / ( gp1 * rho )
    B = ( gm1 / gp1 ) * p

    if( ps > p ):
        # shock:
        return (ps-p) * sqrt( A / ( ps + B ) )
    else:
        # rarefaction:
        return ((2.0*c*(1.-b*rho))/gm1) * ((ps/p)**(gm1/(2.0*gamma))-1.)

#===============================================================================
def fp_K( ps, W ):
    """ Derivative of the function, f_K."""

     # the left (or right) constant variables:
    rho = W[0]
    u   = W[1]
    p   = W[2]

    b = 0 # <<< ideal gas, b = 0.  co-volume gas, b \neq 0.

    # Data dependent constants:
    c = sqrt( abs( gamma*p / rho ) )      # sound speed
    A = 2.0 * (1.0 - b*rho) / ( gp1 * rho )
    B = ( gm1 / gp1 ) * p

    if( ps > p ):
        # shock:
        return sqrt( A / ( ps + B ) )*( 1.0 - ( ps-p )/ (2.0*(B+ps) )  )
    else:
        # rarefaction:
        # return c/(gamma*ps) * (ps/p)**(0.5*gm1/gamma)

        # this line is what is in Toro's book (c.f. eqn 4.37):
        return 1.0 / (rho * c) * (ps/p)**(-0.5*gp1/gamma)

#===============================================================================
def f( ps ):
    return f_K( ps, W_r ) + f_K( ps, W_l ) + ( u_r - u_l )

#===============================================================================
def fp( ps ):
    return fp_K( ps, W_r ) + fp_K( ps, W_l )

#===============================================================================
def NewtonSolve( f, fp, yguess, tol=1e-13, maxiter=10000 ):
    """ Newton iteration to solve for f(y) = 0.
    The user is requested to supply an initial guess, yguess, as well as the
    derivative of the flux function fp.

    This method currently only works on scalar problems.
    """

    yn   = yguess
    fnm1 = f(yn)

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
# TODO - THIS SECTION NEEDS TO BE MASSIVELY CLEANED UP!
#
# Better naming convention, ....  better output data ...

print('Initial data:')
print('    W_l = (%f, %f, %f) ' % (W_l[0],W_l[1],W_l[2]) )
print('    W_r = (%f, %f, %f) ' % (W_r[0],W_r[1],W_r[2]) )
print(' ')

# stupid initial guess, take the arithmetic mean for p:
pstar,niters = NewtonSolve( f, fp, 0.5*(p_l+p_r), tol=1e-14 )

# compute ustar, which is only a function of pstar:
ustar        = 0.5*( u_l + u_r ) + 0.5*( f_K( pstar, W_r ) - f_K( pstar, W_l ) )

print('Newton solve converged in %d iterations' % niters )
print('    p_star = %2.15e' % (pstar) )
print(' ')

# information about left part:
if( pstar > p_l ):
    print("Left Shock:")
    shock_speed_left = u_l - c_l*sqrt( (gp1/(2.*gamma))*(pstar/p_l)+gm1/(2*gamma))
    print("     left shock speed = %f" % shock_speed_left )
else:
    print("Left Rarefaction: ")

    rho_star_l = rho_l * ( pstar / p_l )**(1./gamma)
    cs = sqrt( abs( gamma*pstar / rho_star_l ) )

    c_rarefaction_left_l = u_l   - c_l
    c_rarefaction_left_r = ustar - cs
    print("    Left  foot speed  = %2.15e" % c_rarefaction_left_l )
    print("    Right foot speed  = %2.15e" % c_rarefaction_left_r )

# ustar provides location of contact discontinuity:
print('Contact Discontinuity:')
print('    speed (ustar)     = %2.15e ' % ustar )

# more information about solution on right edge:
if( pstar > p_r ):
    print("Right Shock:")
    shock_speed_right = u_r + c_r*sqrt( (gp1/(2.*gamma))*(pstar/p_r)+gm1/(2*gamma))
    print("    Right shock speed = %2.15e" % shock_speed_right )
else:
    print("Right Rarefaction " )

print(' ')
rho_star_r = \
    rho_r * ( (pstar/p_r) + (gm1/gp1) ) / ( (gm1/gp1)* pstar/ p_r + 1.0 ) 
print("rho_star_right = %f " % rho_star_r )
print("rho_star_left  = ??? " )

#===============================================================================