#===============================================================================
#
# Module written to solve a single Riemann problem.  We are closely following
# the text of Toro, "Riemann Solvers and Numerical Methods for Fluid Dynamics"
#
# The Riemann problem consists of three waves, which we can sort in increasing
# order by { s^1, s^2, s^3 }.  
#
# For Euler, the center wave along characteristic s^2 is always a
# contact discontinuity, which means that the only quantity that can change
# over it is the density.
#
# We'll number the regions as following:
#
# Region 1: W_L  ( Left data defining the Riemann problem )
#
# Region 2: Ws_L ( constant data between s^1 and s^2 )
#
# Region 3: Ws_R ( constant data between s^2 and s^3 )
#
# Region 4: W_R  ( Right data defining the Riemann problem )
#
# The central speed, s^2 is always given by solving for u^*.
#
# In the case of a rarefaction, there is missing data between the two
# sound speeds.  
#
# In the case of a left rarefaction, the speeds enclosing the
# fan are governed by
#
#      S_HL = u_l - c_l;  S_TL = us - cs_l,
#
# In the case of a left Rarefaction, data for region 2 is defined by, 
#
#   rho_s_L = rho_L ( ps / p_L )**(1./gamma).
#
# The full fan for a Left rarefaction is filled out via (eqn 4.56, Toro):
#
#      rho = rho_l   * ( 2/gp1 + (gm1/(c_l*gp1))*( u_l - x/t ) )**(  2/gm1);
#        u = 2/gm1   * ( c_l + gm1/2 *u_l + x/t )
#    press = press_l * ( 2/gp1 + (gm1/(c_l*gp1))*( u_l - x/t ) )**(2*gamma/gm1).
#
# The values for a right rarefaction is filled out via (eqn 4.63, Toro):
#
#      rho = rho_r   * ( 2/gp1 - (gm1/(c_r*gp1))*( u_r - x/t ) )**(  2/gm1);
#        u = 2/gm1   * ( -c_r + gm1/2 *u_r + x/t )
#    press = press_r * ( 2/gp1 - (gm1/(c_r*gp1))*( u_r - x/t ) )**(2*gamma/gm1).
#
#===============================================================================

from math import sqrt

#===============================================================================
# TODO - this could be inserted into a module if so desired
#===============================================================================
def NewtonSolve( f, fp, yguess, tol=1e-13, maxiter=1000 ):
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
    """ These function calls look identical for left/right values.

    The only thing that changes is the W_L vs. W_R.
    
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

    b = 0 #  ( <<< ideal gas, b = 0.  co-volume gas, b \neq 0. )

    # Data dependent constants:
    c = sqrt( abs( gamma*p / rho ) )
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
    """ p_star = ps is found by solving f(ps) = 0. """
    return f_K( ps, W_r ) + f_K( ps, W_l ) + ( u_r - u_l )

#===============================================================================
def fp( ps ):
    """ Derivative of f from above. """
    return fp_K( ps, W_r ) + fp_K( ps, W_l )

#===============================================================================
def LeftFan( x, t, W_l ):
    """ Data in the case of a left rarefaction fan.

    User is responsible for knowing correct bounds of x to supply.
    """

    rho_l = W_l[0]
    u_l   = W_l[1] 
    p_l   = W_l[2]

    rho   = rho_l   * ( 2/gp1 + (gm1/(c_l*gp1))*( u_l - x/t ) )**(  2/gm1);
    u     = 2/gm1   * ( c_l + gm1/2 *u_l + x/t )
    press = press_l * ( 2/gp1 + (gm1/(c_l*gp1))*( u_l - x/t ) )**(2*gamma/gm1)

    return [rho, u, press]

def RightFan( x, t, W_r ):
    """ Data for values inside the RightFan.

    User is responsible for knowing correct bounds of x to supply.
    """

    rho_r = W_r[0]
    u_r   = W_r[1] 
    p_r   = W_r[2]

    rho   = rho_r   * ( 2/gp1 - (gm1/(c_r*gp1))*( u_r - x/t ) )**(  2/gm1);
    u     = 2/gm1   * ( -c_r + gm1/2 *u_r + x/t )
    press = press_r * ( 2/gp1 - (gm1/(c_r*gp1))*( u_r - x/t ) )**(2*gamma/gm1)
 
    return [rho, u, press]

#===============================================================================
# TODO: it would be nice to have a better naming convention, and cleaner
# output data.
#

print('Initial data:')
print('    W_l = (%2.3f, %2.3f, %2.3f) ' % (W_l[0],W_l[1],W_l[2]) )
print('    W_r = (%2.3f, %2.3f, %2.3f) ' % (W_r[0],W_r[1],W_r[2]) )
print(' ')

# stupid initial guess, take the arithmetic mean for p:
pstar,niters = NewtonSolve( f, fp, 0.5*(p_l+p_r), tol=1e-14 )

# compute ustar, which is only a function of pstar:
ustar        = 0.5*( u_l + u_r ) + 0.5*( f_K( pstar, W_r ) - f_K( pstar, W_l ) )

Ws_l = [float('NaN'), ustar, pstar]
Ws_r = [float('NaN'), ustar, pstar]

print('Newton solve converged in %d iterations' % niters )
print(' ')
print('Pressure and velocity inside the entire the star region:')
print('    p_star = %2.15e' % (pstar) )
print('    u_star = %2.15e' % ustar   )
print(' ')

#===============================================================================
# Leftmost characteristic (s^1)
#===============================================================================
if( pstar > p_l ):
    print("Left Shock:")
    shock_speed_left = u_l - c_l*sqrt( (gp1/(2.*gamma))*(pstar/p_l)+gm1/(2*gamma))
    print("     left shock speed = %f" % shock_speed_left )
else:
    print("Left Rarefaction: ")

    rho_star_l = rho_l * ( pstar / p_l )**(1./gamma)
    Ws_l[0] = rho_star_l

    cs = sqrt( abs( gamma*pstar / rho_star_l ) )

    c_rarefaction_left_l = u_l   - c_l
    c_rarefaction_left_r = ustar - cs
    print("    Left  foot speed  = %2.15e" % c_rarefaction_left_l )
    print("    Right foot speed  = %2.15e" % c_rarefaction_left_r )

    # This is the correct value when the left side is a rarefaction
    print("    rho_star_left     = %2.15e " % rho_star_l )

#===============================================================================
# Contact discontinuity  (s^2)
#===============================================================================
print(' ')
print('Contact Discontinuity:')
print('    speed (ustar)     = %2.15e ' % ustar )
print(' ')

#===============================================================================
# Rightmost characteristic (s^3)
#===============================================================================
if( pstar > p_r ):
    print("Right Shock:")
    shock_speed_right = u_r + c_r*sqrt( (gp1/(2.*gamma))*(pstar/p_r)+gm1/(2*gamma))
    print("    Right shock speed = %2.15e" % shock_speed_right )

    # (c.f. equation 4.57, when s^3 produces a shock)
    rho_star_r = \
        rho_r * ( (pstar/p_r) + (gm1/gp1) ) / ( (gm1/gp1)* pstar/ p_r + 1.0 ) 
    print("    rho_star_right    = %2.15e " % rho_star_r )

    Ws_r[0] = rho_star_r

else:
    print("Right Rarefaction " )

print(' ')
print("In summary, we have the following values for each region:")
print('    Region 1: W = ', W_l )
print('    Region 2: W = ', Ws_l )
print('    Region 3: W = ', Ws_r )
print('    Region 4: W = ', W_r )

print("If you saw NaN, that means we haven't included that option yet!")
#===============================================================================
