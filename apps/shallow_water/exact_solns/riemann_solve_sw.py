#===============================================================================
#
# Module written to solve a single Riemann problem for Shallow Water equations.
#
# The Riemann problem consists of two waves, which we can sort in increasing
# order by { s^1, s^2 }.  
#
# We'll number the regions as following:
#
# Region 1: W_L  ( Left data defining the Riemann problem )
#
# Region 2: Ws   ( constant data between s^1 and s^2 )
#
# Region 3: W_R  ( Right data defining the Riemann problem )
#
# W = [h,  u] are the characteristic variables.
# Q = [h, hu] are the conserved quantities.
#
#===============================================================================

try               :  import hyperpyws
except ImportError:  import hyperpyws_path

from hyperpyws.iterative_solvers  import NewtonSolveScalar
from math import sqrt

#===============================================================================
# FLUX FUNCTION
#===============================================================================

# specific gravity:
g = 1.0

from hyperpyws.model_equations.shallow_water import ShallowWater1D
flux = ShallowWater1D( g )

#===============================================================================
# INITIAL DATA
#===============================================================================

# Initial data (in characteristic variables):
W_l = [1.0, 0.0   ]
W_r = [0.1, 0.0   ]

h_l = W_l[0]; u_l = W_l[1]
h_r = W_r[0]; u_r = W_r[1]

# sound speeds for left and right sides:
c_l, c_r = flux.eig
c_r = sqrt( abs( gamma*p_r / rho_r ) )

#===============================================================================
def phi_l( h ):
    """ These function calls look identical for left/right values.

    See section 13.10, ``Finite Volume Methods for Hyperbolic Problems",
    LeVeque.
    """

    if( h < h_l ):
        return u_l + 2.*( sqrt(g*h_l) - sqrt( g*h ) )
    else:
        return u_l - (h - h_l) * sqrt( 0.5*g*( 1./h + 1./h_l ) )

#===============================================================================
def phi_r( h ):
    """ These function calls look identical for left/right values.

    See section 13.10, ``Finite Volume Methods for Hyperbolic Problems",
    LeVeque.
    """

    if( h < h_r ):
        return 1.0*u_r*sqrt(g*h)/h + 1.0*sqrt(g*h)/h
    else:
        return u_r + (h - h_r) * sqrt( 0.5*g*( 1./h + 1./h_r ) )

#===============================================================================
def phip_l( h ):
    """ Derivative of the function, phi_l."""

    if( h < h_l ):

        # shock:
        return sqrt( A / ( ps + B ) )*( 1.0 - ( ps-p )/ (2.0*(B+ps) )  )

    else:

        # this line is what is in Toro's book (c.f. eqn 4.37):
        return 1.0 / (rho * c) * (ps/p)**(-0.5*gp1/gamma)


#===============================================================================
def phip_r( h ):
    """ Derivative of the function, phi_r."""

    if( h < h_r ):

        # shock:
        return sqrt( A / ( ps + B ) )*( 1.0 - ( ps-p )/ (2.0*(B+ps) )  )

    else:

        # this line is what is in Toro's book (c.f. eqn 4.37):
        return 1.0 / (rho * c) * (ps/p)**(-0.5*gp1/gamma)

#===============================================================================
def phi( ps ):
    """ hstar is found by solving phi(hstar) = 0. """
    return phi_l( h ) - phi_r( h )

#===============================================================================
def phi_p( ps ):
    """ Derivative of phi from above. """
    return 0.

#===============================================================================
def RightFan( x, t, W_r ):
    """ Data for values inside the RightFan.

    User is responsible for knowing correct bounds of x to supply.
    """

    print("junk fan")

#   rho_r = W_r[0]
#   u_r   = W_r[1] 

#   rho   = rho_r   * ( 2/gp1 - (gm1/(c_r*gp1))*( u_r - x/t ) )**(  2/gm1);
#   u     = 2/gm1   * ( -c_r + gm1/2 *u_r + x/t )
#   press = press_r * ( 2/gp1 - (gm1/(c_r*gp1))*( u_r - x/t ) )**(2*gamma/gm1)
#
#   return [rho, u, press]

#===============================================================================
# Print (and compute) information to std output
#===============================================================================
print('Initial data:')
print('    W_l = (%2.3f, %2.3f, %2.3f) ' % (W_l[0],W_l[1],W_l[2]) )
print('    W_r = (%2.3f, %2.3f, %2.3f) ' % (W_r[0],W_r[1],W_r[2]) )
print(' ')

# stupid initial guess, take the arithmetic mean for p:
#pstar,niters = NewtonSolveScalar( f, fp, 0.5*(p_l+p_r), tol=1e-14 )

# compute ustar, which is only a function of pstar:
#ustar        = 0.5*( u_l + u_r ) + 0.5*( f_K( pstar, W_r ) - f_K( pstar, W_l ) )

#Ws_l = [float('NaN'), ustar, pstar]
#Ws_r = [float('NaN'), ustar, pstar]

#print('Newton solve converged in %d iterations' % niters )
#print(' ')
#print('Pressure and velocity inside the entire the star region:')
#print('    p_star = %2.15e' % (pstar) )
#print('    u_star = %2.15e' % ustar   )
#print(' ')
#
##===============================================================================
## Leftmost characteristic (s^1)
##===============================================================================
#if( pstar > p_l ):
#    print("Left Shock:")
#    shock_speed_left = u_l - c_l*sqrt( (gp1/(2.*gamma))*(pstar/p_l)+gm1/(2*gamma))
#    print("     left shock speed = %f" % shock_speed_left )
#else:
#    print("Left Rarefaction: ")
#
#    rho_star_l = rho_l * ( pstar / p_l )**(1./gamma)
#    Ws_l[0] = rho_star_l
#
#    cs = sqrt( abs( gamma*pstar / rho_star_l ) )
#
#    c_rarefaction_left_l = u_l   - c_l
#    c_rarefaction_left_r = ustar - cs
#    print("    Left  foot speed  = %2.15e" % c_rarefaction_left_l )
#    print("    Right foot speed  = %2.15e" % c_rarefaction_left_r )
#
#    # This is the correct value when the left side is a rarefaction
#    print("    rho_star_left     = %2.15e " % rho_star_l )
#
##===============================================================================
## Rightmost characteristic (s^3)
##===============================================================================
#if( pstar > p_r ):
#    print("Right Shock:")
#    shock_speed_right = u_r + c_r*sqrt( (gp1/(2.*gamma))*(pstar/p_r)+gm1/(2*gamma))
#    print("    Right shock speed = %2.15e" % shock_speed_right )
#
#    # (c.f. equation 4.57, when s^3 produces a shock)
#    rho_star_r = \
#        rho_r * ( (pstar/p_r) + (gm1/gp1) ) / ( (gm1/gp1)* pstar/ p_r + 1.0 ) 
#    print("    rho_star_right    = %2.15e " % rho_star_r )
#
#    Ws_r[0] = rho_star_r
#
#else:
#    print("Right Rarefaction " )
#
#print(' ')
#print("In summary, we have the following values for each region:")
#print('    Region 1: W = ', W_l )
#print('    Region 2: W = ', Ws_l )
#print('    Region 3: W = ', Ws_r )
#print('    Region 4: W = ', W_r )
#
#print("If you saw NaN, that means we haven't included that option yet!")
##===============================================================================
