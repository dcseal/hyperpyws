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
W_l = [3.0, 0.0   ]
W_r = [1.0, 0.0   ]

h_l = W_l[0]; u_l = W_l[1]
h_r = W_r[0]; u_r = W_r[1]

#c_l, c_r = flux.eig
c_l = h_l - sqrt( g*h_l )
c_r = h_r + sqrt( g*h_r )

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
def phip_l( h ):
    """ Derivative of the function, phi_l."""

    if( h < h_l ):
        return - sqrt( g/h )

    else:
        return -sqrt(2.)*g*( h*h_l + 2.*h**2 + h_l**2 ) #/ ( 4.*h_l * sqrt( g*(1./h + 1./h_l) ) * h**2 )


#===============================================================================
def phi_r( h ):
    """ These function calls look identical for left/right values.

    See section 13.10, ``Finite Volume Methods for Hyperbolic Problems",
    LeVeque.
    """
    if( h < h_r ):
        return u_r - 2.0*( sqrt(g*h_r) - sqrt(g*h) )
    else:
        return u_r + (h - h_r) * sqrt( 0.5*g*( 1./h + 1./h_r ) )

#===============================================================================
def phip_r( h ):
    """ Derivative of the function, phi_r."""

    if( h < h_l ):
        return sqrt( g/h )

    else:
        return sqrt(2.)*g*( h*h_r + 2.*h**2 + h_r**2 ) / ( 4.*h_r * sqrt( g*(1./h + 1./h_r) ) * h**2 )


#===============================================================================
def phi( h ):
    """ hstar is found by solving phi(hstar) = 0. """
    return phi_l( h ) - phi_r( h )

#===============================================================================
def phi_p( h ):
    """ Derivative of phi from above. """
    return phip_l( h ) - phip_r( h )

#===============================================================================
def RightFan( x, t, W_r ):
    """ Data for values inside the RightFan.

    User is responsible for knowing correct bounds of x to supply.
    """

    print("junk fan")

#===============================================================================
def LeftFan( x, t, xs, W ):
    """ Data for values inside the RightFan.

    User is responsible for knowing correct bounds of x to supply.
    """

    u_K = W[0]
    h_K = W[1]

    A = u_K + 2.*sqrt( g*h_K )

    h = 1./(9.*g)*( (u_K + 2.*sqrt(g*h_K)) - (x-xc)/t )**2
    u = u_K - 2.*( sqrt(g*h_K) - sqrt(g*h_K) )
    return h,u


#===============================================================================
# Print (and compute) information to std output
#===============================================================================
print('Initial data:')
print('    W_l = (%2.3f, %2.3f) ' % (W_l[0],W_l[1] ) )
print('    W_r = (%2.3f, %2.3f) ' % (W_r[0],W_r[1] ) )
print(' ')

# stupid initial guess, take the arithmetic mean for p:
hstar, niters = NewtonSolveScalar( phi, phi_p, 0.5*(h_l+h_r), tol=1e-14 )

# TODO: this is correct only in the case of a 1-rarefaction:
# (or is it?)
# ( eqn: (13.57), LeVeque )
ustar = u_l + 2.*( sqrt(g*h_l) - sqrt( g*hstar ) )

Ws = [ustar, hstar]

print('Newton solve converged in %d iterations' % niters )
print(' ')
print('    hstar = %2.15e' % hstar   )
print('    ustar = %2.15e' % ustar   )
print(' ')

#print('ustar + sqrt(g*hstar) = %f' % (ustar + sqrt(g*hstar) ) )
#
##===============================================================================
## Left characteristic (s^1)
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
## Right characteristic (s^2)
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
