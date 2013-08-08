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

try               :  import hyperpyws
except ImportError:  import hyperpyws_path

#===============================================================================
# FLUX FUNCTION
#===============================================================================

from hyperpyws.model_equations.buckley_leverett  import BuckleyLeverett1D
from hyperpyws.iterative_solvers import SecantSolveScalar

M    = 1.0/3.0
flux = BuckleyLeverett1D( M )
f    = flux.f
fp   = flux.eig

#-----------------------------------------------------------------------------
def main():
    """Single function to compute the two shock speeds that are present for
    our numerical example.
   
    This has not been written in a general fashion - certainly we could
    rewrite this so it could be used again in the future.
    """

    # Riemman state that remains static:
    q_fixed  = 1.
    f_qfixed = f( [q_fixed] )[0]
    def G(qs):
        """ zero function that needs to be solved for left-state of our
        example Riemman problem

        R-H conditions:

            f'(qs) = (f(qs) - f(q_fixed)) / (q - q_fixed).

        """
        return ( fp( [qs] )[0] - ( f([qs])[0] - f_qfixed ) / (qs-q_fixed) )


    qs_left_problem = SecantSolveScalar( G, 0.13, 0.130000001 )

    print('left  value qs = %2.15e, speed = %2.15e' % 
        (qs_left_problem, fp([qs_left_problem])[0]  ) )

    q_fixed  = 0.
    f_qfixed = f([q_fixed])[0]

    qs_right_problem = SecantSolveScalar( G, 0.48, 0.480000001 )
    print('right value qs = %2.15e, speed = %2.15e' % 
        (qs_right_problem, fp( [qs_right_problem])[0] ) )

    return( qs_left_problem, qs_right_problem )
