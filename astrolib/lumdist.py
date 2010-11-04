"""
 NAME:
    LUMDIST

 PURPOSE:
    Calculate luminosity distance (in Mpc) of an object given its redshift
 EXPLANATION:
    The luminosity distance in the Friedmann-Robertson-Walker model is
    taken from  Caroll, Press, and Turner (1992, ARAA, 30, 499), p. 511
    Uses a closed form (Mattig equation) to compute the distance when the
    cosmological constant is zero.   Otherwise integrates the function using
    QSIMP.
 CALLING SEQUENCE:
    result = lumdist(z, [H0 = , k = , Omega_M =, Lambda0 = , q0 = ,/SILENT])

 INPUTS:
    z = redshift, positive scalar or vector

 OPTIONAL KEYWORD INPUTS:
    /SILENT - If set, the program will not display adopted cosmological
        parameters at the terminal.
    H0: Hubble parameter  in km/s/Mpc, default is 70

        No more than two of the following four parameters should be
        specified.   None of them need be specified -- the adopted defaults
        are given.
    k - curvature constant, normalized to the closure density.   Default is
        0, indicating a flat universe
    Omega_m -  Matter density, normalized to the closure density, default
        is 0.3.   Must be non-negative
    Lambda0 - Cosmological constant, normalized to the closure density,
        default is 0.7
    q0 - Deceleration parameter, numeric scalar = -R*(R'')/(R')^2, default
        is -0.55

 OUTPUTS:
    The result of the function is the luminosity distance (in Mpc) for each
    input value of z.

 EXAMPLE:
    (1) Plot the distance of a galaxy in Mpc as a function of redshift out
        to z = 5.0, assuming the default cosmology (Omega_m=0.3, Lambda = 0.7,
        H0 = 70 km/s/Mpc)

        IDL> z = findgen(50)/10.
        IDL> plot,z,lumdist(z),xtit='z',ytit='Distance (Mpc)'

        Now overplot the relation for zero cosmological constant and
        Omega_m=0.3
        IDL> oplot,z,lumdist(z,lambda=0,omega=0.3),linestyle=1
 COMMENTS:
    (1) Integrates using the IDL Astronomy Version procedure QSIMP.    (The
    intrinsic IDL QSIMP function is not called because of its ridiculous
    restriction that only scalar arguments can be passed to the integrating
    function.)
    (2) Can fail to converge at high redshift for closed universes with
    non-zero lambda.   This can presumably be fixed by replacing QSIMP with
    an integrator that can handle a singularity
 PROCEDURES CALLED:
    COSMO_PARAM, QSIMP
 REVISION HISTORY:
    Written   W. Landsman        Raytheon ITSS       April 2000
    Avoid integer overflow for more than 32767 redshifts  July 2001
    Use double precision J. Moustakas/W. Landsman April 2008
"""

from numpy import *
from cosmo_param import cosmo_param
import scipy.integrate

def ldist(z, q0=None, lambda0=None):
   
   term1 = (1. + z) ** 2
   term2 = 1. + 2. * (q0 + lambda0) * z
   term3 = z * (2. + z) * lambda0
   denom = (term1 * term2 - term3)
   out = z * 0.
   good = denom > 0.0
   if isinstance(z,ndarray):
      out[good] = 1. / sqrt(denom[good])
   else:
      if good:
         out = 1. / sqrt(denom)
      
   return out


def lumdist(z, h0=None, k=None, lambda0=None, omega_m=None, q0=None, silent=None):

   n_params = 1
   
   # COMPILE_OPT IDL2
   if n_params == 0:   
      print 'Syntax: result = lumdist(z, H0 = ,k=, Lambda0 = ])'
      print 'Returns luminosity distance in Mpc'
      return 0.
   if isinstance(z,ndarray) or isinstance(z,list):
      n = len(array(z))
      z=numpy.array(z)
   else:
      z=array([z])
      n=len(z)
      print z

   omega_m,lambda0,k,q0=cosmo_param(omega_m, lambda0, k, q0)
   
   # Check keywords
   c = 2.99792458e5                  #  speed of light in km/s
   if h0 is None:
      h0 = 70
   if not silent:
      print 'LUMDIST: H0:', h0, ' Omega_m:', omega_m, ' Lambda0', lambda0, ' q0: ', q0, ' k: ', k#, format='(A,I3,A,f5.2,A,f5.2,A,f5.2,A,F5.2)'
   
   # For the case of Lambda = 0, we use the closed form from equation 5.238 of
   # Astrophysical Formulae (Lang 1998).   This avoids terms that almost cancel
   # at small q0*z better than the more familiar Mattig formula.
   #
   if lambda0 == 0:   
      denom = sqrt(1 + 2 * q0 * z) + 1 + q0 * z
      dlum = (c * z / h0) * (1 + z * (1 - q0) / denom)
      return dlum
      
      # For non-zero lambda
   else:   
      dlum = z * 0.0
      for i in arange(0, (n - 1)+(1)):
         if z[i] <= 0.0:   
            dlum[i] = 0.0
         else:   
            #qsimp('LDIST', 0, z[i], lz, q0=q0, lambda0=lambda0)
            lz=scipy.integrate.quad(ldist,0,z[i],args=(q0,lambda0))
            dlum[i] = lz[0]
      
      if k > 0:   
         dlum = sinh(sqrt(k) * dlum) / sqrt(k)
      else:   
         if k < 0:   
            dlum = maximum(sin(sqrt(-k) * dlum) / sqrt(-k), 0)
      return c * (1 + z) * dlum / h0
   

