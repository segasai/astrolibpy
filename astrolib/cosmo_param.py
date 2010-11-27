def cosmo_param(omega_m=None, omega_lambda=None, omega_k=None, q0=None):
   """
    NAME:
        COSMO_PARAM
    PURPOSE:
        Derive full set of cosmological density parameters from a partial set
    EXPLANATION:
        This procedure is called by LUMDIST and GALAGE to allow the user a choice
        in defining any two of four cosmological density parameters.
   
        Given any two of the four input parameters -- (1) the normalized matter
        density Omega_m (2) the normalized cosmological constant, Omega_lambda
        (3) the normalized curvature term, Omega_k and (4) the deceleration
        parameter q0 --  this  program will derive the remaining two.     Here
        "normalized" means divided by the closure density so that
        Omega_m + Omega_lambda + Omega_k = 1.    For a more
        precise definition see Carroll, Press, & Turner (1992, ArAA, 30, 499).
   
        If less than two parameters are defined, this procedure sets default
        values of Omega_k=0 (flat space), Omega_lambda = 0.7, Omega_m = 0.3
        and q0 = -0.55
    CALLING SEQUENCE:
          COSMO_PARAM, Omega_m, Omega_lambda, Omega_k, q0
   
    INPUT-OUTPUTS:
        Omega_M - normalized matter energy density, non-negative numeric scalar
        Omega_Lambda - Normalized cosmological constant, numeric scalar
        Omega_k - normalized curvature parameter, numeric scalar.   This is zero
                  for a flat universe
        q0 - Deceleration parameter, numeric scalar = -R*(R'')/(R')^2
             = 0.5*Omega_m - Omega_lambda
    NOTES:
        If more than two parameters are defined upon input (overspecification),
        then the first two defined parameters in the ordered list Omega_m,
        Omega_lambda, Omega_k, q0 are used to define the cosmology.
    EXAMPLE:
        Suppose one has Omega_m = 0.3, and Omega_k = 0.5 then to determine
        Omega_lambda and q0
   
          IDL> cosmo_param, 0.3, omega_lambda, 0.5, q0
   
          which will return omega_lambda = 0.2 and q0 = -2.45
    REVISION HISTORY:
          W. Landsman         Raytheon ITSS         April 2000
   """

   nk = omega_k is not None
   nlambda = omega_lambda is not None 
   nomega = omega_m is not None
   nq0 = q0 is not None
   # Check which two parameters are defined, and then determine the other two
   
   if nomega and nlambda:
      if not nk:   
         omega_k = 1 - omega_m - omega_lambda
      if not nq0:   
         q0 = omega_m / 2. - omega_lambda
   
   if nomega and nk:
      if not nlambda:   
         omega_lambda = 1. - omega_m - omega_k
      if not nq0:   
         q0 = -1 + omega_k + 3 * omega_m / 2

   if nlambda and nk:   
      if not nomega:   
         omega_m = 1. - omega_lambda - omega_k
      if not nq0:   
         q0 = (1 - omega_k - 3. * omega_lambda) / 2.

   if nomega and nq0:   
      if not nk:   
         omega_k = 1 + q0 - 3 * omega_m / 2.
      if not nlambda:   
         omega_lambda = 1. - omega_m - omega_k

   if nlambda and nq0:   
      if not nk:   
         omega_k = 1 - 2 * q0 - 3 * omega_lambda
      if not nomega:   
         omega_m = 1. - omega_lambda - omega_k
   
   if nk and nq0:   
      if not nomega:   
         omega_m = (1 + q0 - omega_k) * 2 / 3.
      if not nlambda:   
         omega_lambda = 1. - omega_m - omega_k
   
   #Set default values
   if omega_k is None:
      omega_k = 0       #Default is flat space
   if omega_lambda is None:
      omega_lambda = 0.7
   if omega_m is None:
      omega_m = 1 - omega_lambda
   if q0 is None:
      q0 = (1 - omega_k - 3 * omega_lambda) / 2.
   
   return omega_m, omega_lambda, omega_k, q0

