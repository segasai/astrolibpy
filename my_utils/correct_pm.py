import numpy as np
import astropy.coordinates as acoo
import astropy.units as auni

# Use v4.0 defaults
GCPARAMS = acoo.galactocentric_frame_defaults.get_from_registry(
    "v4.0")['parameters']

kms = auni.km / auni.s
masyr = auni.mas / auni.year
kpc = auni.kpc

# --- Pre-calculate constants and solar velocity vector ---

# Conversion factor: (km/s) / (kpc * mas/yr)
# k_masyr_kpc = (kms / (kpc * masyr)).to_value(1) # ~4.74047
k_masyr_kpc = 4.740470463498976

# Get the Sun's velocity vector relative to the GC in Galactocentric Cartesian
# coordinates (U, V, W)
v_sun_gal_xyz = GCPARAMS['galcen_v_sun'].get_d_xyz().to_value(
    kms)  # [U, V, W] in km/s


# We need the Sun's velocity vector expressed in ICRS Cartesian coordinates
# (vx, vy, vz).
# The original Astropy code effectively calculates the apparent motion seen
# by the Sun
# due to an object at rest in the Galactocentric frame. This corresponds to
# transforming the *negative* of the Sun's velocity vector (-U, -V, -W) from
# the Galactocentric frame to the ICRS frame.
# We can get this precise vector by performing the Astropy calculation once.
def _get_apparent_solar_velocity_icrs(params):
    """
    Helper to get the ICRS velocity vector corresponding to solar motion
    relative to the GC frame origin, consistent with the original code's logic.
    """
    # Use a dummy coordinate transformation
    # A point at rest (0 velocity) in Galactocentric frame
    gc_rest = acoo.Galactocentric(
        x=0 * kpc,
        y=0 * kpc,
        z=0 * kpc,  # Position doesn't affect velocity transform
        v_x=0 * kms,
        v_y=0 * kms,
        v_z=0 * kms,
        **params)  # Pass GCPARAMS for frame definition
    # Transform this point *to* ICRS. Its velocity in ICRS represents the
    # apparent velocity induced by the Sun's motion relative to the GC frame.
    icrs_apparent = gc_rest.transform_to(acoo.ICRS())
    # Return the Cartesian velocity components in ICRS
    return icrs_apparent.velocity.d_xyz.to_value(kms)


# Calculate this vector once using the specified GCPARAMS
V_APPARENT_XYZ_ICRS = _get_apparent_solar_velocity_icrs(GCPARAMS)
vx, vy, vz = V_APPARENT_XYZ_ICRS  # km/s


def correct_vel(ra, dec, vel, vlsr=None):
    """
    Corrects the radial velocity for the speed of the Sun using direct
    calculation.

    Arguments:
        ra - RA in deg (array or scalar)
        dec -- Declination in deg (array or scalar)
        vel -- heliocentric rv in km/s (array or scalar)
        vlsr (optional) -- ignored, kept for interface compatibility
        split (optional) -- ignored, calculations are vectorized

    Returns:
        radial velocity corrected for solar reflex motion (km/s)
    """
    if vlsr is not None:
        print('WARNING vlsr is ignored')

    # Ensure inputs are numpy arrays for vectorized operations
    is_scalar = np.isscalar(ra) and np.isscalar(dec) and np.isscalar(vel)
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)
    vel = np.atleast_1d(vel)

    # Convert angles to radians
    ra_rad = np.deg2rad(ra)
    dec_rad = np.deg2rad(dec)

    # Calculate unit vector components in line-of-sight direction (p_hat)
    cos_dec = np.cos(dec_rad)
    sin_dec = np.sin(dec_rad)
    del dec_rad
    cos_ra = np.cos(ra_rad)
    sin_ra = np.sin(ra_rad)
    del ra_rad
    # Unit vector p_hat = (cos(dec)cos(ra), cos(dec)sin(ra), sin(dec))
    px = cos_dec * cos_ra
    py = cos_dec * sin_ra
    pz = sin_dec
    del cos_dec, sin_dec, cos_ra, sin_ra
    # Project apparent solar velocity (V_APPARENT_XYZ_ICRS) onto the line
    # of sight (p_hat)
    # This gives the radial velocity component *caused* by solar motion
    # (like C1.rv in original)
    # delta_rv = vx * px + vy * py + vz * pz  # km/s
    # Use einsum for potentially better performance/clarity with large arrays
    delta_rv = np.dot(V_APPARENT_XYZ_ICRS,
                      np.stack([px, py, pz], axis=0))  # CORRECTED using np.dot
    # Corrected velocity = Observed velocity - component due to Sun's motion
    corrected_vel = vel - delta_rv

    # Return scalar if input was scalar
    if is_scalar:
        return corrected_vel.item()
    return corrected_vel


def correct_pm(ra, dec, pmra, pmdec, dist, vlsr=None):
    """
    Corrects the proper motion for the speed of the Sun using direct
    calculation.

    Arguments:
        ra - RA in deg (array or scalar)
        dec -- Declination in deg (array or scalar)
        pmra -- pm in RA in mas/yr (with cosine term) (array or scalar)
        pmdec -- pm in declination in mas/yr (array or scalar)
        dist -- distance in kpc (array or scalar)
        split (optional) -- ignored, calculations are vectorized
        vlsr (optional) -- ignored, kept for interface compatibility

    Returns:
        (pmra, pmdec) tuple with proper motions corrected for Sun's motion
        (mas/yr)
    """
    if vlsr is not None:
        print('WARNING vlsr is ignored')

    # Ensure inputs are numpy arrays
    is_scalar = (np.isscalar(ra) and np.isscalar(dec) and np.isscalar(pmra)
                 and np.isscalar(pmdec) and np.isscalar(dist))
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)
    pmra = np.atleast_1d(pmra)
    pmdec = np.atleast_1d(pmdec)
    dist = np.atleast_1d(dist)

    # Handle case where dist is scalar but others are arrays
    if dist.size == 1 and ra.size > 1:
        dist = np.full_like(ra, dist.item())
    elif dist.size != ra.size:
        raise ValueError("Shape mismatch between coordinates and distance")

    # Convert angles to radians
    ra_rad = np.deg2rad(ra)
    dec_rad = np.deg2rad(dec)

    # Calculate unit vector components for RA and Dec directions
    sin_ra = np.sin(ra_rad)
    cos_ra = np.cos(ra_rad)
    sin_dec = np.sin(dec_rad)
    cos_dec = np.cos(dec_rad)
    del ra_rad, dec_rad
    # RA direction unit vector (ra_hat = (-sin(ra), cos(ra), 0))
    ra_hat_x = -sin_ra
    ra_hat_y = cos_ra
    ra_hat_z = np.zeros_like(sin_ra)  # Explicitly zero for clarity

    # Dec direction unit vector
    # (dec_hat = (-sin(dec)cos(ra), -sin(dec)sin(ra), cos(dec)))
    dec_hat_x = -sin_dec * cos_ra
    dec_hat_y = -sin_dec * sin_ra
    dec_hat_z = cos_dec
    del sin_ra, cos_ra, cos_dec, sin_dec
    # Project apparent solar velocity (V_APPARENT_XYZ_ICRS) onto RA
    # and Dec directions
    # This gives the tangential velocity components *caused* by solar motion
    # v_ra = vx * ra_hat_x + vy * ra_hat_y + vz * ra_hat_z
    # v_dec = vx * dec_hat_x + vy * dec_hat_y + vz * dec_hat_z
    v_ra = np.dot(V_APPARENT_XYZ_ICRS,
                  np.stack([ra_hat_x, ra_hat_y, ra_hat_z],
                           axis=0))  # CORRECTED using np.dot
    v_dec = np.dot(V_APPARENT_XYZ_ICRS,
                   np.stack([dec_hat_x, dec_hat_y, dec_hat_z],
                            axis=0))  # CORRECTED using np.dot
    # Convert tangential velocities (km/s) to proper motions (mas/yr)
    # delta_pm = v_tangential / (distance_kpc * k_masyr_kpc)
    # Avoid division by zero or very small distances if necessary
    with np.errstate(divide='ignore',
                     invalid='ignore'):  # Suppress warnings for dist=0
        delta_pmra_cosdec = v_ra / (dist * k_masyr_kpc)  # mas/yr
        delta_pmdec = v_dec / (dist * k_masyr_kpc)  # mas/yr

        # Handle potential NaNs/Infs resulting from dist=0 or very small dist
        delta_pmra_cosdec[~np.isfinite(delta_pmra_cosdec
                                       )] = np.nan  # Or np.nan if preferred
        delta_pmdec[~np.isfinite(delta_pmdec
                                 )] = np.nan  # Or np.nan if preferred

    # Corrected PM = Observed PM - component due to Sun's motion
    corrected_pmra = pmra - delta_pmra_cosdec
    corrected_pmdec = pmdec - delta_pmdec

    # Return scalars if input was scalar
    if is_scalar:
        return corrected_pmra.item(), corrected_pmdec.item()
    return corrected_pmra, corrected_pmdec


# --- Example Usage and Comparison (Optional) ---
if __name__ == "__main__":
    # Import original functions for comparison
    from correct_pm import correct_pm as correct_pm_orig
    from correct_pm import correct_vel as correct_vel_orig

    # Sample Data (e.g., 5 stars)
    N_stars = 5000
    ra_deg = np.random.uniform(0, 360, N_stars)
    dec_deg = np.random.uniform(-90, 90, N_stars)
    pmra_masyr = np.random.normal(0, 5, N_stars)  # pm_ra * cos(dec)
    pmdec_masyr = np.random.normal(0, 5, N_stars)
    dist_kpc = np.random.uniform(0.1, 10, N_stars)
    vel_kms = np.random.normal(0, 50, N_stars)

    # --- Test Velocity Correction ---
    print("--- Velocity Correction Test ---")
    corrected_vel_orig = correct_vel_orig(ra_deg, dec_deg, vel_kms)
    corrected_vel_fast = correct_vel(ra_deg, dec_deg, vel_kms)

    print("Original Astropy Corrected Vel:", corrected_vel_orig)
    print("Fast Numpy Corrected Vel:     ", corrected_vel_fast)
    print("Difference:", corrected_vel_orig - corrected_vel_fast)
    print("Max absolute difference:",
          np.max(np.abs(corrected_vel_orig - corrected_vel_fast)))

    # --- Test Proper Motion Correction ---
    print("\n--- Proper Motion Correction Test ---")
    pmra_corr_orig, pmdec_corr_orig = correct_pm_orig(ra_deg, dec_deg,
                                                      pmra_masyr, pmdec_masyr,
                                                      dist_kpc)
    pmra_corr_fast, pmdec_corr_fast = correct_pm(ra_deg, dec_deg, pmra_masyr,
                                                 pmdec_masyr, dist_kpc)

    print("Original Astropy Corrected PMRA:", pmra_corr_orig)
    print("Fast Numpy Corrected PMRA:     ", pmra_corr_fast)
    print("Difference PMRA:", pmra_corr_orig - pmra_corr_fast)
    print("Max absolute difference PMRA:",
          np.max(np.abs(pmra_corr_orig - pmra_corr_fast)))

    print("\nOriginal Astropy Corrected PMDEC:", pmdec_corr_orig)
    print("Fast Numpy Corrected PMDEC:     ", pmdec_corr_fast)
    print("Difference PMDEC:", pmdec_corr_orig - pmdec_corr_fast)
    print("Max absolute difference PMDEC:",
          np.max(np.abs(pmdec_corr_orig - pmdec_corr_fast)))

    # --- Test Scalar Input ---
    print("\n--- Scalar Input Test ---")
    ra_s, dec_s, pmra_s, pmdec_s, dist_s, vel_s = ra_deg[0], dec_deg[
        0], pmra_masyr[0], pmdec_masyr[0], dist_kpc[0], vel_kms[0]

    vel_orig_s = correct_vel_orig(ra_s, dec_s, vel_s)
    vel_fast_s = correct_vel(ra_s, dec_s, vel_s)
    print(f"Scalar Vel Orig: {vel_orig_s:.5f}, Fast: {vel_fast_s:.5f},"
          f"Diff: {vel_orig_s - vel_fast_s:.2e}")

    pm_orig_s = correct_pm_orig(ra_s, dec_s, pmra_s, pmdec_s, dist_s)
    pm_fast_s = correct_pm(ra_s, dec_s, pmra_s, pmdec_s, dist_s)
    print(f"Scalar PMRA Orig: {pm_orig_s[0]:.5f}, Fast: {pm_fast_s[0]:.5f},"
          f" Diff: {pm_orig_s[0] - pm_fast_s[0]:.2e}")
    print(f"Scalar PMDEC Orig: {pm_orig_s[1]:.5f}, Fast: {pm_fast_s[1]:.5f},"
          f" Diff: {pm_orig_s[1] - pm_fast_s[1]:.2e}")

    # --- Test Scalar Input2 ---
    print("\n--- Scalar Input Test ---")
    ra_s, dec_s, pmra_s, pmdec_s, dist_s, vel_s = ra_deg[0], dec_deg[
        0], pmra_masyr[0], pmdec_masyr[0], dist_kpc[0], vel_kms[0]

    pmra_corr_orig, pmdec_corr_orig = correct_pm_orig(ra_deg, dec_deg,
                                                      pmra_masyr, pmdec_masyr,
                                                      dist_s)
    pmra_corr_fast, pmdec_corr_fast = correct_pm(ra_deg, dec_deg, pmra_masyr,
                                                 pmdec_masyr, dist_s)
    print("Original Astropy Corrected PMRA:", pmra_corr_orig)
    print("Fast Numpy Corrected PMRA:     ", pmra_corr_fast)
    print("Difference PMRA:", pmra_corr_orig - pmra_corr_fast)
    print("Max absolute difference PMRA:",
          np.max(np.abs(pmra_corr_orig - pmra_corr_fast)))

    print("\nOriginal Astropy Corrected PMDEC:", pmdec_corr_orig)
    print("Fast Numpy Corrected PMDEC:     ", pmdec_corr_fast)
    print("Difference PMDEC:", pmdec_corr_orig - pmdec_corr_fast)
    print("Max absolute difference PMDEC:",
          np.max(np.abs(pmdec_corr_orig - pmdec_corr_fast)))
