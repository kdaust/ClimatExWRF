#### CODE BASED ON GEMINI OUTPUT BUT ADAPTED ####
#### SHOULD DOUBLE CHECK ####

# Computes both IPW and IVT 

# Ingests qv, p, pb, psfc, ua (from wrf.getvar), and va (from wrf.getvar)

def compute_ipw_ivt(qv, p, pb, psfc, ua, va):
    import numpy as np
    
    # Define paths and target
    TARGET_P = 50000.0  # 50 kPa in Pascals
    G = 9.80665  # Gravity (m/s^2)

    # 1. Compute total full pressure at mass points
    p_tot = p + pb
    nz, ny, nx = p_tot.shape

    # 2. Reconstruct layer interfaces (staggered vertical boundary pressures)
    p_interface = np.zeros((nz + 1, ny, nx))
    p_interface[0, :, :] = psfc  # Bottom interface is the surface pressure

    # Intermediate interfaces are calculated using the midpoints of mass levels
    p_interface[1:-1, :, :] = 0.5 * (p_tot[:-1, :, :] + p_tot[1:, :, :])

    # Top interface (extrapolated or top model layer boundary approximation)
    p_interface[-1, :, :] = p_tot[-1, :, :] - (p_interface[-2, :, :] - p_tot[-1, :, :])


    # 3. Interpolate the exact qv at the 50 kPa target surface
    # We will use vectorization across the horizontal grid to locate the crossover layer
    qv_at_target = np.zeros((ny, nx))
    ua_at_target = np.zeros((ny, nx))
    va_at_target = np.zeros((ny, nx))

    below_target = p_tot < TARGET_P
    has_crossing = below_target.any(axis=0)
    # argmax returns zero for all-False columns: exclude those explicitly.
    z_above = below_target.argmax(axis=0)
    first_level = has_crossing & (z_above == 0)
    qv_at_target[first_level] = qv[0][first_level]
    ua_at_target[first_level] = ua[0][first_level]
    va_at_target[first_level] = va[0][first_level]

    y, x = np.nonzero(has_crossing & (z_above > 0))
    above = z_above[y, x]
    below = above - 1
    log_b = np.log(p_tot[below, y, x])
    log_a = np.log(p_tot[above, y, x])
    weight = (np.log(TARGET_P) - log_b) / (log_a - log_b)
    weight_above = (log_a - np.log(TARGET_P)) / (log_a - log_b)

    qv_at_target[y, x] = qv[below, y, x] + weight * (
        qv[above, y, x] - qv[below, y, x]
    )
    # Preserve the original wind weights; changing the integration method is
    # separate from vectorising it.
    ua_at_target[y, x] = weight * ua[below, y, x] + weight_above * ua[above, y, x]
    va_at_target[y, x] = weight * va[below, y, x] + weight_above * va[above, y, x]


    # 4. Integrate layer by layer up to the exact 50 kPa boundary
    ipw_50kPa = np.zeros((ny, nx))
    ivtx_50kPa = np.zeros((ny, nx))
    ivty_50kPa = np.zeros((ny, nx))

    # Layer by layer
    for z in range(nz):
        # Determine the boundaries of the current layer
        p_bot = p_interface[z, :, :] # bottom interface
        p_top = p_interface[z + 1, :, :] # top interface

        # Conditions for integration:
        # Case A: Entire layer is below 50 kPa boundary (Pressure > 50000 Pa)
        mask_full_layer = p_top >= TARGET_P

        # Case B: Layer straddles the 50 kPa boundary (p_bot > 50000 >= p_top)
        mask_partial_layer = (p_bot > TARGET_P) & (p_top < TARGET_P)

        # Calculate dP (thickness in Pa) for full layers
        dp = np.zeros((ny, nx))
        dp[mask_full_layer] = p_bot[mask_full_layer] - p_top[mask_full_layer]

        # Calculate truncated dP for the boundary-straddling layer
        dp[mask_partial_layer] = p_bot[mask_partial_layer] - TARGET_P

        # Assign correct qv profile
        qv_layer = np.zeros((ny, nx))
        ua_layer = np.zeros((ny, nx))
        va_layer = np.zeros((ny, nx))
        # Full layers take their standard mixing ratio and velocities
        qv_layer[mask_full_layer] = qv[z, mask_full_layer]
        ua_layer[mask_full_layer] = ua[z, mask_full_layer]
        va_layer[mask_full_layer] = va[z, mask_full_layer]


        # Straddling layers take average between layer bottom qv and interpolated target qv
        qv_layer[mask_partial_layer] = 0.5 * (
            qv[z, mask_partial_layer] + qv_at_target[mask_partial_layer]
        )
        
        # Straddling layers take average between layer bottom ua and interpolated target ua
        ua_layer[mask_partial_layer] = 0.5 * (
            ua[z, mask_partial_layer] + ua_at_target[mask_partial_layer]
        )

        # Straddling layers take average between layer bottom va and interpolated target va
        va_layer[mask_partial_layer] = 0.5 * (
            va[z, mask_partial_layer] + va_at_target[mask_partial_layer]
        )


        # Accumulate water mass per unit area (kg/m^2 which is equivalent to mm)
        ipw_50kPa += (qv_layer * dp) / G

        # Accumulate IVTx
        ivtx_50kPa += (qv_layer * ua_layer * dp) / G
        # Accumulate IVTy
        ivty_50kPa += (qv_layer * va_layer * dp) / G


    # Compute IVT magnitude
    ivt_50kPa = np.sqrt(ivtx_50kPa**2 + ivty_50kPa**2)

    """
    print("Integration complete.")
    print(f"Min PW to 50kPa: {pw_50kPa.min():.2f} mm")
    print(f"Max PW to 50kPa: {pw_50kPa.max():.2f} mm")

    #print(f"Min PW total: {pw_getvar.min():.2f} mm")
    #print(f"Max PW total: {pw_getvar.max():.2f} mm")

    print(f"Min IVT to 50kPa: {ivt_50kPa.min():.2f} kg m-1 s-1")
    print(f"Max IVT to 50kPa: {ivt_50kPa.max():.2f} kg m-1 s-1")
    """

    return(ipw_50kPa,ivt_50kPa,ivtx_50kPa,ivty_50kPa)
 
if __name__ == '__main__':
    pass

