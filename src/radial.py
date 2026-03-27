# ==============================================================
#  radial.py — Free vortex radial distribution
#
#  References:
#    Dixon & Hall (2014), Chapter 6
#    Lieblein et al. (1953), NACA RM E53D01
# ==============================================================

import numpy as np


def free_vortex(res, n_stations=50, sigma=1.0):
    """
    Compute radial distribution of blade angles and aerodynamic
    parameters using the free vortex law: r · Cθ = constant.

    Parameters
    ----------
    res         : dict  — Output of meanline_analysis()
    n_stations  : int   — Number of radial stations (hub to tip)
    sigma       : float — Assumed solidity for Lieblein DF [-]

    Returns
    -------
    dict — Arrays of radial quantities (length = n_stations)
    """
    rm   = res['r_mean']
    rt   = res['r_tip']
    rh   = res['r_hub']
    Ca   = res['Ca']
    Ct2m = res['C_theta2']
    omega= res['omega']

    r        = np.linspace(rh, rt, n_stations)
    pct_span = (r - rh) / (rt - rh) * 100

    # Free vortex: Cθ₂(r) = Cθ₂,mean · r_mean / r
    Ct2 = Ct2m * rm / r
    U_r = omega * r

    # Blade angles
    beta1      = np.degrees(np.arctan(U_r / Ca))
    beta2      = np.degrees(np.arctan((U_r - Ct2) / Ca))
    alpha2     = np.degrees(np.arctan(Ct2 / Ca))
    delta_beta = beta1 - beta2

    # De Haller
    W1 = np.sqrt(Ca**2 + U_r**2)
    W2 = np.sqrt(Ca**2 + (U_r - Ct2)**2)
    DH = W2 / W1

    # Lieblein diffusion factor
    DF = 1.0 - DH + Ct2 / (2.0 * sigma * W1)

    return {
        'r':          r,
        'pct_span':   pct_span,
        'U_r':        U_r,
        'Ct2':        Ct2,
        'beta1':      beta1,
        'beta2':      beta2,
        'alpha2':     alpha2,
        'delta_beta': delta_beta,
        'W1':         W1,
        'W2':         W2,
        'DH':         DH,
        'DF':         DF,
        'sigma':      sigma,
    }