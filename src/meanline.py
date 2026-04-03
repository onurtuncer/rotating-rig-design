# ==============================================================
#  meanline.py  —  Single-stage rotor-only meanline analysis
#
#  Original module from the rotating-rig-design repo.
#  Provides the parametric design-space sweep used in the
#  meanline notebooks before IGV / detailed design modules
#  take over.
#
#  See also: igv.py / meanline_with_igv() for the updated
#  analysis that includes IGV pre-swirl.
# ==============================================================

import numpy as np
from src.constants import gamma, R, Cp, T0_in, P0_in


def omega_from_RPM(N):
    return 2.0 * np.pi * N / 60.0


def meanline_analysis(D_tip, N, PR, eta_is, nu=0.70, phi=0.15, h_min=0.05):
    """
    Single-stage rotor meanline analysis at mean radius.
    Zero inlet swirl (no IGV pre-swirl).

    Parameters
    ----------
    D_tip   : float  tip diameter [m]
    N       : float  rotational speed [RPM]
    PR      : float  total-to-total pressure ratio [-]
    eta_is  : float  isentropic efficiency [-]
    nu      : float  hub-to-tip radius ratio [-]
    phi     : float  flow coefficient Ca/U_mean [-]
    h_min   : float  minimum acceptable blade height [m]

    Returns
    -------
    dict with geometry, velocities, thermodynamics, and
    dimensionless aerodynamic coefficients.
    """
    # Geometry
    r_tip  = D_tip / 2.0
    r_hub  = nu * r_tip
    r_mean = 0.5 * (r_tip + r_hub)
    h      = r_tip - r_hub
    A_flow = np.pi * (r_tip**2 - r_hub**2)

    # Velocities
    omega  = omega_from_RPM(N)
    U_mean = omega * r_mean
    U_tip  = omega * r_tip
    M_tip  = U_tip / np.sqrt(gamma * R * T0_in)
    Ca     = phi * U_mean

    # Thermodynamics
    T0_ratio    = PR**((gamma - 1.0) / gamma)
    delta_T0_is = T0_in * (T0_ratio - 1.0)
    delta_T0    = delta_T0_is / eta_is
    W_euler     = Cp * delta_T0

    # Dimensionless coefficients
    psi      = W_euler / U_mean**2
    C_theta2 = W_euler / U_mean

    # Physical validity: C_theta2 must be less than U_mean (psi < 1.0)
    physically_valid = bool(C_theta2 < U_mean)

    # De Haller number
    W1_mag    = np.sqrt(Ca**2 + U_mean**2)
    W2_mag    = np.sqrt(Ca**2 + (U_mean - C_theta2)**2)
    De_Haller = W2_mag / W1_mag

    # Mass flow and shaft power
    rho_est = P0_in / (R * T0_in)
    mdot    = rho_est * Ca * A_flow
    P_shaft = mdot * W_euler / 1000.0

    # Blade angles
    beta1      = np.degrees(np.arctan(U_mean / Ca))
    beta2      = np.degrees(np.arctan((U_mean - C_theta2) / Ca))
    alpha2     = np.degrees(np.arctan(C_theta2 / Ca))
    delta_beta = beta1 - beta2

    return {
        'D_tip_mm':        D_tip * 1000,
        'r_tip':           r_tip,
        'r_hub':           r_hub,
        'r_mean':          r_mean,
        'h_mm':            h * 1000,
        'h_ok':            h >= h_min,
        'physically_valid': physically_valid,
        'nu':              nu,
        'A_flow_m2':       A_flow,
        'N_RPM':           N,
        'omega':           omega,
        'U_mean':          U_mean,
        'U_tip':           U_tip,
        'M_tip':           M_tip,
        'PR':              PR,
        'eta_is':          eta_is,
        'delta_T0':        delta_T0,
        'delta_T0_is':     delta_T0_is,
        'W_euler_kJ':      W_euler / 1000.0,
        'psi':             psi,
        'phi':             phi,
        'Ca':              Ca,
        'De_Haller':       De_Haller,
        'mdot_kg_s':       mdot,
        'P_shaft_kW':      P_shaft,
        'beta1_deg':       beta1,
        'beta2_deg':       beta2,
        'alpha2_deg':      alpha2,
        'delta_beta_deg':  delta_beta,
        'C_theta2':        C_theta2,
    }
