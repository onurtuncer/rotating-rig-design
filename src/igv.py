# ==============================================================
#  igv.py  —  Inlet Guide Vane aerodynamic design
#
#  References:
#    Tian et al. (2018), Experiments in Fluids 59:63
#      SJTU single-stage axial compressor rig
#      B=21 rotor blades (NACA 65), V=13 IGV blades
#      D_tip=598 mm, nu=0.70, N=3000 rpm, phi_des=0.151
#      IGV-rotor axial gap = 1.10 * rotor chord (mid-span)
#
#    Dixon & Hall (2014), Fluid Mechanics and
#      Thermodynamics of Turbomachinery, 7th ed., Ch.3-4
#
#    Cumpsty (2004), Compressor Aerodynamics, Ch.5
#      IGV blade selection and solidity guidelines
#
#    NACA TN 1368 (Lieblein et al., 1956)
#      NACA 65-series profile data and camber-lift relations
#
#  Design philosophy:
#    The IGV is a set of *fixed* stator vanes upstream of
#    the rotor.  At the aerodynamic design point the IGV
#    exit flow angle (= rotor inlet absolute flow angle)
#    equals the prescribed pre-swirl angle alpha1.
#    For the SJTU rig the design intent is alpha1 = 0°
#    (axial inlet), so the IGV acts purely as a flow
#    straightener / reference surface.  Non-zero alpha1
#    would shift the rotor incidence and is the knob used
#    in variable-IGV configurations.
# ==============================================================

import numpy as np
from src.constants import gamma, R, Cp, T0_in, P0_in


# ---------------------------------------------------------------
# 1.  IGV aerodynamic sizing
# ---------------------------------------------------------------

def igv_geometry(
    D_tip,          # [m]   rotor tip diameter (= IGV tip diameter)
    nu,             # [-]   hub-to-tip ratio (same annulus as rotor)
    N_RPM,          # [rpm]
    phi,            # [-]   flow coefficient (= Ca / U_mean)
    alpha1_deg,     # [deg] desired IGV exit / rotor inlet swirl angle
                    #       positive = co-rotating with rotor
    V,              # [-]   number of IGV blades
    chord_igv=None, # [m]   IGV chord (mid-span); if None → estimated
    t_c=0.10,       # [-]   blade max thickness / chord (NACA 65)
    axial_gap_factor=1.10,   # IGV TE → rotor LE gap / rotor chord
    rotor_chord_mid=0.060,   # [m] rotor mid-span chord (60 mm, SJTU)
    eta_igv=0.995,           # total-pressure loss factor across IGV
):
    """
    Compute IGV geometry and velocity triangles at hub, mean, and tip.

    Returns a dict with:
      - annulus dimensions
      - velocity triangles at 3 radii
      - NACA 65-series blade parameters (camber, stagger, solidity)
      - axial positioning
      - loss estimate
    """
    # --- Annulus ---
    r_tip  = D_tip / 2.0
    r_hub  = nu * r_tip
    r_mean = 0.5 * (r_tip + r_hub)
    h_an   = r_tip - r_hub          # blade height
    A_ann  = np.pi * (r_tip**2 - r_hub**2)

    # --- Rotational speed ---
    omega  = 2.0 * np.pi * N_RPM / 60.0
    U_mean = omega * r_mean

    # --- Axial velocity from flow coefficient ---
    # phi = Ca / U_mean  →  Ca constant across span (free-vortex inlet)
    Ca = phi * U_mean

    # --- Mass flow (inlet total conditions, no swirl upstream of IGV) ---
    rho0   = P0_in / (R * T0_in)       # stagnation density
    # Static conditions at IGV inlet (pure axial flow, low Ma)
    Ma_ax  = Ca / np.sqrt(gamma * R * T0_in)   # approximate
    T_in   = T0_in / (1.0 + (gamma - 1) / 2.0 * Ma_ax**2)
    P_in   = P0_in * (T_in / T0_in) ** (gamma / (gamma - 1))
    rho_in = P_in / (R * T_in)
    mdot   = rho_in * Ca * A_ann

    # --- IGV exit velocity triangle at mean radius ---
    alpha1 = np.radians(alpha1_deg)
    # IGV turns the flow from axial to alpha1
    # Ca is constant (incompressible assumption, low Ma)
    # C1 = Ca / cos(alpha1),  C_theta1 = Ca * tan(alpha1)
    C1         = Ca / np.cos(alpha1)
    C_theta1   = Ca * np.tan(alpha1)
    # Rotor inlet relative velocity (at mean)
    W_theta1   = C_theta1 - U_mean      # < 0 for alpha1 < beta1
    W1         = np.sqrt(Ca**2 + W_theta1**2)
    beta1_mean = np.degrees(np.arctan(np.abs(W_theta1) / Ca))
    beta1_mean = -beta1_mean if W_theta1 < 0 else beta1_mean

    # --- IGV inlet (upstream): pure axial ---
    C0     = Ca                          # no swirl upstream
    alpha0 = 0.0                         # [deg]

    # --- IGV flow turning ---
    delta_alpha = alpha1_deg - alpha0    # flow deflection through IGV [deg]
    # Positive delta_alpha = co-swirl (in rotor direction)

    # ---------------------------------------------------------------
    # 2.  Radial distribution (free-vortex: r * C_theta = const)
    #     The IGV exit satisfies free-vortex condition so that the
    #     rotor sees uniform total enthalpy at every radius.
    # ---------------------------------------------------------------
    radii  = np.array([r_hub, r_mean, r_tip])
    labels = ['hub', 'mean', 'tip']
    stations = {}
    for r, lbl in zip(radii, labels):
        Ct1_r  = C_theta1 * r_mean / r   # free-vortex
        alpha1_r = np.degrees(np.arctan(Ct1_r / Ca))
        C1_r   = np.sqrt(Ca**2 + Ct1_r**2)
        U_r    = omega * r
        Wt1_r  = Ct1_r - U_r
        W1_r   = np.sqrt(Ca**2 + Wt1_r**2)
        beta1_r = np.degrees(np.arctan(np.abs(Wt1_r) / Ca))
        beta1_r = -beta1_r if Wt1_r < 0 else beta1_r
        stations[lbl] = {
            'r_m':       r * 1000,          # [mm]
            'Ca_m_s':    Ca,
            'U_m_s':     U_r,
            'alpha1_deg': alpha1_r,
            'C1_m_s':    C1_r,
            'C_theta1':  Ct1_r,
            'W1_m_s':    W1_r,
            'beta1_deg': beta1_r,
        }

    # ---------------------------------------------------------------
    # 3.  NACA 65-series blade design at mean radius
    #
    #     For a stator (IGV):
    #       inlet angle  = alpha0  = 0°   (axial)
    #       exit  angle  = alpha1  (desired swirl)
    #       flow deflection  = delta_alpha
    #
    #     NACA 65 camber-to-lift relation (thin-aerofoil theory,
    #     Lieblein 1956):
    #       CL ≈ 2 * pi * sin(theta_camber)   (isolated aerofoil)
    #       For cascade: CL_cascade ≈ 2 * (s/c) * [sin(alpha_m)]
    #       where alpha_m = mean flow angle = (alpha0 + alpha1) / 2
    #
    #     Lieblein's empirical relation for design incidence:
    #       i_des ≈ -6° + 0.18 * delta_alpha   (for NACA 65)
    #     Deviation rule (Carter, modified):
    #       delta_dev ≈ m * sqrt(camber / sigma)
    #       where m ≈ 0.23 for NACA 65 at mid-chord
    #
    #     We iterate: choose solidity sigma, compute camber from
    #     deflection + deviation, check Lieblein DF < 0.45
    # ---------------------------------------------------------------

    alpha_m     = (alpha0 + alpha1_deg) / 2.0   # mean flow angle [deg]
    # Solidity from Zweifel criterion: Zw ≈ 0.8 for lightly loaded stator
    # Zw = 2 * (s/c) * cos^2(alpha_m) * |tan(alpha0) - tan(alpha1)|
    # → sigma = c/s = 2 * cos^2(alpha_m) * |tan(alpha0)-tan(alpha1)| / Zw
    Zw          = 0.80
    tan_diff    = abs(np.tan(np.radians(alpha0)) - np.tan(np.radians(alpha1_deg)))
    if tan_diff < 1e-6:
        # axial IGV (alpha1 = 0): no aerodynamic loading; use structural min
        sigma_igv = 1.0
    else:
        sigma_igv = 2.0 * np.cos(np.radians(alpha_m))**2 * tan_diff / Zw

    # Chord from solidity: c = sigma * pitch = sigma * 2*pi*r_mean / V
    pitch_igv   = 2.0 * np.pi * r_mean / V
    if chord_igv is None:
        chord_igv = sigma_igv * pitch_igv

    # Recompute actual sigma with chosen chord
    sigma_act   = chord_igv / pitch_igv

    # Carter deviation:  delta_dev = 0.23 * sqrt(camber_deg / sigma_act)
    # Camber = deflection + deviation → iterate once
    camber_deg  = delta_alpha       # first guess (no deviation)
    for _ in range(10):
        dev         = 0.23 * np.sqrt(abs(camber_deg) / sigma_act) * np.sign(camber_deg)
        camber_deg  = delta_alpha + dev   # correction
    deviation_deg = dev

    # Design incidence (Lieblein): i_des ≈ -6 + 0.18 * |delta_alpha|
    # For axial IGV (delta_alpha=0) → i_des = 0 (stagger = inlet angle)
    i_des = -6.0 + 0.18 * abs(delta_alpha)
    if abs(delta_alpha) < 1.0:
        i_des = 0.0

    # Stagger angle (mid-chord line angle to axial)
    # kappa_LE = alpha0 - i_des  (blade leading-edge metal angle)
    # kappa_TE = alpha1 - deviation
    # stagger   = (kappa_LE + kappa_TE) / 2
    kappa_LE    = alpha0 - i_des
    kappa_TE    = alpha1_deg - deviation_deg
    stagger_igv = (kappa_LE + kappa_TE) / 2.0

    # Lieblein Diffusion Factor for IGV (stator form):
    # DF_stator = 1 - (C1/C0) + |dCt| / (2*sigma*C0)
    dCt         = abs(C_theta1 - 0.0)   # change in swirl
    DF_igv      = 1.0 - (C1 / C0) + dCt / (2.0 * sigma_act * C0)

    # Total pressure loss (profile loss approximation)
    # Denton (1993): omega_p ≈ 0.008 * (t/c) * (C1/C0)^3  per blade row
    # Use simplified: zeta_igv = 1 - eta_igv
    dP0_igv     = P0_in * (1.0 - eta_igv)   # [Pa] loss

    # ---------------------------------------------------------------
    # 4.  Axial positioning
    # ---------------------------------------------------------------
    axial_gap   = axial_gap_factor * rotor_chord_mid   # [m]  IGV-TE to Rotor-LE
    # IGV chord projected axially = chord * cos(stagger)
    igv_axial_length = chord_igv * np.cos(np.radians(stagger_igv))
    # Suggested IGV leading-edge axial position (upstream of rotor LE)
    igv_LE_to_rotor_LE = igv_axial_length + axial_gap  # [m]

    # ---------------------------------------------------------------
    # 5.  Tip clearance note (matching Tian 2018)
    # ---------------------------------------------------------------
    tip_clearance = 0.022 * h_an   # 2.2% of blade height [m]

    # ---------------------------------------------------------------
    # 6.  IGV-rotor blade count interaction (Tyler-Sofrin modes)
    # ---------------------------------------------------------------
    # Interaction mode orders: m = n*B + k*V  (Tyler & Sofrin 1962)
    # We flag the lowest-order interaction that could alias
    # given the 36-sensor ring (max unambiguous mode = 18 with 36 sensors)
    B = 21   # rotor blades (SJTU rig)
    interaction_modes = []
    for n in range(1, 4):
        for k in range(-3, 4):
            m = n * B + k * V
            interaction_modes.append((n, k, m))
    interaction_modes.sort(key=lambda x: abs(x[2]))

    return {
        # Annulus
        'D_tip_mm':          D_tip * 1000,
        'r_tip_mm':          r_tip * 1000,
        'r_hub_mm':          r_hub * 1000,
        'r_mean_mm':         r_mean * 1000,
        'h_annulus_mm':      h_an * 1000,
        'A_annulus_m2':      A_ann,
        # Flow
        'N_RPM':             N_RPM,
        'phi':               phi,
        'Ca_m_s':            Ca,
        'U_mean_m_s':        U_mean,
        'mdot_kg_s':         mdot,
        'Ma_axial':          Ma_ax,
        # IGV exit / rotor inlet (mean)
        'alpha1_des_deg':    alpha1_deg,
        'alpha0_deg':        alpha0,
        'delta_alpha_deg':   delta_alpha,
        'C0_m_s':            C0,
        'C1_mean_m_s':       C1,
        'C_theta1_mean':     C_theta1,
        'W1_mean_m_s':       W1,
        'beta1_mean_deg':    beta1_mean,
        # Radial stations
        'stations':          stations,
        # Blade geometry (mid-span)
        'V_blades':          V,
        'chord_igv_mm':      chord_igv * 1000,
        'pitch_igv_mm':      pitch_igv * 1000,
        'sigma_igv':         sigma_act,
        'stagger_igv_deg':   stagger_igv,
        'camber_igv_deg':    camber_deg,
        'deviation_deg':     deviation_deg,
        'i_des_deg':         i_des,
        'kappa_LE_deg':      kappa_LE,
        'kappa_TE_deg':      kappa_TE,
        'DF_igv':            DF_igv,
        't_c':               t_c,
        'max_thickness_mm':  chord_igv * t_c * 1000,
        # Axial positioning
        'axial_gap_mm':      axial_gap * 1000,
        'igv_axial_len_mm':  igv_axial_length * 1000,
        'igv_LE_to_rotor_LE_mm': igv_LE_to_rotor_LE * 1000,
        # Losses
        'dP0_igv_Pa':        dP0_igv,
        'eta_igv':           eta_igv,
        # Tip clearance
        'tip_clearance_mm':  tip_clearance * 1000,
        # Tyler-Sofrin interaction modes (lowest |m| first)
        'interaction_modes': interaction_modes[:8],
    }


# ---------------------------------------------------------------
# 2.  Meanline update: rotor analysis WITH pre-swirl
# ---------------------------------------------------------------

def meanline_with_igv(igv_res, PR, eta_is, sigma_rotor=1.0):
    """
    Compute rotor meanline given IGV exit conditions.

    Parameters
    ----------
    igv_res   : dict  — output of igv_geometry()
    PR        : float — total-to-total pressure ratio across the stage
    eta_is    : float — isentropic efficiency
    sigma_rotor : float — rotor solidity for Lieblein DF

    Returns
    -------
    dict with updated velocity triangles and aerodynamic coefficients
    """
    Ca       = igv_res['Ca_m_s']
    U        = igv_res['U_mean_m_s']
    Ct1      = igv_res['C_theta1_mean']
    W1       = igv_res['W1_mean_m_s']
    alpha1   = igv_res['alpha1_des_deg']

    # Thermodynamics
    T0_ratio    = PR ** ((gamma - 1.0) / gamma)
    dT0_is      = T0_in * (T0_ratio - 1.0)
    dT0         = dT0_is / eta_is
    W_euler     = Cp * dT0

    # Euler work: W = U * (Ct2 - Ct1)
    Ct2      = W_euler / U + Ct1
    C2       = np.sqrt(Ca**2 + Ct2**2)
    Wt2      = Ct2 - U
    W2       = np.sqrt(Ca**2 + Wt2**2)

    beta2    = np.degrees(np.arctan(np.abs(Wt2) / Ca))
    beta2    = -beta2 if Wt2 < 0 else beta2
    alpha2   = np.degrees(np.arctan(Ct2 / Ca))
    beta1    = igv_res['beta1_mean_deg']

    # De Haller
    DH = W2 / W1

    # Lieblein DF (rotor)
    DF = 1.0 - DH + abs(Ct2 - Ct1) / (2.0 * sigma_rotor * W1)

    # Dimensionless
    psi  = W_euler / U**2
    phi  = igv_res['phi']

    # Mass flow and shaft power
    rho_est = P0_in / (R * T0_in)
    mdot    = igv_res['mdot_kg_s']
    P_shaft = mdot * W_euler / 1000.0   # [kW]

    return {
        # Inlet (post-IGV)
        'alpha1_deg':     alpha1,
        'beta1_deg':      beta1,
        'C1_m_s':         igv_res['C1_mean_m_s'],
        'W1_m_s':         W1,
        'Ct1_m_s':        Ct1,
        # Exit
        'alpha2_deg':     alpha2,
        'beta2_deg':      beta2,
        'C2_m_s':         C2,
        'W2_m_s':         W2,
        'Ct2_m_s':        Ct2,
        # Work
        'PR':             PR,
        'eta_is':         eta_is,
        'dT0_K':          dT0,
        'W_euler_kJ_kg':  W_euler / 1000.0,
        'P_shaft_kW':     P_shaft,
        # Coefficients
        'psi':            psi,
        'phi':            phi,
        'De_Haller':      DH,
        'DF_rotor':       DF,
        'delta_beta_deg': abs(beta1) - abs(beta2),
        # Checks
        'psi_ok':         psi < 0.5,
        'DH_ok':          DH > 0.72,
        'DF_ok':          DF < 0.45,
    }


# ---------------------------------------------------------------
# 3.  Pretty-print summary
# ---------------------------------------------------------------

def print_igv_summary(res):
    SEP = "─" * 58
    print(f"\n{'═'*58}")
    print(f"  IGV DESIGN SUMMARY")
    print(f"{'═'*58}")

    print(f"\n{'ANNULUS':}")
    print(SEP)
    print(f"  Tip diameter          : {res['D_tip_mm']:.1f} mm")
    print(f"  Hub diameter          : {res['r_hub_mm']*2:.1f} mm")
    print(f"  Mean radius           : {res['r_mean_mm']:.1f} mm")
    print(f"  Blade height          : {res['h_annulus_mm']:.1f} mm")
    print(f"  Annular area          : {res['A_annulus_m2']*1e4:.2f} cm²")
    print(f"  Tip clearance         : {res['tip_clearance_mm']:.2f} mm  (2.2 % span)")

    print(f"\nFLOW CONDITIONS")
    print(SEP)
    print(f"  Speed                 : {res['N_RPM']:.0f} RPM")
    print(f"  Flow coefficient φ    : {res['phi']:.4f}")
    print(f"  Axial velocity Ca     : {res['Ca_m_s']:.2f} m/s")
    print(f"  Mean blade speed U    : {res['U_mean_m_s']:.2f} m/s")
    print(f"  Axial Mach number     : {res['Ma_axial']:.4f}")
    print(f"  Mass flow ṁ           : {res['mdot_kg_s']:.3f} kg/s")

    print(f"\nIGV AERODYNAMICS (mid-span)")
    print(SEP)
    print(f"  IGV inlet angle α₀    :  {res['alpha0_deg']:.1f}°  (axial)")
    print(f"  IGV exit angle α₁     :  {res['alpha1_des_deg']:.1f}°")
    print(f"  Flow deflection Δα    :  {res['delta_alpha_deg']:.1f}°")
    print(f"  Inlet velocity C₀     : {res['C0_m_s']:.2f} m/s")
    print(f"  Exit velocity C₁      : {res['C1_mean_m_s']:.2f} m/s")
    print(f"  Exit swirl C_θ1       : {res['C_theta1_mean']:.2f} m/s")
    print(f"  Rotor inlet W₁        : {res['W1_mean_m_s']:.2f} m/s")
    print(f"  Rotor inlet β₁        : {res['beta1_mean_deg']:.2f}°")

    print(f"\nBLADE GEOMETRY (NACA 65, mid-span)")
    print(SEP)
    print(f"  Number of IGV blades V: {res['V_blades']}")
    print(f"  Chord (mid-span)      : {res['chord_igv_mm']:.1f} mm")
    print(f"  Pitch (mid-span)      : {res['pitch_igv_mm']:.1f} mm")
    print(f"  Solidity σ            : {res['sigma_igv']:.3f}")
    print(f"  Stagger angle ξ       : {res['stagger_igv_deg']:.2f}°")
    print(f"  Camber angle θ_c      : {res['camber_igv_deg']:.2f}°")
    print(f"  Deviation δ           : {res['deviation_deg']:.2f}°")
    print(f"  Design incidence i    : {res['i_des_deg']:.2f}°")
    print(f"  LE metal angle κ_LE   : {res['kappa_LE_deg']:.2f}°")
    print(f"  TE metal angle κ_TE   : {res['kappa_TE_deg']:.2f}°")
    print(f"  Max thickness t       : {res['max_thickness_mm']:.1f} mm  (t/c={res['t_c']:.2f})")
    print(f"  Lieblein DF           : {res['DF_igv']:.4f}  {'✓' if res['DF_igv']<0.45 else '⚠ EXCEEDS 0.45'}")

    print(f"\nAXIAL POSITIONING")
    print(SEP)
    print(f"  IGV-rotor axial gap   : {res['axial_gap_mm']:.1f} mm  (1.10 × rotor chord)")
    print(f"  IGV axial projection  : {res['igv_axial_len_mm']:.1f} mm")
    print(f"  IGV LE to rotor LE    : {res['igv_LE_to_rotor_LE_mm']:.1f} mm")

    print(f"\nTOTAL-PRESSURE LOSS")
    print(SEP)
    print(f"  IGV η_total           : {res['eta_igv']:.4f}")
    print(f"  ΔP₀ IGV               : {res['dP0_igv_Pa']:.1f} Pa")

    print(f"\nRADIAL DISTRIBUTION (free-vortex IGV exit)")
    print(SEP)
    print(f"  {'Station':<8} {'r [mm]':>8} {'α₁ [°]':>8} {'C₁ [m/s]':>10} {'W₁ [m/s]':>10} {'β₁ [°]':>8}")
    print(f"  {'':─<8} {'':─>8} {'':─>8} {'':─>10} {'':─>10} {'':─>8}")
    for lbl, s in res['stations'].items():
        print(f"  {lbl:<8} {s['r_m']:>8.1f} {s['alpha1_deg']:>8.2f} "
              f"{s['C1_m_s']:>10.2f} {s['W1_m_s']:>10.2f} {s['beta1_deg']:>8.2f}")

    print(f"\nTYLER-SOFRIN INTERACTION MODES  (B={21}, V={res['V_blades']})")
    print(SEP)
    print(f"  {'n':>4} {'k':>4} {'m = nB+kV':>12}  note")
    for (n, k, m) in res['interaction_modes']:
        note = "⚠ aliased on 36-sensor ring" if abs(m) > 18 else "✓ resolved"
        print(f"  {n:>4} {k:>4} {m:>12}  {note}")

    print(f"\n{'═'*58}\n")


def print_rotor_summary(r):
    SEP = "─" * 58
    print(f"\n{'═'*58}")
    print(f"  ROTOR (with IGV pre-swirl)  —  MEANLINE")
    print(f"{'═'*58}")
    print(f"  PR                    : {r['PR']:.3f}")
    print(f"  η_is                  : {r['eta_is']:.3f}")
    print(f"  ΔT₀                   : {r['dT0_K']:.2f} K")
    print(f"  Euler work            : {r['W_euler_kJ_kg']:.3f} kJ/kg")
    print(f"  Shaft power           : {r['P_shaft_kW']:.2f} kW")
    print(f"\n  {'':─<54}")
    print(f"  {'Quantity':<22} {'Inlet (1)':>12} {'Exit (2)':>12}")
    print(f"  {'':─<22} {'':─>12} {'':─>12}")
    print(f"  {'Abs angle α [°]':<22} {r['alpha1_deg']:>12.2f} {r['alpha2_deg']:>12.2f}")
    print(f"  {'Rel angle β [°]':<22} {r['beta1_deg']:>12.2f} {r['beta2_deg']:>12.2f}")
    print(f"  {'Abs vel C [m/s]':<22} {r['C1_m_s']:>12.2f} {r['C2_m_s']:>12.2f}")
    print(f"  {'Rel vel W [m/s]':<22} {r['W1_m_s']:>12.2f} {r['W2_m_s']:>12.2f}")
    print(f"  {'Swirl Cθ [m/s]':<22} {r['Ct1_m_s']:>12.2f} {r['Ct2_m_s']:>12.2f}")
    print(f"\n  Work coefficient ψ    : {r['psi']:.4f}  {'✓' if r['psi_ok'] else '⚠'}")
    print(f"  Flow coefficient φ    : {r['phi']:.4f}")
    print(f"  De Haller W₂/W₁      : {r['De_Haller']:.4f}  {'✓ >0.72' if r['DH_ok'] else '⚠ <0.72'}")
    print(f"  Lieblein DF (rotor)   : {r['DF_rotor']:.4f}  {'✓ <0.45' if r['DF_ok'] else '⚠ >0.45'}")
    print(f"  Blade turning Δβ      : {r['delta_beta_deg']:.2f}°")
    print(f"{'═'*58}\n")


# ---------------------------------------------------------------
# 4.  Entry point — SJTU rig parameters
# ---------------------------------------------------------------

if __name__ == "__main__":
    # SJTU rig (Tian et al. 2018) — design point, zero pre-swirl
    igv_res = igv_geometry(
        D_tip           = 0.598,    # 598 mm
        nu              = 0.70,
        N_RPM           = 3000,
        phi             = 0.151,    # design flow coefficient
        alpha1_deg      = 0.0,      # zero pre-swirl (α₁ = 0°)
        V               = 13,       # 13 IGV blades (Tian 2018)
        chord_igv       = None,     # auto-size from Zweifel
        t_c             = 0.10,     # NACA 65-010 (10% thickness)
        axial_gap_factor= 1.10,     # 110% rotor chord (Tian 2018)
        rotor_chord_mid = 0.060,    # 60 mm mid-span rotor chord
        eta_igv         = 0.995,
    )
    print_igv_summary(igv_res)

    # Rotor meanline with IGV conditioning
    # SJTU rig: LOW-SPEED research fan, not a high-PR stage.
    # From Tian 2018 Fig. 3:  ψ_TP ≈ 0.40 at design φ = 0.151
    # ψ_TP = ΔP0 / (½ρU²)  →  ΔP0 ≈ 1562 Pa  →  PR ≈ 1.0154
    # η_TS ≈ 0.78 (total-to-static, from Fig. 3)
    # We use η_TT ≈ 0.83 for total-to-total (typical offset ~5 pts)
    rotor_res = meanline_with_igv(
        igv_res     = igv_res,
        PR          = 1.0154,
        eta_is      = 0.83,
        sigma_rotor = 1.0,
    )
    print_rotor_summary(rotor_res)

    # --- Sensitivity: non-zero pre-swirl cases ---
    print("\nPRE-SWIRL SENSITIVITY  (RI sweep, Tian 2018 spirit)")
    print("─" * 60)
    print(f"  {'α₁ [°]':>8} {'φ_eff':>8} {'β₁ [°]':>8} {'W₁ [m/s]':>10} {'ψ':>8} {'DH':>8} {'DF':>8}")
    print(f"  {'':─>8} {'':─>8} {'':─>8} {'':─>10} {'':─>8} {'':─>8} {'':─>8}")
    for a1 in [-15, -10, -5, 0, 5, 10, 15]:
        g = igv_geometry(
            D_tip=0.598, nu=0.70, N_RPM=3000, phi=0.151,
            alpha1_deg=a1, V=13, t_c=0.10,
            axial_gap_factor=1.10, rotor_chord_mid=0.060,
        )
        rr = meanline_with_igv(g, PR=1.0154, eta_is=0.83)
        print(f"  {a1:>8.0f} {g['phi']:>8.4f} {rr['beta1_deg']:>8.2f} "
              f"{rr['W1_m_s']:>10.2f} {rr['psi']:>8.4f} "
              f"{rr['De_Haller']:>8.4f} {rr['DF_rotor']:>8.4f}")
    print()