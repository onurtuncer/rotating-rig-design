# ==============================================================
#  igv.py  —  Inlet Guide Vane aerodynamic design
#
#  General-purpose module for any axial compressor rig in the
#  design space:
#    D_tip  700–1000 mm
#    N      3000–4000 RPM
#    PR     1.10–1.20
#    nu     0.55–0.80
#
#  References:
#    Dixon & Hall (2014), Fluid Mechanics and Thermodynamics of
#      Turbomachinery, 7th ed., Ch. 3–4
#
#    Cumpsty (2004), Compressor Aerodynamics, Ch. 5
#      IGV blade selection, solidity, aspect ratio guidelines
#
#    Lieblein, Schwenk & Broderick (1953), NACA RM E53D01
#      NACA 65-series profile data, camber-lift, deviation rules
#
#    Carter (1950), ARC R&M 2804
#      Deviation rule: delta = m * sqrt(theta_c / sigma)
#
#    Tyler & Sofrin (1962), SAE Technical Paper 620532
#      Rotor-stator interaction mode orders m = nB + kV
#
#    Zweifel (1945) — loading criterion for turbomachinery stators
#
#  Design philosophy:
#    The IGV is a set of fixed stator vanes upstream of the rotor.
#    At the aerodynamic design point the IGV exit flow angle equals
#    the prescribed pre-swirl angle alpha1.  The default alpha1 = 0°
#    (axial inlet) makes the IGV a flow-straightening surface with
#    zero aerodynamic loading.  Non-zero alpha1 shifts the rotor
#    inlet incidence and is used for operating-point sweeps.
#
#    Blade count B and chord are auto-sized from annulus geometry
#    when not supplied by the caller.  The auto-sizing targets:
#      Rotor aspect ratio  AR ~ 1.5  (conservative, good stiffness)
#      Rotor solidity      σ  ~ 1.1  (mid-span, single-stage)
#      IGV/rotor count     V/B ~ 0.60, gcd(B,V) = 1  (no aliasing)
#      IGV solidity        σ_igv ~ 0.85 (lightly loaded stator)
# ==============================================================

import math
import numpy as np
from src.constants import gamma, R, Cp, T0_in, P0_in


# ---------------------------------------------------------------
# helpers
# ---------------------------------------------------------------

def _auto_blade_count(r_mean, h, sigma_target=1.1, AR=1.5):
    """
    Estimate rotor blade count B so that mid-span solidity ≈ sigma_target
    given aspect ratio AR = chord / h.

    sigma = c / pitch = (h/AR) / (2*pi*r_mean/B)
          → B = 2*pi*r_mean * sigma * AR / h
    """
    B = int(math.ceil(2 * math.pi * r_mean * sigma_target * AR / h))
    if B % 2 == 0:
        B += 1          # prefer odd (avoids even-order resonances)
    return B


def _auto_igv_count(B):
    """
    Choose IGV count V such that:
      - V / B ≈ 0.60
      - gcd(B, V) = 1  (no shared factors → no locked modes)
      - V is odd
    """
    V = round(B * 0.62)
    if V % 2 == 0:
        V -= 1
    while math.gcd(B, V) > 1:
        V -= 2
        if V < 3:
            V = B - 2
            break
    return max(V, 5)


# ---------------------------------------------------------------
# 1.  IGV aerodynamic sizing
# ---------------------------------------------------------------

def igv_geometry(
    D_tip,                  # [m]   rotor (and IGV) tip diameter
    nu,                     # [-]   hub-to-tip radius ratio
    N_RPM,                  # [rpm] rotational speed
    phi,                    # [-]   flow coefficient  Ca / U_mean
    alpha1_deg,             # [deg] IGV exit / rotor inlet swirl angle
                            #       0° = axial inlet (zero pre-swirl)
                            #       positive = co-rotating with rotor
    B=None,                 # [-]   rotor blade count (None → auto-size)
    V=None,                 # [-]   IGV blade count   (None → auto-size)
    rotor_chord_mid=None,   # [m]   rotor mid-span chord (None → auto AR=1.5)
    chord_igv=None,         # [m]   IGV mid-span chord   (None → auto)
    t_c=0.10,               # [-]   blade t/c ratio (NACA 65-series)
    axial_gap_factor=1.10,  # [-]   IGV-TE to rotor-LE gap / rotor chord
    tip_clearance_frac=0.022, # [-] tip clearance / blade height
    eta_igv=0.995,          # [-]   IGV total-pressure recovery
):
    """
    Compute IGV geometry and velocity triangles at hub, mean, and tip.

    All caller-supplied lengths are in metres.
    Returned dict uses mm for lengths and degrees for angles.
    """
    # --- Annulus ---
    r_tip  = D_tip / 2.0
    r_hub  = nu * r_tip
    r_mean = 0.5 * (r_tip + r_hub)
    h_an   = r_tip - r_hub
    A_ann  = np.pi * (r_tip**2 - r_hub**2)

    # --- Blade counts (auto-size if not supplied) ---
    if B is None:
        B = _auto_blade_count(r_mean, h_an)
    if V is None:
        V = _auto_igv_count(B)

    # --- Rotor chord (auto-size from AR=1.5 if not supplied) ---
    AR_rotor = 1.5
    if rotor_chord_mid is None:
        rotor_chord_mid = h_an / AR_rotor

    # --- Flow ---
    omega  = 2.0 * np.pi * N_RPM / 60.0
    U_mean = omega * r_mean
    Ca     = phi * U_mean

    # Static conditions (low-Ma approximation)
    Ma_ax  = Ca / np.sqrt(gamma * R * T0_in)
    T_in   = T0_in / (1.0 + (gamma - 1) / 2.0 * Ma_ax**2)
    P_in   = P0_in * (T_in / T0_in) ** (gamma / (gamma - 1))
    rho_in = P_in / (R * T_in)
    mdot   = rho_in * Ca * A_ann

    # --- IGV exit velocity triangle (mean radius) ---
    alpha1   = np.radians(alpha1_deg)
    C1       = Ca / np.cos(alpha1)
    C_theta1 = Ca * np.tan(alpha1)
    W_theta1 = C_theta1 - U_mean
    W1       = np.sqrt(Ca**2 + W_theta1**2)
    beta1_mean = np.degrees(np.arctan(np.abs(W_theta1) / Ca))
    if W_theta1 < 0:
        beta1_mean = -beta1_mean

    C0     = Ca
    alpha0 = 0.0
    delta_alpha = alpha1_deg - alpha0

    # --- Radial distribution (free-vortex exit) ---
    radii  = np.array([r_hub, r_mean, r_tip])
    labels = ['hub', 'mean', 'tip']
    stations = {}
    for r, lbl in zip(radii, labels):
        Ct1_r    = C_theta1 * r_mean / r
        alpha1_r = np.degrees(np.arctan(Ct1_r / Ca)) if Ca > 0 else 0.0
        C1_r     = np.sqrt(Ca**2 + Ct1_r**2)
        U_r      = omega * r
        Wt1_r    = Ct1_r - U_r
        W1_r     = np.sqrt(Ca**2 + Wt1_r**2)
        beta1_r  = np.degrees(np.arctan(np.abs(Wt1_r) / Ca))
        if Wt1_r < 0:
            beta1_r = -beta1_r
        stations[lbl] = {
            'r_mm':       r * 1000,
            'Ca_m_s':     Ca,
            'U_m_s':      U_r,
            'alpha1_deg': alpha1_r,
            'C1_m_s':     C1_r,
            'C_theta1':   Ct1_r,
            'W1_m_s':     W1_r,
            'beta1_deg':  beta1_r,
        }

    # --- NACA 65-series blade design at mean radius ---
    # Solidity from Zweifel criterion (Zw ≈ 0.80 for lightly loaded stator)
    alpha_m  = (alpha0 + alpha1_deg) / 2.0
    Zw       = 0.80
    tan_diff = abs(np.tan(np.radians(alpha0)) - np.tan(np.radians(alpha1_deg)))

    pitch_igv = 2.0 * np.pi * r_mean / V
    if chord_igv is None:
        if tan_diff < 1e-6:
            # Zero deflection: no aerodynamic constraint → use geometric minimum
            # sigma ≈ 1.0, chord = pitch
            chord_igv = pitch_igv
        else:
            sigma_zweifel = 2.0 * np.cos(np.radians(alpha_m))**2 * tan_diff / Zw
            chord_igv = sigma_zweifel * pitch_igv

    sigma_act = chord_igv / pitch_igv

    # Carter deviation (iterated)
    camber_deg = delta_alpha
    for _ in range(10):
        dev = 0.23 * np.sqrt(abs(camber_deg) / max(sigma_act, 0.01)) * np.sign(camber_deg)
        camber_deg = delta_alpha + dev
    deviation_deg = dev

    # Design incidence (Lieblein empirical)
    if abs(delta_alpha) < 1.0:
        i_des = 0.0
    else:
        i_des = -6.0 + 0.18 * abs(delta_alpha)

    kappa_LE  = alpha0 - i_des
    kappa_TE  = alpha1_deg - deviation_deg
    stagger_igv = (kappa_LE + kappa_TE) / 2.0

    # Lieblein DF for IGV (stator form)
    dCt   = abs(C_theta1 - 0.0)
    DF_igv = 1.0 - (C1 / C0) + dCt / (2.0 * sigma_act * C0)

    # Total-pressure loss
    dP0_igv = P0_in * (1.0 - eta_igv)

    # --- Axial positioning ---
    axial_gap        = axial_gap_factor * rotor_chord_mid
    igv_axial_length = chord_igv * np.cos(np.radians(stagger_igv))
    igv_LE_to_rotor_LE = igv_axial_length + axial_gap

    # --- Tip clearance ---
    tip_clearance = tip_clearance_frac * h_an

    # --- Tyler–Sofrin interaction modes ---
    ts_modes = []
    for n in range(1, 4):
        for k in range(-3, 4):
            m = n * B + k * V
            ts_modes.append((n, k, m))
    ts_modes.sort(key=lambda x: abs(x[2]))

    return {
        # Annulus
        'D_tip_mm':             D_tip * 1000,
        'r_tip_mm':             r_tip * 1000,
        'r_hub_mm':             r_hub * 1000,
        'r_mean_mm':            r_mean * 1000,
        'h_annulus_mm':         h_an * 1000,
        'A_annulus_m2':         A_ann,
        'nu':                   nu,
        # Flow
        'N_RPM':                N_RPM,
        'phi':                  phi,
        'Ca_m_s':               Ca,
        'U_mean_m_s':           U_mean,
        'mdot_kg_s':            mdot,
        'Ma_axial':             Ma_ax,
        # IGV exit / rotor inlet (mean)
        'alpha1_des_deg':       alpha1_deg,
        'alpha0_deg':           alpha0,
        'delta_alpha_deg':      delta_alpha,
        'C0_m_s':               C0,
        'C1_mean_m_s':          C1,
        'C_theta1_mean':        C_theta1,
        'W1_mean_m_s':          W1,
        'beta1_mean_deg':       beta1_mean,
        # Radial stations
        'stations':             stations,
        # Rotor blade sizing (used by bellmouth / downstream modules)
        'B_blades':             B,
        'rotor_chord_mid_mm':   rotor_chord_mid * 1000,
        # IGV blade geometry (mid-span)
        'V_blades':             V,
        'chord_igv_mm':         chord_igv * 1000,
        'pitch_igv_mm':         pitch_igv * 1000,
        'sigma_igv':            sigma_act,
        'stagger_igv_deg':      stagger_igv,
        'camber_igv_deg':       camber_deg,
        'deviation_deg':        deviation_deg,
        'i_des_deg':            i_des,
        'kappa_LE_deg':         kappa_LE,
        'kappa_TE_deg':         kappa_TE,
        'DF_igv':               DF_igv,
        't_c':                  t_c,
        'max_thickness_mm':     chord_igv * t_c * 1000,
        # Axial positioning
        'axial_gap_mm':         axial_gap * 1000,
        'igv_axial_len_mm':     igv_axial_length * 1000,
        'igv_LE_to_rotor_LE_mm': igv_LE_to_rotor_LE * 1000,
        # Losses
        'dP0_igv_Pa':           dP0_igv,
        'eta_igv':              eta_igv,
        # Tip clearance
        'tip_clearance_mm':     tip_clearance * 1000,
        'tip_clearance_frac':   tip_clearance_frac,
        # Tyler–Sofrin modes
        'ts_modes':             ts_modes[:8],
        # gcd check
        'gcd_BV':               math.gcd(B, V),
    }


# ---------------------------------------------------------------
# 2.  Rotor meanline WITH pre-swirl from IGV
# ---------------------------------------------------------------

def meanline_with_igv(igv_res, PR, eta_is, sigma_rotor=1.1):
    """
    Compute rotor meanline performance given IGV exit conditions.

    Parameters
    ----------
    igv_res     : dict   output of igv_geometry()
    PR          : float  stage total-to-total pressure ratio
    eta_is      : float  isentropic efficiency
    sigma_rotor : float  rotor mid-span solidity (default 1.1)

    Returns
    -------
    dict with velocity triangles, work/flow coefficients, and checks.
    """
    Ca    = igv_res['Ca_m_s']
    U     = igv_res['U_mean_m_s']
    Ct1   = igv_res['C_theta1_mean']
    W1    = igv_res['W1_mean_m_s']
    alpha1 = igv_res['alpha1_des_deg']

    dT0_is  = T0_in * (PR ** ((gamma - 1.0) / gamma) - 1.0)
    dT0     = dT0_is / eta_is
    W_euler = Cp * dT0

    Ct2  = W_euler / U + Ct1
    C2   = np.sqrt(Ca**2 + Ct2**2)
    Wt2  = Ct2 - U
    W2   = np.sqrt(Ca**2 + Wt2**2)

    beta2  = np.degrees(np.arctan(np.abs(Wt2) / Ca))
    if Wt2 > 0:
        beta2 = beta2        # exit swirl > blade speed (unusual)
    else:
        beta2 = -beta2
    alpha2 = np.degrees(np.arctan(Ct2 / Ca))
    beta1  = igv_res['beta1_mean_deg']

    DH = W2 / W1
    DF = 1.0 - DH + abs(Ct2 - Ct1) / (2.0 * sigma_rotor * W1)
    psi = W_euler / U**2
    phi = igv_res['phi']

    mdot   = igv_res['mdot_kg_s']
    P_kW   = mdot * W_euler / 1000.0

    # Rotor ΔP₀ (for loss budget reference in bellmouth)
    dP0_rotor = P0_in * (PR - 1.0)    # isentropic reference

    return {
        'alpha1_deg':     alpha1,
        'beta1_deg':      beta1,
        'C1_m_s':         igv_res['C1_mean_m_s'],
        'W1_m_s':         W1,
        'Ct1_m_s':        Ct1,
        'alpha2_deg':     alpha2,
        'beta2_deg':      beta2,
        'C2_m_s':         C2,
        'W2_m_s':         W2,
        'Ct2_m_s':        Ct2,
        'PR':             PR,
        'eta_is':         eta_is,
        'dT0_K':          dT0,
        'W_euler_kJ_kg':  W_euler / 1000.0,
        'P_shaft_kW':     P_kW,
        'dP0_rotor_Pa':   dP0_rotor,
        'psi':            psi,
        'phi':            phi,
        'De_Haller':      DH,
        'DF_rotor':       DF,
        'delta_beta_deg': abs(beta1) - abs(beta2),
        'psi_ok':         psi < 0.50,
        'DH_ok':          DH >= 0.72,
        'DF_ok':          DF <= 0.45,
    }


# ---------------------------------------------------------------
# 3.  Summary printers
# ---------------------------------------------------------------

def print_igv_summary(res):
    SEP = '─' * 60
    print(f"\n{'='*60}")
    print(f"  IGV DESIGN SUMMARY")
    print(f"{'='*60}")

    print(f"\nANNULUS")
    print(SEP)
    print(f"  Tip diameter          : {res['D_tip_mm']:.1f} mm")
    print(f"  Hub diameter          : {res['r_hub_mm']*2:.1f} mm")
    print(f"  Hub-to-tip ratio ν    : {res['nu']:.3f}")
    print(f"  Mean radius           : {res['r_mean_mm']:.1f} mm")
    print(f"  Blade height          : {res['h_annulus_mm']:.1f} mm")
    print(f"  Annular area          : {res['A_annulus_m2']*1e4:.2f} cm²")
    print(f"  Tip clearance         : {res['tip_clearance_mm']:.2f} mm"
          f"  ({res['tip_clearance_frac']*100:.1f}% span)")

    print(f"\nFLOW CONDITIONS")
    print(SEP)
    print(f"  Speed                 : {res['N_RPM']:.0f} RPM")
    print(f"  Flow coefficient φ    : {res['phi']:.4f}")
    print(f"  Axial velocity Ca     : {res['Ca_m_s']:.2f} m/s")
    print(f"  Mean blade speed U    : {res['U_mean_m_s']:.2f} m/s")
    print(f"  Axial Mach number     : {res['Ma_axial']:.4f}")
    print(f"  Mass flow ṁ           : {res['mdot_kg_s']:.3f} kg/s")

    print(f"\nIGV AERODYNAMICS  (mid-span)")
    print(SEP)
    print(f"  IGV inlet angle α₀    :  {res['alpha0_deg']:.1f}°  (axial)")
    print(f"  IGV exit angle α₁     :  {res['alpha1_des_deg']:.1f}°")
    print(f"  Flow deflection Δα    :  {res['delta_alpha_deg']:.1f}°")
    print(f"  Exit velocity C₁      : {res['C1_mean_m_s']:.2f} m/s")
    print(f"  Exit swirl C_θ1       : {res['C_theta1_mean']:.2f} m/s")
    print(f"  Rotor inlet W₁        : {res['W1_mean_m_s']:.2f} m/s")
    print(f"  Rotor inlet β₁        : {res['beta1_mean_deg']:.2f}°")

    print(f"\nBLADE GEOMETRY  (NACA 65, mid-span)")
    print(SEP)
    print(f"  Rotor blades B        : {res['B_blades']}")
    print(f"  Rotor chord (mid)     : {res['rotor_chord_mid_mm']:.1f} mm")
    print(f"  IGV blades V          : {res['V_blades']}")
    print(f"  gcd(B, V)             : {res['gcd_BV']}  "
          f"{'✓ no shared factors' if res['gcd_BV']==1 else '⚠ shared factors'}")
    print(f"  IGV chord (mid)       : {res['chord_igv_mm']:.1f} mm")
    print(f"  IGV pitch (mid)       : {res['pitch_igv_mm']:.1f} mm")
    print(f"  IGV solidity σ        : {res['sigma_igv']:.3f}")
    print(f"  Stagger angle ξ       : {res['stagger_igv_deg']:.2f}°")
    print(f"  Camber angle θ_c      : {res['camber_igv_deg']:.2f}°")
    print(f"  Deviation δ           : {res['deviation_deg']:.2f}°")
    print(f"  Design incidence i    : {res['i_des_deg']:.2f}°")
    print(f"  Lieblein DF (IGV)     : {res['DF_igv']:.4f}"
          f"  {'✓' if res['DF_igv']<0.45 else '⚠ exceeds 0.45'}")
    print(f"  Max thickness t       : {res['max_thickness_mm']:.1f} mm  (t/c={res['t_c']:.2f})")

    print(f"\nAXIAL POSITIONING")
    print(SEP)
    print(f"  IGV-rotor axial gap   : {res['axial_gap_mm']:.1f} mm"
          f"  (={res['axial_gap_mm']/res['rotor_chord_mid_mm']:.2f} × rotor chord)")
    print(f"  IGV axial projection  : {res['igv_axial_len_mm']:.1f} mm")
    print(f"  IGV LE → rotor LE     : {res['igv_LE_to_rotor_LE_mm']:.1f} mm")

    print(f"\nRADIAL DISTRIBUTION  (free-vortex IGV exit)")
    print(SEP)
    print(f"  {'Station':<8} {'r [mm]':>8} {'α₁ [°]':>8} {'C₁ [m/s]':>10} {'W₁ [m/s]':>10} {'β₁ [°]':>8}")
    print(f"  {'':─<8} {'':─>8} {'':─>8} {'':─>10} {'':─>10} {'':─>8}")
    for lbl, s in res['stations'].items():
        print(f"  {lbl:<8} {s['r_mm']:>8.1f} {s['alpha1_deg']:>8.2f}"
              f" {s['C1_m_s']:>10.2f} {s['W1_m_s']:>10.2f} {s['beta1_deg']:>8.2f}")

    print(f"\nTYLER–SOFRIN INTERACTION MODES  (B={res['B_blades']}, V={res['V_blades']})")
    print(SEP)
    print(f"  {'n':>4} {'k':>4} {'m = nB+kV':>12}")
    for (n, k, m) in res['ts_modes']:
        print(f"  {n:>4} {k:>4} {m:>12}")

    print(f"\n{'='*60}\n")


def print_rotor_summary(r, igv_res=None):
    SEP = '─' * 60
    print(f"\n{'='*60}")
    print(f"  ROTOR MEANLINE  (with IGV pre-swirl)")
    print(f"{'='*60}")
    print(f"  PR                    : {r['PR']:.4f}")
    print(f"  η_is                  : {r['eta_is']:.3f}")
    print(f"  ΔT₀                   : {r['dT0_K']:.2f} K")
    print(f"  Euler work            : {r['W_euler_kJ_kg']:.3f} kJ/kg")
    print(f"  Shaft power           : {r['P_shaft_kW']:.1f} kW")
    print(f"\n  {'':─<54}")
    print(f"  {'Quantity':<24} {'Inlet':>12} {'Exit':>12}")
    print(f"  {'':─<24} {'':─>12} {'':─>12}")
    print(f"  {'α abs [°]':<24} {r['alpha1_deg']:>12.2f} {r['alpha2_deg']:>12.2f}")
    print(f"  {'β rel [°]':<24} {r['beta1_deg']:>12.2f} {r['beta2_deg']:>12.2f}")
    print(f"  {'C abs [m/s]':<24} {r['C1_m_s']:>12.2f} {r['C2_m_s']:>12.2f}")
    print(f"  {'W rel [m/s]':<24} {r['W1_m_s']:>12.2f} {r['W2_m_s']:>12.2f}")
    print(f"  {'Cθ [m/s]':<24} {r['Ct1_m_s']:>12.2f} {r['Ct2_m_s']:>12.2f}")
    print(f"\n  Work coefficient ψ    : {r['psi']:.4f}  {'✓' if r['psi_ok'] else '⚠ >0.50'}")
    print(f"  Flow coefficient φ    : {r['phi']:.4f}")
    print(f"  De Haller W₂/W₁      : {r['De_Haller']:.4f}  {'✓ ≥0.72' if r['DH_ok'] else '⚠ <0.72'}")
    print(f"  Lieblein DF (rotor)   : {r['DF_rotor']:.4f}  {'✓ ≤0.45' if r['DF_ok'] else '⚠ >0.45'}")
    print(f"  Blade turning Δβ      : {r['delta_beta_deg']:.2f}°")
    print(f"{'='*60}\n")


# ---------------------------------------------------------------
# 4.  Entry point — new rig design space
# ---------------------------------------------------------------

if __name__ == "__main__":
    # Confirmed valid design point from meanline parametric sweep:
    #   D_tip = 900 mm,  N = 3500 RPM,  PR = 1.10,  ν = 0.75
    #   φ = 0.725,  η_is = 0.85
    #   DH = 0.736, DF = 0.447, P = 335 kW  (all within limits)
    igv_res = igv_geometry(
        D_tip      = 0.900,
        nu         = 0.75,
        N_RPM      = 3500,
        phi        = 0.725,
        alpha1_deg = 0.0,    # axial inlet (zero pre-swirl)
        # B, V, rotor_chord_mid → all auto-sized from annulus geometry
    )
    print_igv_summary(igv_res)

    rotor_res = meanline_with_igv(igv_res, PR=1.10, eta_is=0.85)
    print_rotor_summary(rotor_res)

    # Pre-swirl sensitivity
    print("PRE-SWIRL SENSITIVITY")
    print("─" * 62)
    print(f"  {'α₁ [°]':>7} {'β₁ [°]':>8} {'W₁ [m/s]':>10} {'ψ':>8} {'DH':>8} {'DF':>8} {'P [kW]':>8}")
    print(f"  {'':─>7} {'':─>8} {'':─>10} {'':─>8} {'':─>8} {'':─>8} {'':─>8}")
    for a1 in [-15, -10, -5, 0, 5, 10, 15]:
        g  = igv_geometry(D_tip=0.900, nu=0.75, N_RPM=3500, phi=0.725, alpha1_deg=a1)
        rr = meanline_with_igv(g, PR=1.10, eta_is=0.85)
        ok = '✓' if rr['DH_ok'] and rr['DF_ok'] else '⚠'
        print(f"  {a1:>7.0f} {rr['beta1_deg']:>8.2f} {rr['W1_m_s']:>10.2f}"
              f" {rr['psi']:>8.4f} {rr['De_Haller']:>8.4f}"
              f" {rr['DF_rotor']:>8.4f} {rr['P_shaft_kW']:>8.1f}  {ok}")
    print()