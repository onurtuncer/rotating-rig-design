# ==============================================================
#  bellmouth.py  —  Annular bellmouth inlet design
#                   outer bellmouth lip + inner centerbody
#
#  General-purpose module for any axial compressor rig in the
#  design space:
#    D_tip  700–1000 mm
#    N      3000–4000 RPM
#    PR     1.10–1.20
#
#  References:
#    ISO 5801 (2007/2017) — Industrial fans: performance testing
#      using standardised airways.  Bellmouth geometry, Cd,
#      measurement plane locations.
#
#    Mehta & Bradshaw (1979), Aeronautical Journal 83(821):737–749.
#      Design rules for low-turbulence settling sections:
#      honeycomb, screens, contraction ratio.
#
#    Bell & Mehta (1988), NASA CR-177488.
#      Contraction design for low-turbulence wind tunnels.
#
#    Hoerner (1965), Fluid Dynamic Drag.
#      Semi-ellipsoidal nose drag coefficient.
#
#  Physical picture
#  ----------------
#    OUTER : bellmouth lip  — ISO 5801 quarter-ellipse contraction
#            from D_lip down to D_tip (casing wall)
#
#    INNER : centerbody / hub spinner  — semi-ellipsoidal nose
#            from a point (nose tip) expanding to r_hub (IGV hub face)
#
#  The annulus between the two surfaces is the working flow path.
#  The centerbody shoulder aligns with the IGV leading-edge plane.
#
#  Coordinate system
#  -----------------
#    x = 0  at rotor leading edge  (positive = downstream)
#    All upstream stations have x < 0
# ==============================================================

import numpy as np
from src.constants import gamma, R, Cp, T0_in, P0_in


# ---------------------------------------------------------------
# 1.  Outer bellmouth contour  (ISO 5801 elliptical)
# ---------------------------------------------------------------

def outer_bellmouth_profile(r_casing, n_points=200, a_b_ratio=0.45):
    """
    Meridional contour (x, r) of the outer bellmouth wall.

    Quarter-ellipse from the lip apex (x=0) to the tangent point
    at x=a, then a straight cylindrical section of length 0.25*D_th.

    Parameters
    ----------
    r_casing  : float  casing (tip) radius [m]
    n_points  : int    contour resolution
    a_b_ratio : float  ISO 5801 axial/radial semi-axis ratio (~0.45)
    """
    D_th = 2.0 * r_casing
    b    = 0.10 * D_th          # radial semi-axis (lip bulge)
    a    = a_b_ratio * b        # axial  semi-axis

    theta = np.linspace(np.pi / 2, 0, n_points)
    x_arc = a * (1.0 - np.sin(theta))
    r_arc = r_casing + b * np.cos(theta)

    L_str = 0.25 * D_th         # ISO 5801 straight measurement section
    x_str = np.linspace(a, a + L_str, 20)
    r_str = np.full_like(x_str, r_casing)

    return {
        'r_casing_m':   r_casing,
        'r_lip_m':      r_casing + b,
        'D_lip_m':      2.0 * (r_casing + b),
        'b_m':          b,
        'a_m':          a,
        'a_b_ratio':    a_b_ratio,
        'L_arc_m':      a,
        'L_straight_m': L_str,
        'L_total_m':    a + L_str,
        'x_contour':    np.concatenate([x_arc, x_str[1:]]),
        'r_contour':    np.concatenate([r_arc, r_str[1:]]),
    }


# ---------------------------------------------------------------
# 2.  Inner centerbody  (semi-ellipsoidal spinner)
# ---------------------------------------------------------------

def centerbody_profile(r_hub, fineness=1.2, n_points=200):
    """
    Meridional contour (x, r) of the hub spinner.

    Profile:  r(x) = r_hub * sqrt(1 - (x/L_nose)²)
    Local x=0 at nose tip, x=L_nose at hub shoulder.

    Parameters
    ----------
    r_hub     : float  hub radius at IGV leading edge [m]
    fineness  : float  L_nose / r_hub  (1.2 = good drag / BL trade-off)
    """
    L_nose = fineness * r_hub
    x_cb   = np.linspace(0.0, L_nose, n_points)
    r_cb   = r_hub * np.sqrt(np.maximum(0.0, 1.0 - (x_cb / L_nose)**2))
    drr    = np.gradient(r_cb, x_cb)
    theta_shoulder = np.degrees(np.arctan(np.abs(drr[-2])))

    return {
        'r_hub_m':           r_hub,
        'L_nose_m':          L_nose,
        'fineness':          fineness,
        'x_contour':         x_cb,
        'r_contour':         r_cb,
        'theta_shoulder_deg': theta_shoulder,
    }


# ---------------------------------------------------------------
# 3.  Full annular bellmouth inlet design
# ---------------------------------------------------------------

def bellmouth_design(
    igv_res,
    contraction_ratio   = 4.0,
    Cd                  = 0.99,
    n_screens           = 1,
    screen_wire_d       = 0.5e-3,
    screen_mesh         = 16.0,
    honeycomb_cell      = 6.35e-3,
    honeycomb_L_D       = 8.0,
    centerbody_fineness = 1.2,
    T_ambient           = 288.15,
    P_ambient           = 101325.0,
):
    """
    Design the complete annular inlet assembly.

    Parameters
    ----------
    igv_res             : dict   output of igv_geometry()
    contraction_ratio   : float  settling-section area / net annular area
    Cd                  : float  bellmouth discharge coefficient
    n_screens           : int    turbulence reduction screens
    screen_wire_d       : float  [m] wire diameter
    screen_mesh         : float  mesh openings per inch
    honeycomb_cell      : float  [m] honeycomb cell diameter
    honeycomb_L_D       : float  honeycomb length / cell diameter
    centerbody_fineness : float  nose length / hub radius
    T_ambient, P_ambient: float  inlet stagnation conditions [K, Pa]

    Returns
    -------
    dict with all geometric and aerodynamic quantities.
    """
    r_tip  = igv_res['r_tip_mm']  / 1000.0
    r_hub  = igv_res['r_hub_mm']  / 1000.0
    h_an   = igv_res['h_annulus_mm'] / 1000.0
    A_ann  = igv_res['A_annulus_m2']
    Ca     = igv_res['Ca_m_s']
    Ma_ax  = igv_res['Ma_axial']

    # IGV leading-edge position (x < 0, upstream of rotor LE)
    x_igv_LE = -(igv_res['igv_LE_to_rotor_LE_mm'] / 1000.0)

    # --- Outer bellmouth ---
    outer        = outer_bellmouth_profile(r_tip)
    L_bell_total = outer['L_total_m']
    x_bell_exit  = x_igv_LE          # bellmouth exit aligns with IGV LE
    x_bell_entry = x_bell_exit - L_bell_total

    # --- Inner centerbody ---
    cb     = centerbody_profile(r_hub, fineness=centerbody_fineness)
    L_nose = cb['L_nose_m']
    x_cb_shoulder = x_igv_LE         # shoulder at IGV LE
    x_cb_nose_tip = x_cb_shoulder - L_nose

    # --- Settling section ---
    A_settle_outer = contraction_ratio * A_ann
    r_settle       = np.sqrt(A_settle_outer / np.pi)
    D_settle       = 2.0 * r_settle
    CR             = A_settle_outer / A_ann

    Ca_settle = Ca * A_ann / A_settle_outer
    Ma_settle = Ca_settle / np.sqrt(gamma * R * T_ambient)

    mesh_pitch  = 0.0254 / screen_mesh
    L_inter     = 3.0 * mesh_pitch
    L_honeycomb = honeycomb_L_D * honeycomb_cell
    L_settle_min = 0.5 * D_settle + L_honeycomb + n_screens * L_inter
    L_settle    = np.ceil(L_settle_min / 0.050) * 0.050

    x_settle_entry = x_bell_entry - L_settle
    nose_ok = x_cb_nose_tip > x_settle_entry

    # --- Throat conditions ---
    rho_th   = P_ambient / (R * T_ambient)
    T_throat = T_ambient - Ca**2 / (2.0 * Cp)
    P_throat = P_ambient * (T_throat / T_ambient) ** (gamma / (gamma - 1))
    q_th     = 0.5 * rho_th * Ca**2
    dP_stat  = P_ambient - P_throat
    mdot_act = Cd * rho_th * Ca * A_ann

    # --- Loss budget ---
    q_settle = 0.5 * rho_th * Ca_settle**2

    # Bellmouth surface friction
    dP0_bell = 0.005 * q_th

    # Centerbody nose drag (Hoerner 1965, semi-ellipsoid fineness 1.2)
    A_frontal = np.pi * r_hub**2
    dP0_cb    = 0.04 * q_settle * A_frontal / A_ann

    # Screens (Mehta & Bradshaw 1979)
    M_pitch   = mesh_pitch + screen_wire_d
    beta_scr  = (1.0 - screen_wire_d / M_pitch)**2
    K_scr     = (1.0 - beta_scr**2) / beta_scr**2
    dP0_screen = n_screens * K_scr * q_settle

    # Honeycomb
    K_honey   = 0.30 * np.sqrt(honeycomb_L_D)
    dP0_honey = K_honey * q_settle

    dP0_total  = dP0_bell + dP0_cb + dP0_screen + dP0_honey
    zeta_total = dP0_total / q_th

    # Fraction of rotor ΔP₀ — use value from rotor result if available,
    # otherwise fall back to isentropic estimate from igv_res metadata
    # (caller can pass rotor_res to get exact value)

    # --- Turbulence ---
    Tu_free  = 0.005
    Tu_honey = 0.10 * Tu_free
    Tu_scr   = Tu_honey
    for _ in range(n_screens):
        Tu_scr *= np.sqrt(K_scr / (K_scr + 1.0))
    Tu_throat = Tu_scr / np.sqrt(CR)

    # --- Boundary layers at IGV LE ---
    nu_air = 1.5e-5
    Re_hub = Ca * L_nose / nu_air
    delta_hub  = 0.37 * L_nose * Re_hub**(-0.2)
    dstar_hub  = delta_hub / 8.0

    Re_cas = Ca * outer['L_total_m'] / nu_air
    delta_cas  = 0.37 * outer['L_total_m'] * Re_cas**(-0.2)
    dstar_cas  = delta_cas / 8.0

    block_hub  = 2.0 * np.pi * r_hub * dstar_hub
    block_cas  = 2.0 * np.pi * r_tip * dstar_cas
    blockage   = (block_hub + block_cas) / A_ann

    # --- ISO 5801 measurement plane ---
    x_meas = x_bell_exit - outer['L_straight_m'] + 0.125 * 2.0 * r_tip

    return {
        # Annulus
        'r_tip_m':              r_tip,
        'r_hub_m':              r_hub,
        'h_annulus_mm':         h_an * 1000,
        'A_ann_m2':             A_ann,
        'nu':                   igv_res['nu'],
        # Outer bellmouth
        'outer':                outer,
        'D_lip_mm':             outer['D_lip_m'] * 1000,
        'L_bell_arc_mm':        outer['L_arc_m'] * 1000,
        'L_bell_straight_mm':   outer['L_straight_m'] * 1000,
        'L_bell_total_mm':      outer['L_total_m'] * 1000,
        'a_b_ratio':            outer['a_b_ratio'],
        # Centerbody
        'centerbody':           cb,
        'L_nose_mm':            L_nose * 1000,
        'fineness':             centerbody_fineness,
        'x_cb_nose_tip_mm':     x_cb_nose_tip * 1000,
        'x_cb_shoulder_mm':     x_cb_shoulder * 1000,
        'nose_inside_settling': nose_ok,
        # Settling section
        'contraction_ratio':    CR,
        'D_settle_mm':          D_settle * 1000,
        'Ca_settle_m_s':        Ca_settle,
        'Ma_settle':            Ma_settle,
        'L_settle_mm':          L_settle * 1000,
        'L_honeycomb_mm':       L_honeycomb * 1000,
        'honeycomb_cell_mm':    honeycomb_cell * 1000,
        'honeycomb_L_D':        honeycomb_L_D,
        'K_honey':              K_honey,
        'n_screens':            n_screens,
        'screen_wire_d_mm':     screen_wire_d * 1000,
        'screen_mesh':          screen_mesh,
        'beta_screen':          beta_scr,
        'K_screen':             K_scr,
        # Flow
        'Ca_m_s':               Ca,
        'Ma_ax':                Ma_ax,
        'mdot_kg_s':            mdot_act,
        'Cd':                   Cd,
        'T_throat_K':           T_throat,
        'P_throat_Pa':          P_throat,
        'q_throat_Pa':          q_th,
        'dP_static_Pa':         dP_stat,
        # Losses
        'dP0_bell_Pa':          dP0_bell,
        'dP0_centerbody_Pa':    dP0_cb,
        'dP0_screen_Pa':        dP0_screen,
        'dP0_honey_Pa':         dP0_honey,
        'dP0_total_Pa':         dP0_total,
        'zeta_total':           zeta_total,
        # Turbulence
        'Tu_freestream':        Tu_free,
        'Tu_throat':            Tu_throat,
        # Boundary layer
        'delta_hub_mm':         delta_hub * 1000,
        'dstar_hub_mm':         dstar_hub * 1000,
        'delta_cas_mm':         delta_cas * 1000,
        'dstar_cas_mm':         dstar_cas * 1000,
        'blockage_frac_pct':    blockage * 100,
        # Axial station map (mm, x=0 at rotor LE)
        'x_rotor_LE_mm':        0.0,
        'x_igv_TE_mm':          -igv_res['axial_gap_mm'],
        'x_igv_LE_mm':          x_igv_LE * 1000,
        'x_bell_exit_mm':       x_bell_exit * 1000,
        'x_bell_entry_mm':      x_bell_entry * 1000,
        'x_cb_nose_tip_mm':     x_cb_nose_tip * 1000,
        'x_settle_entry_mm':    x_settle_entry * 1000,
        'x_meas_mm':            x_meas * 1000,
        'L_inlet_total_mm':     (L_settle + L_bell_total) * 1000,
    }


# ---------------------------------------------------------------
# 4.  Pretty-print summary
# ---------------------------------------------------------------

def print_bellmouth_summary(b, rotor_dP0_Pa=None):
    """
    Print full bellmouth design summary.

    Parameters
    ----------
    b            : dict   output of bellmouth_design()
    rotor_dP0_Pa : float  rotor total pressure rise [Pa] for loss budget
                          percentage; if None, computed from isentropic PR
                          stored in the igv_res that was passed at design time.
                          Can be obtained from meanline_with_igv result as
                          rotor_res['dP0_rotor_Pa'].
    """
    SEP = '─' * 60

    def pct_q(v):
        return v / b['q_throat_Pa'] * 100

    def pct_rotor(v):
        if rotor_dP0_Pa and rotor_dP0_Pa > 0:
            return f"  ({v/rotor_dP0_Pa*100:.3f}% of rotor ΔP₀)"
        return ''

    print(f"\n{'='*60}")
    print(f"  ANNULAR BELLMOUTH INLET  —  DESIGN SUMMARY")
    print(f"{'='*60}")

    print(f"\nANNULUS AT THROAT")
    print(SEP)
    print(f"  Casing radius r_tip   : {b['r_tip_m']*1000:.1f} mm")
    print(f"  Hub radius r_hub      : {b['r_hub_m']*1000:.1f} mm")
    print(f"  Hub-to-tip ratio ν    : {b['nu']:.3f}")
    print(f"  Blade height          : {b['h_annulus_mm']:.1f} mm")
    print(f"  Net annular area      : {b['A_ann_m2']*1e4:.2f} cm²")
    print(f"  Axial velocity Ca     : {b['Ca_m_s']:.3f} m/s")
    print(f"  Throat Mach number    : {b['Ma_ax']:.5f}")
    print(f"  Dynamic pressure q    : {b['q_throat_Pa']:.2f} Pa")

    print(f"\nOUTER BELLMOUTH  (ISO 5801 Annex B, elliptical)")
    print(SEP)
    print(f"  Lip diameter D_lip    : {b['D_lip_mm']:.1f} mm")
    print(f"  Casing diameter D_th  : {b['r_tip_m']*2000:.1f} mm")
    print(f"  Radial overhang b     : {b['outer']['b_m']*1000:.1f} mm  (0.10 × D_th)")
    print(f"  Axial semi-axis a     : {b['L_bell_arc_mm']:.1f} mm")
    print(f"  a/b                   : {b['a_b_ratio']:.3f}  (ISO 5801 target ~0.45)")
    print(f"  Straight section      : {b['L_bell_straight_mm']:.1f} mm  (0.25 × D_th)")
    print(f"  Bellmouth total       : {b['L_bell_total_mm']:.1f} mm")

    print(f"\nINNER CENTERBODY  (semi-ellipsoidal spinner)")
    print(SEP)
    print(f"  Hub radius (at IGV)   : {b['r_hub_m']*1000:.1f} mm")
    print(f"  Max diameter          : {b['r_hub_m']*2000:.1f} mm")
    print(f"  Fineness L/r_hub      : {b['fineness']:.2f}")
    print(f"  Nose cone length      : {b['L_nose_mm']:.1f} mm")
    print(f"  Shoulder angle        : {b['centerbody']['theta_shoulder_deg']:.1f}°  "
          f"(tangent to hub wall — no kink)")
    print(f"  Nose tip position     : x = {b['x_cb_nose_tip_mm']:+.1f} mm")
    print(f"  Hub shoulder position : x = {b['x_cb_shoulder_mm']:+.1f} mm  (= IGV LE)")
    chk = "OK" if b['nose_inside_settling'] else "WARNING — check geometry"
    print(f"  Nose tip in settling  : {chk}")

    print(f"\nSETTLING SECTION")
    print(SEP)
    print(f"  Duct diameter         : {b['D_settle_mm']:.1f} mm")
    print(f"  Contraction ratio CR  : {b['contraction_ratio']:.2f}  (annulus area basis)")
    print(f"  Velocity in settling  : {b['Ca_settle_m_s']:.4f} m/s")
    print(f"  Ma in settling        : {b['Ma_settle']:.6f}")
    print(f"  Length                : {b['L_settle_mm']:.0f} mm")
    print(f"  Honeycomb cell ∅      : {b['honeycomb_cell_mm']:.2f} mm  (1/4\")")
    print(f"  Honeycomb length      : {b['L_honeycomb_mm']:.1f} mm  (L/D={b['honeycomb_L_D']:.0f})")
    scr = f"{b['n_screens']} × {b['screen_wire_d_mm']:.1f} mm wire, {b['screen_mesh']:.0f} mesh/in"
    print(f"  Screens               : {scr}  (β={b['beta_screen']:.3f}, K={b['K_screen']:.3f})")

    print(f"\nTOTAL PRESSURE LOSS BUDGET")
    print(SEP)
    print(f"  {'Source':<32} {'ΔP₀ [Pa]':>9}  {'% of q_th':>9}")
    print(f"  {'':─<32} {'':─>9}  {'':─>9}")
    rows = [
        ("Outer bellmouth friction",  b['dP0_bell_Pa']),
        ("Centerbody nose drag",      b['dP0_centerbody_Pa']),
        (f"Screens ×{b['n_screens']}", b['dP0_screen_Pa']),
        ("Honeycomb",                 b['dP0_honey_Pa']),
    ]
    for name, dp in rows:
        print(f"  {name:<32} {dp:>9.2f}  {pct_q(dp):>8.3f}%")
    print(f"  {'':─<32} {'':─>9}")
    print(f"  {'TOTAL':<32} {b['dP0_total_Pa']:>9.2f}  {b['zeta_total']*100:>8.3f}%"
          + pct_rotor(b['dP0_total_Pa']))
    print(f"  Discharge coeff Cd    : {b['Cd']:.4f}")

    print(f"\nBOUNDARY LAYER AT IGV LE  (flat-plate turbulent estimate)")
    print(SEP)
    print(f"  Hub  δ / δ*           : {b['delta_hub_mm']:.2f} mm / {b['dstar_hub_mm']:.3f} mm")
    print(f"  Casing δ / δ*         : {b['delta_cas_mm']:.2f} mm / {b['dstar_cas_mm']:.3f} mm")
    chk2 = "OK (<1%)" if b['blockage_frac_pct'] < 1.0 else "Review (>1%)"
    print(f"  Area blockage         : {b['blockage_frac_pct']:.3f}%  {chk2}")

    print(f"\nTURBULENCE INTENSITY")
    print(SEP)
    print(f"  Lab freestream Tu     : {b['Tu_freestream']*100:.2f}%")
    print(f"  At rotor face (est.)  : {b['Tu_throat']*100:.4f}%")

    print(f"\nAXIAL STATION MAP  (x = 0 at rotor LE,  upstream = negative)")
    print(SEP)
    stations = [
        ("Settling section entry",      b['x_settle_entry_mm']),
        ("Centerbody nose tip",         b['x_cb_nose_tip_mm']),
        ("Bellmouth lip (outer entry)", b['x_bell_entry_mm']),
        ("ISO 5801 meas. plane",        b['x_meas_mm']),
        ("Bellmouth exit / IGV LE",     b['x_bell_exit_mm']),
        ("IGV trailing edge",           b['x_igv_TE_mm']),
        ("Rotor leading edge",          b['x_rotor_LE_mm']),
    ]
    for name, x in sorted(stations, key=lambda s: s[1]):
        print(f"  {name:<34}: x = {x:+8.1f} mm")

    print(f"\nASSEMBLY TOTALS")
    print(SEP)
    print(f"  Total inlet length    : {b['L_inlet_total_mm']:.0f} mm")
    print(f"  Settling duct ∅       : {b['D_settle_mm']:.1f} mm")
    print(f"  Outer lip ∅           : {b['D_lip_mm']:.1f} mm")
    print(f"  Centerbody max ∅      : {b['r_hub_m']*2000:.1f} mm  (at IGV LE)")
    print(f"\n{'='*60}\n")


# ---------------------------------------------------------------
# 5.  Entry point
# ---------------------------------------------------------------

if __name__ == "__main__":
    from src.igv import igv_geometry, meanline_with_igv

    igv_res = igv_geometry(
        D_tip      = 0.900,
        nu         = 0.75,
        N_RPM      = 3500,
        phi        = 0.725,
        alpha1_deg = 0.0,
    )

    rotor_res = meanline_with_igv(igv_res, PR=1.10, eta_is=0.85)

    bell = bellmouth_design(igv_res, contraction_ratio=4.0, centerbody_fineness=1.2)
    print_bellmouth_summary(bell, rotor_dP0_Pa=rotor_res['dP0_rotor_Pa'])