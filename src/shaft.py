# ==============================================================
#  shaft.py  —  Shaft sizing, bearing selection, and Campbell diagram
#
#  Covers the mechanical design of the rotating assembly for a
#  single-stage axial compressor rig in the design space:
#    D_tip  700–1000 mm
#    N      3000–4000 RPM  (direct motor drive, no gearbox)
#    PR     1.10–1.20
#
#  References:
#    Shigley & Mischke (2014), Mechanical Engineering Design, 10th ed.
#      Shaft design for combined bending and torsion, Ch. 6-7.
#
#    SKF General Catalogue (2018).
#      Angular contact ball bearing selection, L10 life, equivalent
#      dynamic load for combined radial + axial loading.
#
#    ISO 21940-11 (2016) — Mechanical vibration, rotor balancing.
#      Balance quality grades G1.0, G2.5 for rotating machinery.
#
#    ISO 281 (2007) — Rolling bearings — dynamic load ratings and
#      rating life.
#
#    Rankine (1869) — static deflection method for critical speed.
#      f_cr = (1/2π) √(g/δ_st)   [Hz]
#
#    Dunkerley (1894) — reciprocal sum for multi-mass systems.
#
#  Layout convention
#  -----------------
#  Between-bearings layout (rotor between the two bearing planes):
#
#    Motor ── Coupling ── [Bearing A] ── Rotor ── [Bearing B]
#                          (fixed)                  (free)
#
#  Bearing A (drive end, fixed): angular contact ball, locating,
#    carries both radial and axial (thrust) loads.
#  Bearing B (non-drive end, free): deep groove ball, non-locating,
#    carries radial load only, allows thermal axial expansion.
#
#  Coordinate system
#  -----------------
#  x along shaft axis, positive toward drive end.
#  Bearing A at x=0, rotor at x=a, Bearing B at x=L_span.
# ==============================================================

import numpy as np
from src.constants import gamma, R, Cp, T0_in, P0_in


# ---------------------------------------------------------------
# helpers
# ---------------------------------------------------------------

def _round_up_std(d_m, step_mm=5):
    """Round shaft diameter up to next standard step (mm)."""
    return np.ceil(d_m * 1000 / step_mm) * step_mm / 1000


# ---------------------------------------------------------------
# 1.  Rotor mass estimate
# ---------------------------------------------------------------

def rotor_mass(
    r_tip,           # [m]  tip radius
    r_hub,           # [m]  hub radius  (same as bellmouth hub)
    h_blade,         # [m]  blade height = r_tip - r_hub
    B,               # [-]  number of rotor blades
    chord_mid,       # [m]  rotor mid-span chord
    r_bore=None,     # [m]  shaft bore radius (None → 0.07 * r_tip)
    t_rim=0.050,     # [m]  disc rim axial width
    t_web=0.020,     # [m]  disc web (back-plate) axial thickness
    rho_disc=2810.0, # [kg/m³] disc material density (Al 7075 default)
    rho_blade=2810.0,# [kg/m³] blade material density (Al 7075 default)
    t_c=0.10,        # [-]  blade thickness ratio
):
    """
    Estimate rotor assembly mass using a thin-rim disc model.

    Disc geometry:
      - Rim   : annular ring, outer radius r_tip, width 50 mm axially,
                radial extent 50 mm (r_tip − 0.05 to r_tip)
      - Web   : thin back-plate from rim inner radius to bore radius,
                axial thickness t_web
    This is lighter and stiffer than a solid disc and gives better
    critical-speed margin.

    Returns dict with mass breakdown and polar moment of inertia.
    """
    if r_bore is None:
        r_bore = max(0.07 * r_tip, 0.025)  # min 25 mm bore

    r_rim_inner = r_tip - 0.050   # 50 mm radial rim width

    # Rim volume
    V_rim  = np.pi * (r_tip**2 - r_rim_inner**2) * t_rim
    m_rim  = rho_disc * V_rim

    # Web volume
    V_web  = np.pi * (r_rim_inner**2 - r_bore**2) * t_web
    m_web  = rho_disc * V_web

    m_disc = m_rim + m_web

    # Blades (rectangular planform approximation)
    t_avg   = 0.5 * t_c * chord_mid   # average thickness [m]
    V_blade = chord_mid * h_blade * t_avg
    m_blade = rho_blade * V_blade
    m_blades_total = B * m_blade

    m_total = m_disc + m_blades_total

    # Polar moment of inertia  J ≈ 0.5 m r²  (annular disc, approx)
    J_zz = 0.5 * m_disc * (r_tip**2 + r_bore**2) + \
           m_blades_total * (0.5*(r_tip + r_hub))**2

    return {
        'r_tip_m':         r_tip,
        'r_hub_m':         r_hub,
        'r_bore_m':        r_bore,
        'r_rim_inner_m':   r_rim_inner,
        't_rim_m':         t_rim,
        't_web_m':         t_web,
        'm_rim_kg':        m_rim,
        'm_web_kg':        m_web,
        'm_disc_kg':       m_disc,
        'm_blade_kg':      m_blade,
        'm_blades_kg':     m_blades_total,
        'm_total_kg':      m_total,
        'J_zz_kgm2':       J_zz,
        'B':               B,
        'rho_disc':        rho_disc,
    }


# ---------------------------------------------------------------
# 2.  Shaft sizing
# ---------------------------------------------------------------

def shaft_sizing(
    rotor_res,        # dict from rotor_mass()
    igv_res,          # dict from igv_geometry()
    meanline_res,     # dict from meanline_with_igv()
    L_span=0.700,     # [m]  bearing span (A to B)
    a_frac=0.50,      # [-]  rotor axial position as fraction of span
    E_shaft=200e9,    # [Pa] shaft Young's modulus (steel)
    sigma_allow=125e6,# [Pa] allowable stress (42CrMo4 QT, SF=2.0)
    G_grade=1.0,      # [mm/s] ISO balance grade (G1.0 = precision lab rig)
):
    """
    Size shaft diameter, compute static deflection, and estimate
    the first critical speed by Rankine's static deflection method.

    Parameters
    ----------
    L_span       : float  bearing span [m]
    a_frac       : float  rotor position from Bearing A / span (0.5 = centre)
    sigma_allow  : float  allowable combined stress [Pa]  (yield/SF)
    G_grade      : float  ISO 21940 balance grade [mm/s]

    Returns
    -------
    dict with shaft diameter, critical speed, bearing loads.
    """
    m_r   = rotor_res['m_total_kg']
    W_r   = m_r * 9.81                        # rotor weight [N]
    omega = igv_res['omega'] if 'omega' in igv_res else \
            2.0 * np.pi * igv_res['N_RPM'] / 60.0
    N_RPM = igv_res['N_RPM']
    P_W   = meanline_res['P_shaft_kW'] * 1000  # shaft power [W]
    T_t   = P_W / omega                         # torque [N·m]

    # Rotor position
    a     = a_frac * L_span
    b     = L_span - a

    # Static bearing reactions (simply supported)
    R_A   = W_r * b / L_span                   # Bearing A (fixed)
    R_B   = W_r * a / L_span                   # Bearing B (free)
    M_max = W_r * a * b / L_span               # max bending moment

    # Minimum shaft diameter — governed by GREATER of strength or stiffness
    # Strength: d >= [16/(pi*sig_allow) * sqrt(4M^2 + 3T^2)]^(1/3)
    d_strength = (16 / (np.pi * sigma_allow) *
                  np.sqrt(4 * M_max**2 + 3 * T_t**2)) ** (1.0/3.0)
    # Stiffness: enforce N_cr > 1.30 * N_op
    N_cr_target = N_RPM * 1.30
    delta_max   = 9.81 / (2.0 * np.pi * N_cr_target / 60.0) ** 2
    I_min_stiff = W_r * a**2 * b**2 / (3.0 * E_shaft * L_span * delta_max)
    d_stiffness = (64.0 * I_min_stiff / np.pi) ** 0.25
    d_adopt     = _round_up_std(max(d_strength, d_stiffness))

    # Second moment of area
    I_shaft = np.pi * d_adopt**4 / 64

    # Static deflection at rotor (Mohr's formula, simply supported)
    delta_st = (W_r * a**2 * b**2) / (3 * E_shaft * I_shaft * L_span)

    # First critical speed (Rankine)
    f_cr_Hz = (1 / (2 * np.pi)) * np.sqrt(9.81 / delta_st)
    N_cr    = f_cr_Hz * 60                     # RPM

    # Separation margin
    margin  = (N_cr - N_RPM) / N_cr            # positive → subcritical

    # Imbalance force (ISO 21940)
    e_perm   = (G_grade * 1e-3) / omega        # [m]  permissible eccentricity
    F_imb    = m_r * e_perm * omega**2         # [N]
    F_imb_design = m_r * G_grade * 1e-3 * omega  # equivalent: m*G*ω

    # Total radial load per bearing (worst case: imbalance adds to static)
    F_rad_A  = R_A + F_imb_design
    F_rad_B  = R_B + F_imb_design

    # Aerodynamic thrust (conservative: full static pressure rise × annulus)
    PR       = meanline_res['PR']
    dP_total = P0_in * (PR - 1.0)
    dP_stat  = 0.8 * dP_total          # approx static pressure rise
    A_ann    = igv_res['A_annulus_m2']
    F_thrust = dP_stat * A_ann

    return {
        # Shaft
        'd_strength_mm':  d_strength * 1000,
        'd_adopt_mm':     d_adopt * 1000,
        'I_shaft_m4':     I_shaft,
        'L_span_m':       L_span,
        'a_m':            a,
        'b_m':            b,
        'sigma_allow_Pa': sigma_allow,
        # Loads
        'W_rotor_N':      W_r,
        'T_torque_Nm':    T_t,
        'M_bending_Nm':   M_max,
        'R_A_N':          R_A,
        'R_B_N':          R_B,
        # Imbalance
        'G_grade':        G_grade,
        'e_perm_um':      e_perm * 1e6,
        'F_imbalance_N':  F_imb_design,
        'F_radial_A_N':   F_rad_A,
        'F_radial_B_N':   F_rad_B,
        'F_thrust_N':     F_thrust,
        # Critical speed
        'delta_st_um':    delta_st * 1e6,
        'N_cr_RPM':       N_cr,
        'f_cr_Hz':        f_cr_Hz,
        'N_op_RPM':       N_RPM,
        'subcritical':    N_cr > N_RPM,
        'margin_pct':     margin * 100,
        'ratio_N_Ncr':    N_RPM / N_cr,
        # Other
        'omega_rad_s':    omega,
    }


# ---------------------------------------------------------------
# 3.  Bearing selection
# ---------------------------------------------------------------

# Catalogue data: SKF angular contact ball bearings, 40° contact angle
# Columns: (bore [mm], OD [mm], B [mm], C [kN], C0 [kN])
_SKF_72xx = {
    # (bore mm, OD mm, B mm, C kN, C0 kN)  — SKF angular contact, 40° contact angle
    '7211 BEP': (55, 100, 21, 38.5, 30.5),
    '7311 BEP': (55, 120, 29, 72.0, 52.0),
    '7212 BEP': (60, 110, 22, 47.5, 38.0),
    '7312 BEP': (60, 130, 31, 83.2, 62.0),
    '7213 BEP': (65, 120, 23, 51.2, 55.0),
    '7313 BEP': (65, 140, 33, 93.0, 78.0),
    '7214 BEP': (70, 125, 24, 56.0, 60.0),
    '7314 BEP': (70, 150, 35, 104.0, 87.0),
    '7215 BEP': (75, 130, 25, 62.4, 67.0),
    '7315 BEP': (75, 160, 35, 118.0, 102.0),
    '7216 BEP': (80, 140, 26, 71.5, 80.0),
    '7316 BEP': (80, 170, 39, 132.0, 120.0),
    '7217 BEP': (85, 150, 28, 83.2, 93.5),
    '7317 BEP': (85, 180, 41, 153.0, 140.0),
    '7218 BEP': (90, 160, 30, 96.5, 112.0),
    '7318 BEP': (90, 190, 43, 170.0, 160.0),
    '7219 BEP': (95, 170, 32, 112.0, 132.0),
    '7319 BEP': (95, 200, 45, 193.0, 183.0),
    '7220 BEP': (100, 180, 34, 125.0, 150.0),
    '7320 BEP': (100, 215, 47, 216.0, 210.0),
}

_SKF_6xxx = {
    # (bore mm, OD mm, B mm, C kN, C0 kN)  — SKF deep groove ball bearings
    '6211':  (55, 100, 21, 35.1, 19.6),
    '6311':  (55, 120, 29, 55.9, 34.0),
    '6212':  (60, 110, 22, 41.0, 23.2),
    '6312':  (60, 130, 31, 64.4, 40.0),
    '6213':  (65, 120, 23, 44.0, 27.0),
    '6313':  (65, 140, 33, 72.8, 44.0),
    '6214':  (70, 125, 24, 48.8, 30.5),
    '6314':  (70, 150, 35, 86.4, 53.0),
    '6215':  (75, 130, 25, 52.7, 33.5),
    '6315':  (75, 160, 37, 96.0, 58.5),
    '6216':  (80, 140, 26, 62.4, 40.5),
    '6316':  (80, 170, 39, 108.0, 68.0),
    '6217':  (85, 150, 28, 72.8, 47.5),
    '6317':  (85, 180, 41, 122.0, 80.0),
    '6218':  (90, 160, 30, 80.6, 53.0),
    '6318':  (90, 190, 43, 137.0, 90.0),
    '6220':  (100, 180, 34, 98.0, 70.0),
    '6320':  (100, 215, 47, 174.0, 118.0),
}


def bearing_selection(
    shaft_res,          # dict from shaft_sizing()
    L10_target_h=20000, # [h] target L10 life
    p=3,                # bearing life exponent (3 = ball)
):
    """
    Select angular contact (locating) and deep groove (free) bearings
    based on equivalent dynamic load and L10 life.

    Uses SKF load-rating catalogue data embedded above.
    """
    N_RPM    = shaft_res['N_op_RPM']
    F_rad_A  = shaft_res['F_radial_A_N']   # locating bearing (A, fixed)
    F_rad_B  = shaft_res['F_radial_B_N']   # free bearing (B)
    F_thrust = shaft_res['F_thrust_N']
    d_mm     = shaft_res['d_adopt_mm']

    # L10 target in millions of revolutions
    L10_target_M = L10_target_h * N_RPM * 60 / 1e6

    # ---- Locating bearing (A): angular contact, carries thrust ----
    # Equivalent load  P = X Fr + Y Fa
    # SKF 7xx BEP (40° contact): e = 1.14 (load ratio threshold)
    # If Fa/(Fr) > e: X=0.35, Y=0.57; else X=1, Y=0
    # Use conservative single-row formula
    ratio_A = F_thrust / max(F_rad_A, 1.0)
    if ratio_A > 1.14:
        P_A = 0.35 * F_rad_A + 0.57 * F_thrust
    elif ratio_A > 0.68:
        P_A = 0.44 * F_rad_A + 0.72 * F_thrust
    else:
        P_A = F_rad_A

    C_req_A = P_A * L10_target_M ** (1.0 / p)

    # ---- Free bearing (B): deep groove, radial only ----
    P_B     = F_rad_B
    C_req_B = P_B * L10_target_M ** (1.0 / p)

    # Select smallest adequate bearing from catalogue
    def _select(catalogue, d_mm, C_req, P_load, N_RPM):
        candidates = []
        for name, (bore, OD, B_w, C_kN, C0_kN) in catalogue.items():
            if bore != d_mm:
                continue
            C  = C_kN * 1000
            if C < C_req:
                continue
            L10_M = (C / P_load) ** p if P_load > 0 else 1e9
            L10_h = L10_M * 1e6 / (N_RPM * 60)
            candidates.append({
                'name': name, 'bore': bore, 'OD': OD, 'B': B_w,
                'C_kN': C_kN, 'C0_kN': C0_kN,
                'P_N': P_load, 'L10_h': L10_h,
            })
        # Pick smallest OD among adequate
        if candidates:
            return sorted(candidates, key=lambda x: x['OD'])[0]
        # Fall back: pick the highest capacity in catalogue
        fb = sorted(catalogue.items(), key=lambda x: x[1][3], reverse=True)[0]
        name, (bore, OD, B_w, C_kN, C0_kN) = fb
        L10_M = (C_kN*1000/max(P_load,1))**p
        L10_h = L10_M*1e6/(N_RPM*60)
        return {'name': name, 'bore': bore, 'OD': OD, 'B': B_w,
                'C_kN': C_kN, 'C0_kN': C0_kN, 'P_N': P_load, 'L10_h': L10_h,
                'note': 'fallback - verify catalogue'}

    brg_A = _select(_SKF_72xx, d_mm, C_req_A, P_A, N_RPM)
    brg_B = _select(_SKF_6xxx, d_mm, C_req_B, P_B, N_RPM)

    return {
        'L10_target_h':  L10_target_h,
        'L10_target_M':  L10_target_M,
        # Locating (A)
        'P_A_N':         P_A,
        'C_req_A_kN':    C_req_A / 1000,
        'bearing_A':     brg_A,
        # Free (B)
        'P_B_N':         P_B,
        'C_req_B_kN':    C_req_B / 1000,
        'bearing_B':     brg_B,
        # Fa/Fr
        'Fa_Fr_ratio':   ratio_A,
    }


# ---------------------------------------------------------------
# 4.  Campbell diagram data
# ---------------------------------------------------------------

def campbell_data(shaft_res, igv_res, N_range=None):
    """
    Compute engine-order lines and critical speed crossings for
    a Campbell diagram.

    Returns data suitable for plotting or tabular display.
    """
    N_op  = shaft_res['N_op_RPM']
    N_cr  = shaft_res['N_cr_RPM']
    B     = igv_res['B_blades']
    V     = igv_res['V_blades']

    if N_range is None:
        N_range = np.linspace(0, N_op * 1.3, 400)

    # Engine orders of interest
    EO_list = [1, 2, 3, 4, 5, 6, 8, 10, 12, V, B]
    EO_list = sorted(set(EO_list))

    # Natural frequencies (first bending mode only for now)
    # f_n scales linearly with speed only for gyroscopic effects;
    # for a simple shaft model, f_n is constant (Rankine)
    f_n1_Hz = shaft_res['f_cr_Hz']

    crossings = []
    for EO in EO_list:
        # Resonance when EO × (N/60) = f_n1_Hz
        N_cross = f_n1_Hz * 60 / EO
        if 0 < N_cross < N_op * 1.3:
            note = ''
            if EO == B:
                note = 'rotor BPF'
            elif EO == V:
                note = 'IGV BPF'
            elif EO == 1:
                note = 'unbalance'
            elif EO == 2:
                note = '2× per rev'
            crossings.append({
                'EO':        EO,
                'N_cross':   N_cross,
                'f_cross_Hz': f_n1_Hz,
                'below_op':  N_cross < N_op,
                'note':      note,
            })

    return {
        'N_range':    N_range,
        'N_op_RPM':   N_op,
        'N_cr_RPM':   N_cr,
        'f_n1_Hz':    f_n1_Hz,
        'EO_list':    EO_list,
        'crossings':  crossings,
        'B':          B,
        'V':          V,
    }


# ---------------------------------------------------------------
# 5.  Summary printers
# ---------------------------------------------------------------

def print_shaft_summary(sr, rr):
    SEP = '─' * 60
    print(f"\n{'='*60}")
    print(f"  SHAFT & ROTOR MECHANICAL SUMMARY")
    print(f"{'='*60}")

    print(f"\nROTOR MASS  (thin-rim Al disc + Al blades)")
    print(SEP)
    print(f"  Disc rim          : {rr['m_rim_kg']:.1f} kg")
    print(f"  Disc web          : {rr['m_web_kg']:.1f} kg")
    print(f"  Blades (B={rr['B']})  : {rr['m_blades_kg']:.2f} kg")
    print(f"  Total rotor       : {rr['m_total_kg']:.1f} kg")
    print(f"  Rotor J_zz        : {rr['J_zz_kgm2']:.3f} kg·m²")
    print(f"  Bore radius       : {rr['r_bore_m']*1000:.0f} mm")

    print(f"\nSHAFT  (42CrMo4, QT)")
    print(SEP)
    print(f"  Required by strength : {sr['d_strength_mm']:.1f} mm")
    print(f"  Adopted diameter     : {sr['d_adopt_mm']:.0f} mm")
    print(f"  Allowable stress     : {sr['sigma_allow_Pa']/1e6:.0f} MPa  (SF=2.0)")
    print(f"  Bearing span         : {sr['L_span_m']*1000:.0f} mm")

    print(f"\nLOADS")
    print(SEP)
    print(f"  Rotor weight         : {sr['W_rotor_N']:.1f} N")
    print(f"  Shaft torque         : {sr['T_torque_Nm']:.1f} N·m")
    print(f"  Bending moment (max) : {sr['M_bending_Nm']:.1f} N·m")
    print(f"  Aerodynamic thrust   : {sr['F_thrust_N']:.0f} N  ({sr['F_thrust_N']/9.81:.0f} kgf)")
    print(f"  Balance grade        : G{sr['G_grade']:.1f}  (ISO 21940)")
    print(f"  Permissible e        : {sr['e_perm_um']:.2f} μm")
    print(f"  Residual imb. force  : {sr['F_imbalance_N']:.1f} N")

    print(f"\nBEARING LOADS")
    print(SEP)
    print(f"  Bearing A (fixed)    : Fr = {sr['F_radial_A_N']:.0f} N  +  Fa = {sr['F_thrust_N']:.0f} N")
    print(f"  Bearing B (free)     : Fr = {sr['F_radial_B_N']:.0f} N  (radial only)")

    print(f"\nCRITICAL SPEED  (Rankine, first bending mode)")
    print(SEP)
    print(f"  Static deflection    : {sr['delta_st_um']:.1f} μm")
    print(f"  First critical       : {sr['N_cr_RPM']:.0f} RPM  ({sr['f_cr_Hz']:.1f} Hz)")
    print(f"  Operating speed      : {sr['N_op_RPM']:.0f} RPM")
    print(f"  Ratio N_op / N_cr    : {sr['ratio_N_Ncr']:.3f}")
    status = f"SUBCRITICAL  ({sr['margin_pct']:.0f}% margin)" if sr['subcritical'] \
             else f"SUPERCRITICAL  (ratio = {sr['ratio_N_Ncr']:.2f})"
    print(f"  Status               : {status}")
    print(f"{'='*60}\n")


def print_bearing_summary(br):
    SEP = '─' * 60
    print(f"\n{'='*60}")
    print(f"  BEARING SELECTION SUMMARY")
    print(f"{'='*60}")
    print(f"  L10 life target      : {br['L10_target_h']:,} h")

    for label, side, key, brg_key in [
        ('BEARING A  (fixed / locating, angular contact)', 'A', 'A', 'bearing_A'),
        ('BEARING B  (free / non-locating, deep groove)',  'B', 'B', 'bearing_B'),
    ]:
        print(f"\n{label}")
        print(SEP)
        bg = br[brg_key]
        print(f"  Selected             : {bg['name']}")
        print(f"  Bore × OD × B        : {bg['bore']} × {bg['OD']} × {bg['B']} mm")
        print(f"  Dynamic rating C     : {bg['C_kN']:.1f} kN")
        print(f"  Static rating C0     : {bg['C0_kN']:.1f} kN")
        print(f"  Equivalent load P    : {bg['P_N']:.0f} N")
        print(f"  C required           : {br[f'C_req_{side}_kN']:.1f} kN")
        print(f"  Calculated L10       : {bg['L10_h']:,.0f} h  "
              f"{'✓' if bg['L10_h'] >= br['L10_target_h'] else '⚠'}")
    print(f"\n{'='*60}\n")


def print_campbell_summary(cd):
    SEP = '─' * 60
    print(f"\n{'='*60}")
    print(f"  CAMPBELL DIAGRAM  —  ENGINE ORDER CROSSINGS")
    print(f"{'='*60}")
    print(f"  Operating speed      : {cd['N_op_RPM']:.0f} RPM")
    print(f"  First critical (est) : {cd['N_cr_RPM']:.0f} RPM  ({cd['f_n1_Hz']:.1f} Hz)")
    print(f"  Rotor blades B       : {cd['B']}")
    print(f"  IGV blades V         : {cd['V']}")
    print(SEP)
    print(f"  {'EO':>4}  {'N_crossing [RPM]':>18}  {'Relation to N_op':>18}  Note")
    print(f"  {'':─>4}  {'':─>18}  {'':─>18}  {'':─>20}")
    for c in sorted(cd['crossings'], key=lambda x: x['N_cross'], reverse=True):
        rel = f"{'below op' if c['below_op'] else 'above op'}"
        warn = '⚠' if c['below_op'] and c['EO'] <= 6 else ' '
        note = c['note']
        print(f"  {c['EO']:>4}×  {c['N_cross']:>16.0f}  {rel:>18}  {note} {warn}")
    print(f"\n  Note: crossings above operating speed are passed through")
    print(f"  during run-up/down. Design for adequate damping at each.")
    print(f"{'='*60}\n")


# ---------------------------------------------------------------
# 6.  Entry point
# ---------------------------------------------------------------

if __name__ == "__main__":
    from src.igv import igv_geometry, meanline_with_igv

    igv_res = igv_geometry(
        D_tip=0.900, nu=0.75, N_RPM=3500, phi=0.725, alpha1_deg=0.0,
    )
    ml_res = meanline_with_igv(igv_res, PR=1.10, eta_is=0.85)

    # Rotor mass
    rr = rotor_mass(
        r_tip      = igv_res['r_tip_mm'] / 1000,
        r_hub      = igv_res['r_hub_mm'] / 1000,
        h_blade    = igv_res['h_annulus_mm'] / 1000,
        B          = igv_res['B_blades'],
        chord_mid  = igv_res['rotor_chord_mid_mm'] / 1000,
    )

    # Shaft sizing
    sr = shaft_sizing(rr, igv_res, ml_res, L_span=0.700, G_grade=1.0)

    # Bearing selection
    br = bearing_selection(sr, L10_target_h=20000)

    # Campbell diagram
    cd = campbell_data(sr, igv_res)

    # Print summaries
    print_shaft_summary(sr, rr)
    print_bearing_summary(br)
    print_campbell_summary(cd)