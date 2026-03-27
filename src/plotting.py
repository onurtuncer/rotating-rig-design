# ==============================================================
#  plotting.py — All matplotlib figures
# ==============================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

_COLOR_LIST  = ['#1D9E75', '#185FA5', '#D85A30', '#533AB7', '#854F0B']
_MARKER_LIST = ['o', 's', '^', 'D', 'v']


def _build_maps(df):
    """Build RPM→colour and PR→marker dicts from the actual DataFrame values."""
    rpms = sorted(df['N_RPM'].unique().tolist())
    prs  = sorted(df['PR'].unique().tolist())
    colors  = {n: _COLOR_LIST[i % len(_COLOR_LIST)]  for i, n in enumerate(rpms)}
    markers = {p: _MARKER_LIST[i % len(_MARKER_LIST)] for i, p in enumerate(prs)}
    return colors, markers


def plot_design_maps(df, eta_is, nu_target, h_min, save_path=None):
    """
    4-panel parametric design map.
      Panel 1 — Shaft power vs tip diameter
      Panel 2 — Smith chart (ψ vs φ)
      Panel 3 — Blade height vs tip diameter
      Panel 4 — De Haller number vs tip diameter
    """
    RPM_range = sorted(df['N_RPM'].unique().tolist())
    PR_range  = sorted(df['PR'].unique().tolist())
    COLORS_RPM, MARKERS_PR = _build_maps(df)
    mid_PR = PR_range[len(PR_range) // 2]   # solid line for the middle PR value

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        'Single-Stage Axial Compressor — Meanline Parametric Maps\n'
        f'η_is = {eta_is:.0%}  |  ν = {nu_target:.2f}',
        fontsize=14, fontweight='bold', y=1.01,
    )

    # Panel 1: Shaft power
    ax = axes[0, 0]
    for N in RPM_range:
        for PR in PR_range:
            sub = df[(df['N_RPM'] == N) & (df['PR'] == PR)].sort_values('D_tip_mm')
            ax.plot(sub['D_tip_mm'], sub['P_shaft_kW'],
                    color=COLORS_RPM[N], marker=MARKERS_PR[PR],
                    linestyle='-' if PR == mid_PR else '--',
                    alpha=0.85, linewidth=1.8,
                    label=f'{N} RPM' if PR == mid_PR else '')
    for kw, col, lbl in [(75, 'red', '75 kW'), (110, 'orange', '110 kW'), (150, 'purple', '150 kW')]:
        ax.axhline(y=kw, color=col, linestyle=':', linewidth=1.5, label=lbl)
    ax.set_xlabel('Tip Diameter [mm]')
    ax.set_ylabel('Shaft Power [kW]')
    ax.set_title('Shaft Power Requirement')
    ax.legend(fontsize=9)
    ax.set_xlim(650, 1050)

    # Panel 2: Smith chart
    ax = axes[0, 1]
    phi_bg = np.linspace(0.3, 0.75, 100)
    psi_hi = np.clip(0.55 - 0.5 * (phi_bg - 0.5)**2 + 0.2, 0, 0.55)
    ax.fill_between(phi_bg, 0.15, psi_hi, alpha=0.12, color='green', label='High-efficiency region')
    for _, row in df.iterrows():
        ax.scatter(row['phi'], row['psi'],
                   c=COLORS_RPM.get(row['N_RPM'], 'gray'),
                   marker=MARKERS_PR.get(row['PR'], 'o'),
                   s=60, alpha=0.8, edgecolors='white', linewidths=0.5)
    ax.legend(handles=[mpatches.Patch(color=c, label=f'{n} RPM') for n, c in COLORS_RPM.items()],
              fontsize=9, loc='upper right')
    ax.set_xlabel('Flow coefficient  φ = Ca/U')
    ax.set_ylabel('Work coefficient  ψ = W/U²')
    ax.set_title('Smith Chart — Design Points')
    ax.set_xlim(0.2, 0.9)
    ax.set_ylim(0, 0.7)

    # Panel 3: Blade height
    ax = axes[1, 0]
    for N in RPM_range:
        sub = df[df['N_RPM'] == N].groupby('D_tip_mm')['h_mm'].mean().reset_index()
        ax.plot(sub['D_tip_mm'], sub['h_mm'], color=COLORS_RPM[N], linewidth=2.2, label=f'{N} RPM')
    ax.axhline(y=h_min * 1000, color='red', linestyle='--', linewidth=2, label=f'h_min = {h_min*100:.0f} cm')
    ax.fill_between([650, 1050], 0, h_min * 1000, alpha=0.12, color='red', label='Forbidden zone')
    ax.set_xlabel('Tip Diameter [mm]')
    ax.set_ylabel('Blade Height h [mm]')
    ax.set_title('Blade Height Constraint')
    ax.legend(fontsize=9)
    ax.set_xlim(650, 1050)

    # Panel 4: De Haller
    ax = axes[1, 1]
    for N in RPM_range:
        for PR in PR_range:
            sub = df[(df['N_RPM'] == N) & (df['PR'] == PR)].sort_values('D_tip_mm')
            ax.plot(sub['D_tip_mm'], sub['De_Haller'],
                    color=COLORS_RPM[N], marker=MARKERS_PR[PR],
                    linestyle='-' if PR == mid_PR else '--',
                    alpha=0.85, linewidth=1.8)
    ax.axhline(y=0.72, color='red', linestyle='--', linewidth=2, label='De Haller = 0.72')
    ax.fill_between([650, 1050], 0, 0.72, alpha=0.1, color='red', label='Stall risk')
    ax.set_xlabel('Tip Diameter [mm]')
    ax.set_ylabel('De Haller Number  W₂/W₁')
    ax.set_title('Stall Criterion')
    ax.legend(fontsize=9)
    ax.set_xlim(650, 1050)
    ax.set_ylim(0.5, 1.1)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()


def plot_velocity_triangles(res, save_path=None):
    """
    Rotor inlet and exit velocity triangles at mean radius.
    Red = U (peripheral), Blue = C (absolute), Green = W (relative).
    """
    fig, ax = plt.subplots(figsize=(12, 5))

    U   = res['U_mean']
    Ca  = res['Ca']
    Ct2 = res['C_theta2']
    s   = 1.0 / U

    def draw_arrow(origin, vec, color):
        ax.annotate('', xy=origin + vec * s, xytext=origin,
                    arrowprops=dict(arrowstyle='->', color=color, lw=2.2))

    def draw_triangle(ox, oy, U_vec, C_vec, label):
        W_vec = C_vec - U_vec
        o = np.array([ox, oy])
        draw_arrow(o, U_vec, '#D85A30')
        ax.text(ox + U_vec[0]*s/2, oy - 0.09,
                f'U = {np.linalg.norm(U_vec):.0f} m/s',
                ha='center', color='#D85A30', fontsize=10)
        draw_arrow(o, C_vec, '#185FA5')
        ax.text(ox + C_vec[0]*s + 0.06, oy + C_vec[1]*s / 2,
                f'C = {np.linalg.norm(C_vec):.0f} m/s',
                color='#185FA5', fontsize=10)
        ax.annotate('', xy=o + U_vec*s + W_vec*s, xytext=o + U_vec*s,
                    arrowprops=dict(arrowstyle='->', color='#0F6E56', lw=2.2))
        ax.text(ox + U_vec[0]*s + W_vec[0]*s/2 - 0.08,
                oy + W_vec[1]*s/2 + 0.06,
                f'W = {np.linalg.norm(W_vec):.0f} m/s',
                color='#0F6E56', fontsize=10)
        ax.text(ox + 0.5, oy + C_vec[1]*s + 0.13,
                label, ha='center', fontsize=12, fontweight='bold')

    draw_triangle(0.1, 0, np.array([U, 0.0]), np.array([0.0, Ca]),  'Inlet  (1)')
    draw_triangle(1.6, 0, np.array([U, 0.0]), np.array([Ct2, Ca]), 'Exit  (2)')

    ax.set_title(
        f'Velocity Triangles — rₘ = {res["r_mean"]*1000:.0f} mm\n'
        f'β₁={res["beta1_deg"]:.1f}°  β₂={res["beta2_deg"]:.1f}°  '
        f'α₂={res["alpha2_deg"]:.1f}°  Δβ={res["delta_beta_deg"]:.1f}°', fontsize=12)
    ax.plot([], [], color='#D85A30', label='U — Peripheral speed')
    ax.plot([], [], color='#185FA5', label='C — Absolute velocity')
    ax.plot([], [], color='#0F6E56', label='W — Relative velocity')
    ax.legend(loc='lower right', fontsize=9)
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_xlim(-0.1, 2.9)
    ax.set_ylim(-0.3, 1.2)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()


def plot_radial(res, rad, save_path=None):
    """
    3-panel radial distribution plot from free_vortex() output.
    Panel 1 — Blade angles
    Panel 2 — Deflection and De Haller (dual x-axis)
    Panel 3 — Lieblein diffusion factor
    """
    pct  = rad['pct_span']
    fig, axes = plt.subplots(1, 3, figsize=(14, 7), sharey=True)
    fig.suptitle(
        f'Radial Distribution — Free Vortex\n'
        f'D_tip={res["D_tip_mm"]:.0f} mm  rₘ={res["r_mean"]*1000:.0f} mm  '
        f'N={res["N_RPM"]:.0f} RPM  PR={res["PR"]:.2f}',
        fontsize=13, fontweight='bold')

    ax = axes[0]
    ax.plot(rad['beta1'],  pct, color='#D85A30', lw=2.2, label='β₁ rotor inlet')
    ax.plot(rad['beta2'],  pct, color='#0F6E56', lw=2.2, label='β₂ rotor exit')
    ax.plot(rad['alpha2'], pct, color='#185FA5', lw=2.2, linestyle='--', label='α₂ stator inlet')
    ax.axhline(50, color='gray', lw=0.8, linestyle=':')
    ax.set_xlabel('Angle [°]')
    ax.set_ylabel('% Span')
    ax.set_title('Blade Angles')
    ax.legend(fontsize=9)

    ax  = axes[1]
    ax2 = ax.twiny()
    ax.plot(rad['delta_beta'], pct, color='#533AB7', lw=2.5, label='Δβ')
    ax2.plot(rad['DH'], pct, color='#185FA5', lw=2.0, linestyle='--', label='De Haller')
    ax2.axvline(0.72, color='red', lw=1.5, linestyle=':', label='DH=0.72')
    ax2.set_xlabel('De Haller W₂/W₁', color='#185FA5')
    ax2.tick_params(axis='x', colors='#185FA5')
    ax.set_xlabel('Deflection Δβ [°]')
    ax.set_title('Deflection & De Haller')
    ax.axhline(50, color='gray', lw=0.8, linestyle=':')
    lines1, labs1 = ax.get_legend_handles_labels()
    lines2, labs2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labs1 + labs2, fontsize=9)

    ax = axes[2]
    ax.plot(rad['DF'], pct, color='#854F0B', lw=2.5)
    ax.axvline(0.45, color='red', lw=1.5, linestyle='--', label='DF=0.45 limit')
    ax.fill_betweenx(pct, 0.45, rad['DF'].max() + 0.05, where=rad['DF'] > 0.45, alpha=0.12, color='red')
    ax.set_xlabel(f'Diffusion Factor DF  (σ={rad["sigma"]:.1f})')
    ax.set_title('Lieblein DF')
    ax.legend(fontsize=9)
    ax.axhline(50, color='gray', lw=0.8, linestyle=':')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()