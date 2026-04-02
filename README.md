# Compressor Rotating Rig Design

Preliminary design of a single-stage axial compressor rotating rig for turbomachinery research. The rig is designed for direct motor drive (no gearbox), a pressure ratio of 1.1–1.2, and a tip diameter in the 700–1000 mm range. An inlet guide vane (IGV) row is included upstream of the rotor and can be removed or made variable for different test programmes.

The workflow follows the classical meanline-first approach:

```
Meanline parametric sweep  →  pick mean radius  →  IGV design
→  bellmouth inlet  →  shaft & bearings  →  blade design  →  instrumentation
```

---

## Design objectives

| Parameter | Value |
|---|---|
| Configuration | Single-stage: IGV (optional) + rotor + stator |
| Tip diameter | 700–1000 mm |
| Pressure ratio | 1.10–1.20 |
| Blade height | ≥ 50 mm |
| Speed | 3 000–4 000 RPM, direct motor drive |
| Drive | 4-pole or 6-pole induction motor + VFD, no gearbox |

---

## Reference design point

The parametric meanline sweep over the full design space (D_tip, N, PR, ν, φ) identified the following point as a well-balanced starting configuration. All module outputs below are computed from this point.

| Parameter | Symbol | Value |
|---|---|---|
| Tip diameter | D_tip | 900 mm |
| Hub diameter | D_hub | 675 mm |
| Hub-to-tip ratio | ν | 0.75 |
| Blade height | h | 112.5 mm |
| Rotational speed | N | 3 500 RPM |
| Flow coefficient | φ | 0.725 |
| Pressure ratio | PR | 1.10 |
| Isentropic efficiency | η_is | 0.85 |
| Tip blade speed | U_tip | 164.9 m/s |
| Tip Mach number | M_tip | 0.485 |
| Axial velocity | Ca | 104.6 m/s |
| Mean blade speed | U_mean | 144.3 m/s |
| Mass flow | ṁ | 34.0 kg/s |
| Shaft power | P | 320 kW |
| Temperature rise | ΔT₀ | 9.36 K |
| Work coefficient | ψ | 0.451 |
| De Haller number | W₂/W₁ | 0.736 ✓ |
| Lieblein DF | DF | 0.430 ✓ |
| Rotor inlet angle | β₁ | −54.1° |
| Rotor exit angle | β₂ | −37.1° |
| Stator inlet angle | α₂ | 31.9° |

Smith-chart position (ψ ≈ 0.45, φ ≈ 0.73) sits in the high-efficiency island for single-stage axial compressors. The 320 kW power requirement matches a standard IEC 315 or IEC 355 frame 4-pole induction motor running at 3 500 RPM via VFD.

---

## Repository structure

```
rotating-rig-design/
├── README.md
├── requirements.txt
└── src/
    ├── __init__.py          package exports
    ├── constants.py         ISA sea-level conditions, air properties
    ├── meanline.py          rotor-only meanline (no IGV)
    ├── radial.py            free-vortex radial distribution, Lieblein DF
    ├── plotting.py          parametric design maps, velocity triangles
    ├── igv.py               IGV aerodynamic design  ★
    ├── bellmouth.py         annular bellmouth inlet design  ★
    └── shaft.py             shaft sizing, bearings, Campbell diagram  ★
```

Modules marked ★ are part of the current design programme.

---

## Module descriptions

### `constants.py`

ISA sea-level inlet total conditions: T₀ = 288.15 K, P₀ = 101 325 Pa, R = 287.05 J/(kg·K), γ = 1.4, Cₚ = 1 004.5 J/(kg·K).

### `meanline.py`

Rotor-only meanline analysis (zero pre-swirl inlet). Given tip diameter, speed, pressure ratio, and isentropic efficiency, returns velocity triangles at mean radius, work coefficient ψ, flow coefficient φ, De Haller number, Lieblein DF (at assumed solidity), shaft power, and mass flow. The primary tool for the parametric design-space sweep.

### `radial.py`

Free-vortex radial distribution from hub to tip. Computes β₁, β₂, α₂, blade turning Δβ, De Haller W₂/W₁, and Lieblein DF at each radial station. Flags stations where DF > 0.45.

### `plotting.py`

Four-panel parametric design map sweeping D_tip, N, and PR: shaft power, Smith chart (ψ vs φ), blade height, and De Haller number. Also produces velocity triangle diagrams and radial distribution plots.

### `igv.py` ★

Complete IGV aerodynamic design. Auto-sizes rotor blade count B and chord from annulus geometry targeting aspect ratio AR = 1.5 and mid-span solidity σ ≈ 1.1. Selects IGV count V with gcd(B, V) = 1 to avoid locked Tyler–Sofrin interaction modes. Applies the Zweifel loading criterion and Carter deviation rule to set IGV camber and stagger. Computes free-vortex velocity triangles at hub, mean, and tip, and evaluates Lieblein DF for the IGV row.

| Parameter | Value |
|---|---|
| Rotor blades B | 37 |
| Rotor chord (mid-span) | 75 mm |
| IGV blades V | 23 |
| gcd(B, V) | 1 (no shared factors) |
| IGV chord (mid-span) | 100 mm |
| IGV–rotor axial gap | 82 mm (1.10 × rotor chord) |
| Tip clearance | 2.47 mm (2.2 % of span) |

### `bellmouth.py` ★

Full annular inlet design with two matched stream surfaces.

**Outer bellmouth** — ISO 5801 Annex B quarter-ellipse contraction (a/b = 0.45) followed by a straight cylindrical measurement section (0.25 × D_th). Wall static taps in the straight section give the ISO 5801 mass-flow measurement.

**Inner centerbody** — semi-ellipsoidal hub spinner (fineness L/r_hub = 1.2, Hoerner Cd_nose ≈ 0.04). The shoulder aligns with the IGV leading-edge plane so the hub wall is continuous into the blade passage.

**Settling section** — sized to Mehta & Bradshaw (1979): contraction ratio 4, 1/4″ honeycomb (L/D = 8), one 16-mesh turbulence screen (0.5 mm wire).

| Parameter | Value |
|---|---|
| Bellmouth lip diameter | 1 080 mm |
| Settling duct diameter | 1 191 mm |
| Centerbody nose length | 405 mm |
| Total inlet assembly length | 966 mm |
| Total inlet ΔP₀ | 1 245 Pa (12.3 % of rotor ΔP₀) |
| Turbulence intensity at rotor face | 0.020 % |
| Annular BL blockage at IGV LE | 1.2 % |
| Discharge coefficient Cd | 0.99 |

### `shaft.py` ★

Shaft sizing, bearing selection, and Campbell diagram for the rotating assembly.

**Rotor mass** is estimated with a thin-rim disc model (50 mm annular rim + 20 mm web back-plate, Al 7075-T6), saving roughly 40 kg versus a solid disc and improving critical-speed margin.

**Shaft sizing** is governed by stiffness (N_cr > 1.30 × N_op), not strength. The strength calculation alone requires only 39.6 mm; keeping the first critical speed 30 % above the operating speed requires 55 mm.

**Bearing layout** is fixed–free: an angular contact pair at the drive end (locating, takes axial thrust), and a deep-groove ball bearing at the non-drive end (non-locating, radial only).

| Parameter | Value |
|---|---|
| Rotor total mass | 50.1 kg |
| Rotor polar moment J_zz | 5.3 kg·m² |
| Shaft material | 42CrMo4, quenched and tempered |
| Shaft diameter | 55 mm |
| Bearing span | 700 mm |
| First critical speed | 4 782 RPM (N_op / N_cr = 0.73) |
| Subcritical separation margin | 27 % |
| Balance grade | G1.0 (ISO 21940-11) |
| Permissible eccentricity | 2.73 μm |
| Locating bearing (A, drive end) | SKF 7211 BEP — 55 × 100 × 21 mm, C = 38.5 kN |
| Free bearing (B, non-drive end) | SKF 6211 — 55 × 100 × 21 mm, C = 35.1 kN |
| Calculated L10 (bearing A) | 103 750 h |
| Aerodynamic axial thrust | 2 256 N (230 kgf) |

**Campbell diagram** — all engine-order crossings with the first bending mode (4 782 RPM) fall below the operating speed and are encountered only during run-up and coast-down. The critical crossings to manage are:

| Engine order | Crossing speed | Note |
|---|---|---|
| 1× | 4 782 RPM | above operating — run-through only |
| 2× | 2 391 RPM | below operating — traverse on every start/stop |
| 3× | 1 594 RPM | below operating |

Programme ±200 RPM skip bands in the VFD for the 2× and 1× crossings. With G1.0 balancing the amplitudes are small and traversal is rapid.

---

## Axial station map

All positions referenced to the rotor leading edge (x = 0, positive downstream).

```
x = −1 156 mm   Settling section entry         (duct D = 1 191 mm)
x =   −595 mm   Centerbody nose tip
x =   −456 mm   Bellmouth lip (outer entry)    (D_lip = 1 080 mm)
x =   −303 mm   ISO 5801 measurement plane     (wall static taps)
x =   −190 mm   Bellmouth exit / IGV LE
x =    −82 mm   IGV trailing edge
x =      0 mm   Rotor leading edge             ◀ reference plane
```

Total inlet assembly (settling section entry to IGV LE): **966 mm**.

---

## Quick start

```python
from src.igv import igv_geometry, meanline_with_igv, print_igv_summary, print_rotor_summary
from src.bellmouth import bellmouth_design, print_bellmouth_summary
from src.shaft import (rotor_mass, shaft_sizing, bearing_selection, campbell_data,
                       print_shaft_summary, print_bearing_summary, print_campbell_summary)

# 1. Aerodynamic design — blade counts and chords auto-sized
igv = igv_geometry(D_tip=0.900, nu=0.75, N_RPM=3500, phi=0.725, alpha1_deg=0.0)
rotor = meanline_with_igv(igv, PR=1.10, eta_is=0.85)
print_igv_summary(igv)
print_rotor_summary(rotor)

# 2. Inlet design
bell = bellmouth_design(igv, contraction_ratio=4.0, centerbody_fineness=1.2)
print_bellmouth_summary(bell, rotor_dP0_Pa=rotor['dP0_rotor_Pa'])

# 3. Shaft and bearings
rr = rotor_mass(
    r_tip=igv['r_tip_mm']/1000,  r_hub=igv['r_hub_mm']/1000,
    h_blade=igv['h_annulus_mm']/1000,
    B=igv['B_blades'],  chord_mid=igv['rotor_chord_mid_mm']/1000,
)
sr = shaft_sizing(rr, igv, rotor, L_span=0.700, G_grade=1.0)
br = bearing_selection(sr, L10_target_h=20000)
cd = campbell_data(sr, igv)
print_shaft_summary(sr, rr)
print_bearing_summary(br)
print_campbell_summary(cd)
```

---

## Modules planned

| Module | Content |
|---|---|
| `balancing.py` | ISO G-grade derivation, balance plane positions, trial-weight procedure, correction mass limits |
| `blade.py` | NACA 65-series rotor and stator blade design: camber line, thickness distribution, stagger, stacking axis |
| `drive.py` | Motor sizing, VFD selection, flexible coupling, skip-band schedule |
| `instrumentation.py` | Casing static taps, dynamic pressure transducer ring, five-hole probe traverse, DAQ specification |
| `test_matrix.py` | Operating-point sweep procedure, throttle schedule, performance map acquisition |

---

## Dependencies

Install with:

```
pip install -r requirements.txt
```

---

## References

- Dixon, S. L., & Hall, C. A. (2014). *Fluid Mechanics and Thermodynamics of Turbomachinery* (7th ed.). Butterworth-Heinemann.
- Cumpsty, N. A. (2004). *Compressor Aerodynamics*. Krieger Publishing.
- Lieblein, S., Schwenk, F. C., & Broderick, R. L. (1953). Diffusion factor for estimating losses and limiting blade loadings in axial-flow-compressor blade elements. NACA RM E53D01.
- Carter, A. D. S. (1950). The low-speed performance of related aerofoils in cascade. ARC R&M 2816.
- Zweifel, O. (1945). The spacing of turbomachine blading, especially with large angular deflection. *Brown Boveri Review*, 32(12), 436–444.
- Tyler, J. M., & Sofrin, T. G. (1962). Axial flow compressor noise studies. SAE Technical Paper 620532.
- Mehta, R. D., & Bradshaw, P. (1979). Design rules for small low speed wind tunnels. *Aeronautical Journal*, 83(827), 443–449.
- ISO 5801:2017. *Industrial fans — Performance testing using standardized airways*.
- ISO 21940-11:2016. *Mechanical vibration — Rotor balancing — Procedures and tolerances for rotors with rigid behaviour*.
- ISO 281:2007. *Rolling bearings — Dynamic load ratings and rating life*.
- Shigley, J. E., & Mischke, C. R. (2014). *Mechanical Engineering Design* (10th ed.). McGraw-Hill.
- Hoerner, S. F. (1965). *Fluid Dynamic Drag*. Published by the author.
- SKF Group (2018). *SKF General Catalogue*. SKF.



## 👤 Author

**Prof.Dr. Onur Tuncer**  
Aerospace Engineer, Researcher & C++ Systems Developer  
Email: **onur.tuncer@itu.edu.tr**

<p align="left">
  <img src="assets/itu_logo.png" width="180" alt="Istanbul Technical University"/>
</p>
