# ==============================================================
#  constants.py — Air properties and ISA standard conditions
#
#  Reference:
#    ICAO (2018). Manual of the ICAO Standard Atmosphere (3rd ed.),
#    Doc 7488. International Civil Aviation Organization.
# ==============================================================

# Thermodynamic properties of air
gamma = 1.4
R     = 287.05                      # Specific gas constant [J/(kg·K)]
Cp    = gamma * R / (gamma - 1)     # ≈ 1004.5 J/(kg·K)

# ISA sea-level inlet total conditions
T0_in  = 288.15     # [K]
P0_in  = 101325.0   # [Pa]
rho_in = P0_in / (R * T0_in)   # [kg/m³]