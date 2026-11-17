"""Model data fixed by the brief (section 3), in one place so every script uses the same values."""

SPAN: float = 12.0            # plan side [m]
RISE: float = 2.4             # span / 5 [m]
E: float = 33e9               # C30/37 [Pa]
NU: float = 0.2
RHO: float = 2500.0           # [kg/m^3]
G: float = 9.81               # [m/s^2]
A_H: float = 0.20             # lateral acceleration / g; sensitivity at 0.15 and 0.30
T_REF: float = 0.060          # V60 reference: uniform 60 mm shell [m]
T_BASE_RIBBED: float = 0.050  # base shell of variants B-D [m]
RIB_WIDTH: float = 0.150      # [m]
T_MIN: float = 0.040          # thickness bounds of variant E [m]
T_MAX: float = 0.150
