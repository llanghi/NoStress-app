import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

# --------------------------------------------------------
# Functions
# --------------------------------------------------------
def fracture_stability(sigma_n, tau, pf, IFC, CS):
    """Return Delta p (Pa) for given IFC and CS."""
    pf_fail = sigma_n + (CS - tau) / IFC 
    dp = pf_fail - pf 
    return max(dp, 0)

def column_height(dp, rho_w, rho_g):
    """Return column height (m), rounded to nearest 50."""
    denom = (rho_w - rho_g) * 9.81
    if denom <= 0:
        return np.nan
    h = dp / denom
    return int(round(h / 50.0) * 50)

# --------------------------------------------------------
# Streamlit UI
# --------------------------------------------------------
st.title("Fault Stability & Column Height")

st.sidebar.header("Input Parameters")

# Stress gradients (MPa/km → Pa/m)
sigma_v_grad = st.sidebar.number_input("Vertical stress gradient (MPa/km)", min_value=10.0, max_value=50.0, value=20.0) * 1e6 / 1000
sigma_H_grad = st.sidebar.number_input("Max horizontal stress gradient (MPa/km)", min_value=10.0, max_value=50.0, value=25.0) * 1e6 / 1000
sigma_h_grad = st.sidebar.number_input("Min horizontal stress gradient (MPa/km)", min_value=10.0, max_value=50.0, value=23.0) * 1e6 / 1000
pf_grad      = st.sidebar.number_input("Fluid pressure gradient (MPa/km)", 10.0) * 1e6 / 1000

# Geometry & fluid
depth  = st.sidebar.number_input("Depth (m)", min_value=200.0, value=2000.0)
dip    = st.sidebar.number_input("Fault dip (deg)", min_value=20.0, value=60.0)
strike = st.sidebar.number_input(
    "Fault strike azimuth (deg from North)", 
    min_value=0.0, max_value=360.0, value=90.0
)
az_sigmaH = st.sidebar.number_input(
    "Max horizontal stress azimuth (deg from North)", 
    min_value=0.0, max_value=360.0, value=90.0
)
rho_w  = st.sidebar.number_input("Water density (kg/m³)", min_value=850.0, value=1030.0)
rho_g  = st.sidebar.number_input("Gas density (kg/m³)", min_value=100.0, value=200.0)

# Fault property sliders
IFC = st.sidebar.slider("Internal Friction Coefficient", 0.2, 0.9, 0.6, 0.01)
CS  = st.sidebar.slider("Cohesive Strength (MPa)", 0.0, 10.0, 5.0, 0.1) * 1e6

# --------------------------------------------------------
# Stress calculation (simplified)
# --------------------------------------------------------
sigma_v = sigma_v_grad * depth
sigma_H = sigma_H_grad * depth
sigma_h = sigma_h_grad * depth
pf      = pf_grad * depth

# Assume principal stresses: σ1=σH, σ2=σv, σ3=σh
sigma1, sigma2, sigma3 = sigma_H, sigma_v, sigma_h

# Fault orientation (simplified normal vector in SHmax coords)
dip_rad = np.radians(dip)
dip_dir = strike + 90
dip_dir_rad = np.radians(dip_dir)
az_H_rad = np.radians(az_sigmaH)

# Direction cosines relative to SHmax coordinate system
l = np.sin(dip_rad) * np.cos(dip_dir_rad - az_H_rad)
m = np.sin(dip_rad) * np.sin(dip_dir_rad - az_H_rad)
n = np.cos(dip_rad)
norm = np.sqrt(l**2 + m**2 + n**2)
l, m, n = l/norm, m/norm, n/norm

# Normal and shear stress
sigma_n = sigma1*l**2 + sigma2*m**2 + sigma3*n**2
tau = np.sqrt((sigma1*l)**2 + (sigma2*m)**2 + (sigma3*n)**2 - sigma_n**2)

# Fracture stability & column height
dp = fracture_stability(sigma_n, tau, pf, IFC, CS)
h = column_height(dp, rho_w, rho_g)

# st.subheader("Results")
st.write(f"**Normal stress (σn):** {sigma_n/1e6:.2f} MPa")
st.write(f"**Shear stress (τ):** {tau/1e6:.2f} MPa")
st.write(f"**Fracture stability (Δp):** {dp/1e6:.2f} MPa")
st.success(f"**Max column height: {h:.0f} m (rounded to nearest 50m)**")

# Show unrounded value
h_unrounded = dp / ((rho_w - rho_g) * 9.81)
st.success(f"Max column height: {h_unrounded:.0f} m (unrounded)")

# --------------------------------------------------------
# Lookup table heatmap
# --------------------------------------------------------
# st.subheader("Lookup Table: Column vs Fault Strength")

IFC_vals = np.linspace(0.2, 0.9, 50)
CS_vals  = np.linspace(0, 10e6, 50)
H = np.zeros((len(CS_vals), len(IFC_vals)))

for i, cs in enumerate(CS_vals):
    for j, ifc in enumerate(IFC_vals):
        dp_ij = fracture_stability(sigma_n, tau, pf, ifc, cs)
        H[i, j] = column_height(dp_ij, rho_w, rho_g)

# Use consistent levels every 200 m
levels = np.arange(0, np.nanmax(H)+200, 200)

fig, ax = plt.subplots(figsize=(8,6))
contourf = ax.contourf(IFC_vals, CS_vals/1e6, H, levels=levels, cmap="viridis")

cbar = fig.colorbar(contourf, ax=ax, label="Column height (m)")
cbar.set_ticks(levels)  # ticks every 200 m
cbar.ax.set_yticklabels([f"{int(t)}" if t % 1000 == 0 else "" for t in levels]) # show labels only every 1000 m

# Isolines every 200 m
cs = ax.contour(IFC_vals, CS_vals/1e6, H, levels=levels, colors='k', linewidths=0.5)
ax.clabel(cs, inline=True, fontsize=8, fmt='%d m')

# Mark selected point
ax.plot(IFC, CS/1e6, 'ro', markersize=8, label="column height")
ax.text(IFC+0.01, CS/1e6, f"{h} m", color="black", fontsize=11, weight="bold")
ax.legend(loc="upper left")

ax.set_xlabel("Internal Friction Coefficient")
ax.set_ylabel("Cohesion (MPa)")
#ax.set_title("column height for vs IFC & CS")
ax.set_xlim(0.2, 0.9)  # force X axis from 0.2 to 0.9

st.pyplot(fig)
