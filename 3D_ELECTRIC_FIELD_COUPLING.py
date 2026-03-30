"""
================================================================================
  3-D ELECTRIC FIELD SOLVER — MULTI-EMITTER ELECTROSPRAY (SQUARE ARRAY)
================================================================================

PHYSICS
-------
  Governing equation:  ∇²φ = 0  (Laplace, charge-free inter-electrode region)
  This is exact for the electrostatic problem between metallic conductors.

  Geometry:
    • N_x × N_y emitters on a square lattice (pitch p_x, p_y)
    • Each emitter = truncated metallic cone:  tip half-angle α, apex radius r_tip
    • Emitter body held at  V_emit  (high voltage)
    • Extractor plate: grounded plane at  z = L_emit + d_gap

  Inter-emitter coupling (the key physics):
    All emitters are simultaneously present as Dirichlet boundaries.
    Neighbours raise local φ → reduce effective ΔV → weaker E at inner emitters.
    This "shielding" is the dominant non-uniformity in multi-emitter arrays.

NUMERICAL METHOD
----------------
  1. Uniform Cartesian 3-D grid (FD stencil)
  2. Immersed-boundary masking: cone geometry embedded in Cartesian grid
  3. Red-black SOR iterative solver  (ω tuned for optimal spectral radius)
     — proper convergence check: residual measured ONLY on interior free nodes,
       excluding the 1-cell halo adjacent to Dirichlet boundaries
  4. Tip-field analytical correction (hyperboloid model) to recover
     the sub-grid tip enhancement missed by coarse FD
  5. Post-processing:
     • E = −∇φ  via second-order central differences
     • Axial field profiles per emitter
     • Coupling / shielding metrics
     • Full parametric sweep  (pitch / gap ratio)

PLOTS GENERATED
---------------
  01  Axial Ez and potential along every emitter axis + coupling heatmap
  02  XY cross-sections at 4 z-planes (potential + |E| log-scale)
  03  XZ side-view (potential, field-lines, axial Ez component)
  04  Inter-emitter coupling deep-dive  (grouped by corner/edge/centre)
  05  E-field 3-D isosurface approximation via 2-D contour stack
  06  Parametric coupling sweep:  shielding vs. pitch/gap ratio
  07  Electric field gradient  dEz/dz  (ion focusing / stability)
  08  Radial E-field profiles at tip height (transverse spreading)
  09  Equipotential lines — XZ cross-section  (classical field-line diagram)
  10  Convergence history + full parameter table

================================================================================
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from matplotlib.colors import LogNorm, Normalize, BoundaryNorm
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings, time, os, sys
warnings.filterwarnings('ignore')

np.random.seed(42)

# ══════════════════════════════════════════════════════════════════════════════
#  CONFIGURABLE PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════

P = dict(
    # ── Array layout ──────────────────────────────────────────────────────────
    n_x               = 3,           # emitters along x-axis
    n_y               = 3,           # emitters along y-axis
    pitch_x           = 2.0e-3,      # centre-to-centre x spacing  [m]
    pitch_y           = 2.0e-3,      # centre-to-centre y spacing  [m]

    # ── Emitter geometry ──────────────────────────────────────────────────────
    emitter_length    = 1.5e-3,      # needle/cone height  [m]
    tip_half_angle_deg= 15.0,        # cone half-angle at apex  [°]
    tip_radius        = 25e-6,       # apex radius of curvature  [m]

    # ── Electrode voltages & gap ──────────────────────────────────────────────
    V_emitter         = 2500.0,      # emitter potential  [V]
    V_extractor       = 0.0,         # extractor (ground)  [V]
    gap               = 1.0e-3,      # tip-to-extractor distance  [m]

    # ── Flow (reference, not used in electrostatics) ──────────────────────────
    flow_per_emitter  = 5.0e-11,     # volumetric flow rate per emitter  [m³/s]

    # ── Solver ────────────────────────────────────────────────────────────────
    N_grid            = 65,          # grid points along longest axis
    omega             = 1.84,        # SOR relaxation factor  (1 < ω < 2)
    max_iter          = 5000,        # SOR iteration cap
    tol               = 1e-5,        # relative residual (interior free nodes only)

    # ── Parametric sweep ──────────────────────────────────────────────────────
    sweep_ratios      = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0],  # pitch/gap ratios
)

# ══════════════════════════════════════════════════════════════════════════════
#  DARK-THEME COLOUR PALETTE
# ══════════════════════════════════════════════════════════════════════════════

BG  = '#0b0f17'      # figure background
PAN = '#131925'      # axes face
BDR = '#2a3140'      # spine / grid colour
TXT = '#dce6f5'      # primary text
SUB = '#8090a8'      # secondary text / tick labels
ACC = '#4fa3e8'      # accent blue
HOT = '#f06a50'      # hot / corner
YLW = '#e8c44a'      # edge
GRN = '#4ec98e'      # centre / cool
PRP = '#b07ef0'      # purple accent
MAG = '#e868b0'      # magenta accent

def sax(ax, title=None, xl=None, yl=None, fs=8.5, grid=True):
    ax.set_facecolor(PAN)
    ax.tick_params(colors=SUB, labelsize=fs-1.5)
    for sp in ax.spines.values():
        sp.set_edgecolor(BDR)
    if grid:
        ax.grid(True, color=BDR, lw=0.4, alpha=0.6)
        ax.set_axisbelow(True)
    if title: ax.set_title(title, color=TXT, fontsize=fs, pad=5, fontweight='semibold')
    if xl:    ax.set_xlabel(xl,   color=SUB, fontsize=fs-1)
    if yl:    ax.set_ylabel(yl,   color=SUB, fontsize=fs-1)

def cbar(im, ax, label='', fs=7.5):
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label(label, color=SUB, fontsize=fs)
    cb.ax.tick_params(colors=SUB, labelsize=fs-1)
    return cb

def fig_bg(fig):
    fig.patch.set_facecolor(BG)

GCOL = {'corner': HOT, 'edge': YLW, 'centre': GRN}


# ══════════════════════════════════════════════════════════════════════════════
#  GRID CONSTRUCTION
# ══════════════════════════════════════════════════════════════════════════════

def build_grid(P):
    Ne_x, Ne_y = P['n_x'], P['n_y']
    px, py     = P['pitch_x'], P['pitch_y']
    L          = P['emitter_length']
    gap        = P['gap']
    margin_xy  = max(px, py) * 0.70
    margin_z   = gap * 0.30

    xlo = -(Ne_x-1)*px/2 - margin_xy;  xhi = (Ne_x-1)*px/2 + margin_xy
    ylo = -(Ne_y-1)*py/2 - margin_xy;  yhi = (Ne_y-1)*py/2 + margin_xy
    zlo = 0.0;                          zhi = L + gap + margin_z

    Lmax = max(xhi-xlo, yhi-ylo, zhi-zlo)
    N    = P['N_grid']
    ds   = Lmax / N
    Nx   = max(16, round((xhi-xlo)/ds))
    Ny   = max(16, round((yhi-ylo)/ds))
    Nz   = max(16, round((zhi-zlo)/ds))

    x = np.linspace(xlo, xhi, Nx)
    y = np.linspace(ylo, yhi, Ny)
    z = np.linspace(zlo, zhi, Nz)
    dx, dy, dz = x[1]-x[0], y[1]-y[0], z[1]-z[0]

    print(f"  Grid: {Nx}×{Ny}×{Nz}  ({Nx*Ny*Nz/1e6:.2f}M pts)  "
          f"Δs = {ds*1e6:.1f} µm")
    return x, y, z, dx, dy, dz


# ══════════════════════════════════════════════════════════════════════════════
#  GEOMETRY — EMITTER CONE MASK
# ══════════════════════════════════════════════════════════════════════════════

def build_geometry(x, y, z, P):
    Ne_x, Ne_y = P['n_x'], P['n_y']
    px, py     = P['pitch_x'], P['pitch_y']
    L          = P['emitter_length']
    gap        = P['gap']
    alpha      = np.deg2rad(P['tip_half_angle_deg'])
    r_tip      = P['tip_radius']

    xs_c = np.linspace(-(Ne_x-1)*px/2, (Ne_x-1)*px/2, Ne_x)
    ys_c = np.linspace(-(Ne_y-1)*py/2, (Ne_y-1)*py/2, Ne_y)
    centres = [(float(cx), float(cy)) for cx in xs_c for cy in ys_c]

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    body = np.zeros(X.shape, dtype=bool)
    for cx, cy in centres:
        r_xy   = np.sqrt((X-cx)**2 + (Y-cy)**2)
        r_cone = r_tip + np.maximum(0.0, L - Z) * np.tan(alpha)
        body  |= (r_xy <= r_cone) & (Z >= 0.0) & (Z <= L)

    # Extractor plate: one z-layer at z = L+gap
    iz_ext = int(np.argmin(np.abs(z - (L+gap))))
    extr   = np.zeros(X.shape, dtype=bool)
    extr[:, :, iz_ext] = True

    iz_tip  = int(np.argmin(np.abs(z - L)))
    iz_base = int(np.argmin(np.abs(z - 0.0)))

    # ── Classify emitters ────────────────────────────────────────────────────
    cx_u = np.unique([c[0] for c in centres])
    cy_u = np.unique([c[1] for c in centres])

    def classify(cx, cy):
        ex = np.isclose(cx, cx_u[0]) or np.isclose(cx, cx_u[-1])
        ey = np.isclose(cy, cy_u[0]) or np.isclose(cy, cy_u[-1])
        return 'corner' if (ex and ey) else ('edge' if (ex or ey) else 'centre')

    cls_map = {c: classify(*c) for c in centres}

    print(f"  Emitters: {len(centres)} ({Ne_x}×{Ne_y})  |  "
          f"pitch=({px*1e3:.2f},{py*1e3:.2f})mm  |  "
          f"gap={gap*1e3:.2f}mm  |  L={L*1e3:.2f}mm")

    return (X, Y, Z, body, extr, centres, cls_map,
            cx_u, cy_u, iz_tip, iz_ext, iz_base)


# ══════════════════════════════════════════════════════════════════════════════
#  SOR ITERATIVE SOLVER  (red-black, corrected residual)
# ══════════════════════════════════════════════════════════════════════════════

def sor_solve(phi0, fixed, dx, dy, dz, P, label=''):
    phi  = phi0.copy().astype(np.float64)
    om   = P['omega']
    maxi = P['max_iter']
    tol  = P['tol']
    Nx, Ny, Nz = phi.shape
    dx2, dy2, dz2 = dx**2, dy**2, dz**2
    D = 2.0*(1/dx2 + 1/dy2 + 1/dz2)

    # Pre-build interior free-node index mask (excluding halo around fixed)
    from scipy.ndimage import binary_dilation
    fixed_halo = binary_dilation(fixed)   # 1-cell dilation
    interior   = ~fixed_halo              # nodes far from BCs for residual check
    interior[[0,-1],:,:] = False
    interior[:,[0,-1],:] = False
    interior[:,:,[0,-1]] = False

    hist = []
    t0   = time.time()
    tag  = f"  [{label}] " if label else "  "

    for it in range(1, maxi+1):
        for color in (0, 1):
            ix = np.arange(1, Nx-1)
            iy = np.arange(1, Ny-1)
            iz = np.arange(1, Nz-1)
            IX, IY, IZ = np.meshgrid(ix, iy, iz, indexing='ij')
            cm  = ((IX+IY+IZ) % 2 == color)
            i   = IX[cm]; j = IY[cm]; k = IZ[cm]
            fm  = ~fixed[i, j, k]
            i, j, k = i[fm], j[fm], k[fm]
            if len(i) == 0:
                continue
            rhs = ((phi[i-1,j,k]+phi[i+1,j,k])/dx2 +
                   (phi[i,j-1,k]+phi[i,j+1,k])/dy2 +
                   (phi[i,j,k-1]+phi[i,j,k+1])/dz2) / D
            phi[i,j,k] += om * (rhs - phi[i,j,k])

        # Neumann on domain faces (∂φ/∂n = 0)
        phi[ 0,:,:] = phi[ 1,:,:]
        phi[-1,:,:] = phi[-2,:,:]
        phi[:, 0,:] = phi[:, 1,:]
        phi[:,-1,:] = phi[:,-2,:]
        phi[:,:, 0] = phi[:,:, 1]
        phi[:,:,-1] = phi[:,:,-2]
        phi[fixed]  = phi0[fixed]   # re-impose Dirichlet

        if it % 100 == 0 or it == 1:
            # Residual ONLY on interior free nodes (away from BC halo)
            lap = (
                (phi[2:,1:-1,1:-1]-2*phi[1:-1,1:-1,1:-1]+phi[:-2,1:-1,1:-1])/dx2 +
                (phi[1:-1,2:,1:-1]-2*phi[1:-1,1:-1,1:-1]+phi[1:-1,:-2,1:-1])/dy2 +
                (phi[1:-1,1:-1,2:]-2*phi[1:-1,1:-1,1:-1]+phi[1:-1,1:-1,:-2])/dz2
            )
            mask_int = interior[1:-1,1:-1,1:-1]
            if mask_int.any():
                res_abs = np.abs(lap[mask_int]).max()
                phi_ref = np.abs(phi[interior]).max() + 1e-30
                rrel    = res_abs / (phi_ref * D)
            else:
                rrel = 0.0
            hist.append((it, rrel))
            elapsed = time.time()-t0
            print(f"{tag}iter {it:5d}  res={rrel:.3e}  [{elapsed:.1f}s]")
            if rrel < tol:
                print(f"{tag}✓ Converged  (iter={it}, res={rrel:.2e})")
                break
    else:
        print(f"{tag}⚠ Max iter reached (res={hist[-1][1]:.2e})")

    return phi, hist


# ══════════════════════════════════════════════════════════════════════════════
#  ELECTRIC FIELD  E = −∇φ
# ══════════════════════════════════════════════════════════════════════════════

def compute_efield(phi, dx, dy, dz):
    Ex = -np.gradient(phi, dx, axis=0)
    Ey = -np.gradient(phi, dy, axis=1)
    Ez = -np.gradient(phi, dz, axis=2)
    Em = np.sqrt(Ex**2 + Ey**2 + Ez**2)
    return Ex, Ey, Ez, Em


# ══════════════════════════════════════════════════════════════════════════════
#  TIP-FIELD ANALYTICAL CORRECTION
#  For a hyperboloid tip in infinite parallel-plate geometry:
#    E_tip = V / (r_tip · ln(4·d/r_tip))
#  This recovers the enhancement factor the coarse grid misses.
# ══════════════════════════════════════════════════════════════════════════════

def tip_field_corrected(P):
    V  = P['V_emitter'] - P['V_extractor']
    r  = P['tip_radius']
    d  = P['gap']
    # Hyperboloid / needle-plane analytical formula
    E_tip = V / (r * np.log(4*d / r))
    return E_tip


# ══════════════════════════════════════════════════════════════════════════════
#  AXIAL PROFILES  (per emitter)
# ══════════════════════════════════════════════════════════════════════════════

def extract_profiles(phi, Ex, Ey, Ez, Em, x, y, z, centres, cls_map, P):
    L   = P['emitter_length']
    gap = P['gap']
    profiles = {}
    for (cx, cy) in centres:
        ix = int(np.argmin(np.abs(x-cx)))
        iy = int(np.argmin(np.abs(y-cy)))
        # axial gradient  dEz/dz
        dEz_dz = np.gradient(Ez[ix,iy,:], z)
        profiles[(cx,cy)] = dict(
            ix=ix, iy=iy,
            label=f"({cx*1e3:.1f},{cy*1e3:.1f})mm",
            cls=cls_map[(cx,cy)],
            phi    = phi[ix,iy,:],
            Ez     = Ez[ ix,iy,:],
            Em     = Em[ ix,iy,:],
            dEz_dz = dEz_dz,
        )
    return profiles


# ══════════════════════════════════════════════════════════════════════════════
#  PARAMETRIC COUPLING SWEEP
#  For each pitch/gap ratio, re-run a lightweight solve and record shielding.
# ══════════════════════════════════════════════════════════════════════════════

def parametric_sweep(P_base, ratios):
    """
    Vary pitch_x = pitch_y = ratio × gap while keeping everything else fixed.
    Returns list of (ratio, shielding_pct, E_corner_MV, E_centre_MV).
    """
    results = []
    P_sw = P_base.copy()
    # Lighter grid for sweep
    P_sw['N_grid']   = 44
    P_sw['max_iter'] = 2500
    P_sw['tol']      = 5e-4
    P_sw['omega']    = 1.78

    for ratio in ratios:
        pitch = ratio * P_base['gap']
        P_sw['pitch_x'] = pitch
        P_sw['pitch_y'] = pitch
        print(f"\n  ── Sweep: pitch/gap={ratio:.1f}  pitch={pitch*1e3:.2f}mm ──")

        x, y, z, dx, dy, dz = build_grid(P_sw)
        (X3,Y3,Z3,body,extr,centres,cls_map,
         cx_u,cy_u,iz_tip,iz_ext,iz_base) = build_geometry(x,y,z,P_sw)

        L_ref = P_sw['emitter_length']+P_sw['gap']
        phi0  = P_sw['V_emitter'] * np.maximum(0, 1 - z[np.newaxis,np.newaxis,:]/L_ref)
        phi0  = np.broadcast_to(phi0,(len(x),len(y),len(z))).copy()
        phi0[body] = P_sw['V_emitter']
        phi0[extr] = P_sw['V_extractor']
        fixed = body | extr

        phi, _ = sor_solve(phi0, fixed, dx, dy, dz, P_sw, label=f"r={ratio:.1f}")
        _, _, _, Em = compute_efield(phi, dx, dy, dz)

        E_tips = {}
        for (cx,cy) in centres:
            ix_ = int(np.argmin(np.abs(x-cx)))
            iy_ = int(np.argmin(np.abs(y-cy)))
            iz_ = iz_tip
            E_tips[(cx,cy)] = Em[ix_,iy_,iz_]

        all_E   = list(E_tips.values())
        corners = [E_tips[c] for c in centres if cls_map[c]=='corner']
        edges   = [E_tips[c] for c in centres if cls_map[c]=='edge']
        cntrs   = [E_tips[c] for c in centres if cls_map[c]=='centre']
        shld    = (max(all_E)-min(all_E))/max(all_E)*100 if len(all_E)>1 else 0
        results.append(dict(
            ratio    = ratio,
            shielding= shld,
            E_corner = np.mean(corners)*1e-6 if corners else np.mean(all_E)*1e-6,
            E_edge   = np.mean(edges)*1e-6   if edges   else np.nan,
            E_centre = np.mean(cntrs)*1e-6   if cntrs   else np.mean(all_E)*1e-6,
            E_max    = max(all_E)*1e-6,
            E_min    = min(all_E)*1e-6,
        ))
        print(f"  → shielding={shld:.2f}%  "
              f"corner={results[-1]['E_corner']:.3f}  "
              f"centre={results[-1]['E_centre']:.3f} MV/m")

    return results


# ══════════════════════════════════════════════════════════════════════════════
#  PLOTTING FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

OUT = '/mnt/user-data/outputs'

# ─────────────────────────────────────────────────────────────────────────────
# FIG 01  Axial field profiles + coupling heatmap
# ─────────────────────────────────────────────────────────────────────────────
def fig01_axial(profiles, centres, cls_map, cx_u, cy_u, E_tip_corr, P):
    fig = plt.figure(figsize=(18,11), facecolor=BG)
    fig_bg(fig)
    fig.suptitle(
        f"Fig 1 — Axial E-Field & Potential Profiles  "
        f"({P['n_x']}×{P['n_y']} array, V={P['V_emitter']:.0f} V, "
        f"gap={P['gap']*1e3:.1f} mm, pitch={P['pitch_x']*1e3:.1f} mm)",
        color=TXT, fontsize=12, fontweight='bold', y=0.98)

    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.44, wspace=0.38,
                           left=0.06, right=0.97, top=0.93, bottom=0.07)
    ax_ez  = fig.add_subplot(gs[0, :3])
    ax_phi = fig.add_subplot(gs[0,  3])
    ax_map = fig.add_subplot(gs[1, :2])
    ax_tip = fig.add_subplot(gs[1,  2])
    ax_leg = fig.add_subplot(gs[1,  3])

    L, gap = P['emitter_length'], P['gap']
    cpal   = plt.cm.plasma(np.linspace(0.08, 0.92, len(centres)))

    # ── Main axial Ez plot ───────────────────────────────────────────────────
    for idx, (cx,cy) in enumerate(centres):
        pf  = profiles[(cx,cy)]
        zm  = z_arr * 1e3
        Ezp = -pf['Ez'] * 1e-6
        col = cpal[idx]
        lw  = 2.2 if pf['cls']=='corner' else (1.6 if pf['cls']=='edge' else 1.2)
        ls  = '-'  if pf['cls']=='corner' else ('--' if pf['cls']=='edge' else ':')
        ax_ez.plot(zm, Ezp, color=col, lw=lw, ls=ls, alpha=0.88)
        ax_phi.plot(zm, pf['phi'], color=col, lw=lw, ls=ls, alpha=0.88)
        ax_tip.plot(zm[(z_arr>=L-0.05e-3)&(z_arr<=L+gap*1.05)],
                    Ezp[(z_arr>=L-0.05e-3)&(z_arr<=L+gap*1.05)],
                    color=col, lw=lw+0.3, alpha=0.9)

    # Tip / extractor lines
    for a in [ax_ez, ax_tip]:
        a.axvline(L*1e3,       color=ACC, lw=1.0, ls='--', alpha=0.75, zorder=2)
        a.axvline((L+gap)*1e3, color=SUB, lw=0.8, ls='--', alpha=0.65, zorder=2)
        yl = a.get_ylim()
        a.text(L*1e3+0.01,       yl[1]*0.97, 'tip',    color=ACC, fontsize=7, va='top')
        a.text((L+gap)*1e3+0.01, yl[1]*0.97, 'extr.',  color=SUB, fontsize=7, va='top')

    # Analytical tip-field marker
    ax_ez.axhline(E_tip_corr*1e-6, color=MAG, lw=1.0, ls='-.',
                  alpha=0.8, label=f"Analytical E_tip = {E_tip_corr*1e-6:.1f} MV/m")
    ax_ez.legend(facecolor=PAN, edgecolor=BDR, labelcolor=TXT, fontsize=8)

    sax(ax_ez,  "−Ez along emitter axis", "z  [mm]", "−Ez  [MV/m]")
    sax(ax_phi, "Potential φ(z)",          "z  [mm]", "φ  [V]")
    sax(ax_tip, "Gap region zoom",         "z  [mm]", "−Ez  [MV/m]")

    # ── Coupling heatmap ─────────────────────────────────────────────────────
    Ne_x, Ne_y = P['n_x'], P['n_y']
    E_grid = np.zeros((Ne_x, Ne_y))
    px, py = P['pitch_x'], P['pitch_y']
    for ix_e, cxv in enumerate(cx_u):
        for iy_e, cyv in enumerate(cy_u):
            iz_ = int(np.argmin(np.abs(z_arr - L)))
            E_grid[ix_e, iy_e] = profiles[(cxv,cyv)]['Em'][iz_] * 1e-6

    ext_h = [cx_u[0]*1e3-px*0.5e3, cx_u[-1]*1e3+px*0.5e3,
             cy_u[0]*1e3-py*0.5e3, cy_u[-1]*1e3+py*0.5e3]
    im = ax_map.imshow(E_grid.T, origin='lower', cmap='plasma',
                       aspect='equal', interpolation='nearest', extent=ext_h)
    cbar(im, ax_map, '|E|_tip  [MV/m]')
    for ix_e, cxv in enumerate(cx_u):
        for iy_e, cyv in enumerate(cy_u):
            ax_map.text(cxv*1e3, cyv*1e3, f"{E_grid[ix_e,iy_e]:.3f}",
                        ha='center', va='center', color='white',
                        fontsize=8.5, fontweight='bold')
    all_E = E_grid.flatten()
    shld  = (all_E.max()-all_E.min())/all_E.max()*100
    sax(ax_map, f"|E|_tip heatmap — shielding={shld:.1f}%",
        "x  [mm]", "y  [mm]")

    # ── Legend panel ─────────────────────────────────────────────────────────
    ax_leg.set_facecolor(PAN); ax_leg.axis('off')
    items = [
        (HOT, '-',  2.2, 'Corner emitters'),
        (YLW, '--', 1.6, 'Edge emitters'),
        (GRN, ':',  1.2, 'Centre emitters'),
        (MAG, '-.',  1.0, 'Analytical E_tip (hyperboloid)'),
        (ACC, '--',  1.0, 'Tip plane  z = L'),
        (SUB, '--',  0.8, 'Extractor  z = L+gap'),
    ]
    for k,(col,ls,lw_,lbl) in enumerate(items):
        ax_leg.plot([0.05,0.32],[0.88-k*0.14,0.88-k*0.14],
                    color=col, lw=lw_, ls=ls, transform=ax_leg.transAxes)
        ax_leg.text(0.38, 0.88-k*0.14, lbl, color=TXT, fontsize=8,
                    va='center', transform=ax_leg.transAxes)
    ax_leg.set_title("Legend", color=TXT, fontsize=9, pad=5)

    plt.savefig(f'{OUT}/01_axial_profiles.png', dpi=150,
                bbox_inches='tight', facecolor=BG)
    plt.close(fig)
    print("  ✓ Fig 01 saved")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 02  XY cross-sections at 4 z-planes
# ─────────────────────────────────────────────────────────────────────────────
def fig02_xy_slices(phi, Ex, Ey, Ez, Em, x, y, z, centres, P, iz_tip, iz_ext):
    L, gap = P['emitter_length'], P['gap']
    iz_s = [
        int(np.argmin(np.abs(z - L*0.55))),
        iz_tip,
        int(np.argmin(np.abs(z - (L + gap*0.45)))),
        iz_ext,
    ]
    zlabels = [
        f"z={z[iz_s[0]]*1e3:.2f}mm (emitter body)",
        f"z={z[iz_s[1]]*1e3:.2f}mm (tip plane)",
        f"z={z[iz_s[2]]*1e3:.2f}mm (mid-gap)",
        f"z={z[iz_s[3]]*1e3:.2f}mm (extractor)",
    ]
    ext = [x[0]*1e3, x[-1]*1e3, y[0]*1e3, y[-1]*1e3]

    fig = plt.figure(figsize=(20, 11), facecolor=BG)
    fig_bg(fig)
    fig.suptitle("Fig 2 — Electric Field: XY Cross-Sections at Key Heights",
                 color=TXT, fontsize=12, fontweight='bold', y=0.98)
    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.34, wspace=0.28,
                           left=0.05, right=0.97, top=0.93, bottom=0.06)

    for col, (iz_sl, zlbl) in enumerate(zip(iz_s, zlabels)):
        ax_p = fig.add_subplot(gs[0, col])
        ax_E = fig.add_subplot(gs[1, col])

        phi_s  = phi[:, :, iz_sl]
        Em_s   = Em[ :, :, iz_sl] * 1e-6
        Ex_s   = Ex[ :, :, iz_sl]
        Ey_s   = Ey[ :, :, iz_sl]

        # Potential
        im_p = ax_p.imshow(phi_s.T, origin='lower', cmap='RdBu_r',
                           extent=ext, aspect='equal',
                           vmin=P['V_extractor'], vmax=P['V_emitter'])
        cbar(im_p, ax_p, 'φ [V]')

        # Equipotentials overlay
        levels = np.linspace(P['V_extractor']+100, P['V_emitter']-100, 12)
        ax_p.contour(np.linspace(x[0]*1e3,x[-1]*1e3,phi_s.shape[0]),
                     np.linspace(y[0]*1e3,y[-1]*1e3,phi_s.shape[1]),
                     phi_s.T, levels=levels, colors='white', linewidths=0.4,
                     alpha=0.5)
        for (cx,cy) in centres:
            ax_p.plot(cx*1e3, cy*1e3, 'w+', ms=6, mew=1.5, alpha=0.9)
        sax(ax_p, f"φ  —  {zlbl}", "x [mm]", "y [mm]", grid=False)

        # |E| log-scale + streamlines
        Em_c   = np.clip(Em_s, 1e-3, None)
        vmin_E = max(Em_c.min(), 5e-3)
        vmax_E = np.percentile(Em_c, 99.5)
        im_E   = ax_E.imshow(Em_c.T, origin='lower', cmap='inferno',
                             extent=ext, aspect='equal',
                             norm=LogNorm(vmin=vmin_E, vmax=vmax_E))
        cbar(im_E, ax_E, '|E|  [MV/m]')

        sk = max(1, len(x)//15)
        xs_ = x[::sk]*1e3; ys_ = y[::sk]*1e3
        Exd = Ex_s[::sk,::sk].T; Eyd = Ey_s[::sk,::sk].T
        spd = np.sqrt(Exd**2+Eyd**2)+1e-30
        try:
            ax_E.streamplot(xs_, ys_, Exd/spd, Eyd/spd,
                            color='white', linewidth=0.55, density=0.75,
                            arrowsize=0.9, arrowstyle='->')
        except: pass
        for (cx,cy) in centres:
            ax_E.plot(cx*1e3, cy*1e3, 'w+', ms=6, mew=1.5, alpha=0.9)
        sax(ax_E, f"|E| log + field lines  —  {zlbl.split('(')[1][:-1]}",
            "x [mm]", "y [mm]", grid=False)

    plt.savefig(f'{OUT}/02_xy_slices.png', dpi=150,
                bbox_inches='tight', facecolor=BG)
    plt.close(fig)
    print("  ✓ Fig 02 saved")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 03  XZ cross-section
# ─────────────────────────────────────────────────────────────────────────────
def fig03_xz(phi, Ex, Ey, Ez, Em, x, y, z, centres, P):
    L, gap = P['emitter_length'], P['gap']
    iy_c   = int(np.argmin(np.abs(y - 0.0)))
    ext    = [z[0]*1e3, z[-1]*1e3, x[0]*1e3, x[-1]*1e3]

    phi_xz = phi[:, iy_c, :]
    Em_xz  = Em[ :, iy_c, :] * 1e-6
    Ez_xz  = Ez[ :, iy_c, :]
    Ex_xz  = Ex[ :, iy_c, :]

    fig, axes = plt.subplots(1, 3, figsize=(19, 7.5), facecolor=BG)
    fig_bg(fig)
    fig.suptitle("Fig 3 — XZ Cross-Section Through Array Centre (y = 0)",
                 color=TXT, fontsize=12, fontweight='bold')

    # (a) Potential + equipotentials
    ax = axes[0]; ax.set_facecolor(PAN)
    im = ax.imshow(phi_xz, origin='lower', cmap='RdBu_r', extent=ext,
                   aspect='auto', vmin=P['V_extractor'], vmax=P['V_emitter'])
    cbar(im, ax, 'φ  [V]')
    levels = np.linspace(P['V_extractor']+100, P['V_emitter']-100, 20)
    zz = np.linspace(z[0]*1e3, z[-1]*1e3, phi_xz.shape[1])
    xx = np.linspace(x[0]*1e3, x[-1]*1e3, phi_xz.shape[0])
    ax.contour(zz, xx, phi_xz, levels=levels, colors='white',
               linewidths=0.45, alpha=0.65)
    ax.axvline(L*1e3,       color='white', lw=0.9, ls='--', alpha=0.6)
    ax.axvline((L+gap)*1e3, color=ACC,     lw=0.9, ls='--', alpha=0.7)
    ax.text(L*1e3+0.02,       x[0]*1e3*0.9+x[-1]*1e3*0.1, 'tip',   color='white', fontsize=7.5)
    ax.text((L+gap)*1e3+0.02, x[0]*1e3*0.9+x[-1]*1e3*0.1, 'extr.', color=ACC,     fontsize=7.5)
    sax(ax, "Potential φ(x,z) + equipotentials", "z [mm]", "x [mm]", grid=False)

    # (b) |E| log + field lines
    ax = axes[1]; ax.set_facecolor(PAN)
    Em_c = np.clip(Em_xz, 5e-3, None)
    im   = ax.imshow(Em_c, origin='lower', cmap='inferno', extent=ext,
                     aspect='auto',
                     norm=LogNorm(vmin=max(Em_c.min(), 5e-3), vmax=Em_c.max()))
    cbar(im, ax, '|E|  [MV/m]')
    skz = max(1, len(z)//22); skx = max(1, len(x)//15)
    zs_ = z[::skz]*1e3; xs_ = x[::skx]*1e3
    Ezd = Ez_xz[::skx,::skz]; Exd = Ex_xz[::skx,::skz]
    spd = np.sqrt(Ezd**2+Exd**2)+1e-30
    try:
        ax.streamplot(zs_, xs_, Ezd/spd, Exd/spd,
                      color='white', linewidth=0.6, density=1.1,
                      arrowsize=0.9, arrowstyle='->')
    except: pass
    ax.axvline(L*1e3,       color='white', lw=0.9, ls='--', alpha=0.55)
    ax.axvline((L+gap)*1e3, color=ACC,     lw=0.9, ls='--', alpha=0.7)
    sax(ax, "|E| (log scale) + field lines", "z [mm]", "x [mm]", grid=False)

    # (c) Axial Ez (signed) — shows direction and structure
    ax = axes[2]; ax.set_facecolor(PAN)
    Ezp = -Ez_xz * 1e-6
    vl  = np.percentile(np.abs(Ezp), 98.5)
    im  = ax.imshow(Ezp, origin='lower', cmap='seismic', extent=ext,
                    aspect='auto', vmin=-vl, vmax=vl)
    cbar(im, ax, '−Ez  [MV/m]')
    ax.contour(zz, xx, Ezp, levels=np.linspace(-vl*0.9, vl*0.9, 18),
               colors='white', linewidths=0.35, alpha=0.45)
    ax.axvline(L*1e3,       color='white', lw=0.9, ls='--', alpha=0.6)
    ax.axvline((L+gap)*1e3, color=ACC,     lw=0.9, ls='--', alpha=0.7)
    sax(ax, "Axial field component −Ez(x,z)", "z [mm]", "x [mm]", grid=False)

    plt.tight_layout(rect=[0,0,1,0.94])
    plt.savefig(f'{OUT}/03_xz_cross_section.png', dpi=150,
                bbox_inches='tight', facecolor=BG)
    plt.close(fig)
    print("  ✓ Fig 03 saved")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 04  Coupling deep-dive (grouped analysis)
# ─────────────────────────────────────────────────────────────────────────────
def fig04_coupling(profiles, centres, cls_map, cx_u, cy_u, P):
    L, gap = P['emitter_length'], P['gap']
    groups = {'corner':[], 'edge':[], 'centre':[]}
    for (cx,cy), pf in profiles.items():
        groups[pf['cls']].append(pf)

    fig = plt.figure(figsize=(18, 12), facecolor=BG)
    fig_bg(fig)
    fig.suptitle("Fig 4 — Inter-Emitter Coupling: Shielding & Field Non-Uniformity",
                 color=TXT, fontsize=12, fontweight='bold', y=0.98)
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.44, wspace=0.36,
                           left=0.07, right=0.97, top=0.93, bottom=0.07)

    ax1 = fig.add_subplot(gs[0, :2])  # mean Ez per class
    ax2 = fig.add_subplot(gs[0,  2])  # gap-region zoom
    ax3 = fig.add_subplot(gs[1,  0])  # |E|_tip vs radial distance
    ax4 = fig.add_subplot(gs[1,  1])  # potential along central row
    ax5 = fig.add_subplot(gs[1,  2])  # dEz/dz (axial gradient)

    zm = z_arr * 1e3
    mask_gap = (z_arr >= L-0.1e-3) & (z_arr <= L+gap*1.05)

    # ── (ax1) Mean Ez by class + spread ─────────────────────────────────────
    for g, profs in groups.items():
        if not profs: continue
        col  = GCOL[g]
        Ezs  = np.array([-pf_['Ez'] for pf_ in profs]) * 1e-6
        mean = Ezs.mean(0)
        ax1.plot(zm, mean, color=col, lw=2.2,
                 label=f"{g.capitalize()} (n={len(profs)})")
        if len(profs) > 1:
            ax1.fill_between(zm, Ezs.min(0), Ezs.max(0), color=col, alpha=0.13)

    ax1.axvline(L*1e3,       color=ACC, lw=1.0, ls='--', alpha=0.7)
    ax1.axvline((L+gap)*1e3, color=SUB, lw=0.8, ls='--', alpha=0.6)
    ax1.legend(facecolor=PAN, edgecolor=BDR, labelcolor=TXT, fontsize=9.5)
    sax(ax1, "Mean −Ez by emitter class  (shading = intra-class spread)",
        "z  [mm]", "−Ez  [MV/m]")

    # ── (ax2) Gap-only zoom ──────────────────────────────────────────────────
    for g, profs in groups.items():
        if not profs: continue
        col = GCOL[g]
        for pf_ in profs:
            ax2.plot(zm[mask_gap], (-pf_['Ez'][mask_gap])*1e-6,
                     color=col, lw=1.4, alpha=0.75)
    ax2.axvline(L*1e3,       color=ACC, lw=1.0, ls='--', alpha=0.7)
    ax2.axvline((L+gap)*1e3, color=SUB, lw=0.8, ls='--', alpha=0.6)
    sax(ax2, "Gap region: all emitters", "z  [mm]", "−Ez  [MV/m]")

    # ── (ax3) |E|_tip vs radial distance from centre ─────────────────────────
    iz_L = int(np.argmin(np.abs(z_arr - L)))
    for (cx,cy), pf in profiles.items():
        d  = np.sqrt(cx**2+cy**2)*1e3
        Et = pf['Em'][iz_L]*1e-6
        ax3.scatter(d, Et, color=GCOL[pf['cls']], s=80, zorder=5,
                    edgecolors='white', lw=0.6)
        ax3.annotate(pf['label'], (d,Et), color=SUB, fontsize=6.5,
                     textcoords='offset points', xytext=(4,3))
    dists_ = [np.sqrt(c[0]**2+c[1]**2)*1e3 for c in profiles]
    Ets_   = [profiles[c]['Em'][iz_L]*1e-6 for c in profiles]
    if len(dists_) > 2:
        coef = np.polyfit(dists_, Ets_, 1)
        dl   = np.linspace(0, max(dists_)*1.05, 60)
        ax3.plot(dl, np.polyval(coef,dl), color=SUB, lw=1.2, ls='--', alpha=0.6)
    sax(ax3, "|E|_tip vs. radial position",
        "r from array centre  [mm]", "|E|_tip  [MV/m]")

    # ── (ax4) Potential along central row (y≈0) at tip height ───────────────
    row = sorted([(cx, pf) for (cx,cy), pf in profiles.items()
                  if np.isclose(cy, 0, atol=P['pitch_y']*0.45)],
                 key=lambda t: t[0])
    if row:
        xs_r  = [t[0]*1e3 for t in row]
        phi_r = [t[1]['phi'][iz_L] for t in row]
        Ez_r  = [-t[1]['Ez'][iz_L]*1e-6 for t in row]
        ax4.plot(xs_r, phi_r, 'o-', color=ACC, lw=2.2, ms=8, label='φ  [V]')
        ax4b = ax4.twinx()
        ax4b.plot(xs_r, Ez_r, 's--', color=HOT, lw=2.2, ms=8, label='−Ez  [MV/m]')
        ax4b.tick_params(colors=SUB, labelsize=7)
        ax4b.set_ylabel("−Ez  [MV/m]", color=SUB, fontsize=8)
        ax4b.spines['right'].set_edgecolor(BDR)
        ax4.legend(loc='upper left',  facecolor=PAN, edgecolor=BDR,
                   labelcolor=TXT, fontsize=8)
        ax4b.legend(loc='upper right', facecolor=PAN, edgecolor=BDR,
                    labelcolor=TXT, fontsize=8)
    sax(ax4, "φ & −Ez along central row at tip height", "x  [mm]", "φ  [V]")

    # ── (ax5) dEz/dz gradient (focusing strength) per class ─────────────────
    for g, profs in groups.items():
        if not profs: continue
        col = GCOL[g]
        grad_all = np.array([pf_['dEz_dz'] for pf_ in profs]) * 1e-9   # GV/m²
        mean_g   = grad_all.mean(0)
        ax5.plot(zm, mean_g, color=col, lw=2.2, label=g.capitalize())
        if len(profs) > 1:
            ax5.fill_between(zm, grad_all.min(0), grad_all.max(0),
                             color=col, alpha=0.12)
    ax5.axvline(L*1e3,       color=ACC, lw=1.0, ls='--', alpha=0.7)
    ax5.axvline((L+gap)*1e3, color=SUB, lw=0.8, ls='--', alpha=0.6)
    ax5.axhline(0, color=BDR, lw=0.7)
    ax5.legend(facecolor=PAN, edgecolor=BDR, labelcolor=TXT, fontsize=9)
    sax(ax5, "Axial E-field gradient dEz/dz  (ion focusing)",
        "z  [mm]", "dEz/dz  [GV/m²]")

    note = (
        f"Shielding: neighbours raise φ → weaken ΔV → reduce E at inner emitters. "
        f"pitch/gap = {P['pitch_x']/P['gap']:.1f} | "
        f"Corner/centre E ratio = "
        f"{max(Ets_)/min(Ets_):.3f}"
    )
    fig.text(0.5, 0.005, note, ha='center', color=SUB, fontsize=8, style='italic')

    plt.savefig(f'{OUT}/04_coupling_analysis.png', dpi=150,
                bbox_inches='tight', facecolor=BG)
    plt.close(fig)
    print("  ✓ Fig 04 saved")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 05  2-D contour stack approximating 3-D isosurface
# ─────────────────────────────────────────────────────────────────────────────
def fig05_isosurface_stack(phi, Em, x, y, z, centres, P):
    L, gap = P['emitter_length'], P['gap']
    # Choose z-planes spanning emitter to extractor
    z_planes = np.linspace(L*0.8, L+gap, 8)
    iz_planes = [int(np.argmin(np.abs(z-zp))) for zp in z_planes]

    fig = plt.figure(figsize=(20, 10), facecolor=BG)
    fig_bg(fig)
    fig.suptitle("Fig 5 — 3-D Field Structure: |E| Contour Stack (z-slices)",
                 color=TXT, fontsize=12, fontweight='bold', y=0.98)

    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.32, wspace=0.22,
                           left=0.04, right=0.97, top=0.93, bottom=0.05)

    ext = [x[0]*1e3, x[-1]*1e3, y[0]*1e3, y[-1]*1e3]
    # Global colour scale
    Em_all = np.array([Em[:,:,iz]*1e-6 for iz in iz_planes])
    vmin   = max(Em_all.min(), 1e-2)
    vmax   = np.percentile(Em_all, 99.5)

    for k, (iz_s, zp) in enumerate(zip(iz_planes, z_planes)):
        row, col_ = divmod(k, 4)
        ax = fig.add_subplot(gs[row, col_])
        ax.set_facecolor(PAN)
        Em_s = np.clip(Em[:,:,iz_s]*1e-6, vmin, None)
        im   = ax.imshow(Em_s.T, origin='lower', cmap='inferno',
                         extent=ext, aspect='equal',
                         norm=LogNorm(vmin=vmin, vmax=vmax))
        # |E| iso-contours
        levels = np.logspace(np.log10(vmin), np.log10(vmax), 8)
        ax.contour(np.linspace(x[0]*1e3,x[-1]*1e3,Em_s.shape[0]),
                   np.linspace(y[0]*1e3,y[-1]*1e3,Em_s.shape[1]),
                   Em_s.T, levels=levels, colors='white',
                   linewidths=0.4, alpha=0.55)
        # Emitter tips
        for (cx,cy) in centres:
            ax.plot(cx*1e3, cy*1e3, '+', color=ACC, ms=6, mew=1.5)
        z_rel = (zp - L)*1e3
        region = "body" if zp < L else ("gap" if zp < L+gap else "extr.")
        ax.set_title(f"z={zp*1e3:.2f}mm  ({region})",
                     color=TXT, fontsize=8, pad=3)
        ax.tick_params(colors=SUB, labelsize=6)
        for sp in ax.spines.values(): sp.set_edgecolor(BDR)

    # Shared colorbar
    fig.subplots_adjust(right=0.93)
    cax = fig.add_axes([0.945, 0.08, 0.012, 0.83])
    sm  = plt.cm.ScalarMappable(cmap='inferno', norm=LogNorm(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label('|E|  [MV/m]', color=SUB, fontsize=8.5)
    cb.ax.tick_params(colors=SUB, labelsize=7.5)

    plt.savefig(f'{OUT}/05_contour_stack.png', dpi=150,
                bbox_inches='tight', facecolor=BG)
    plt.close(fig)
    print("  ✓ Fig 05 saved")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 06  Parametric sweep: shielding vs. pitch/gap
# ─────────────────────────────────────────────────────────────────────────────
def fig06_parametric(sweep_results, P_base):
    fig = plt.figure(figsize=(16, 9), facecolor=BG)
    fig_bg(fig)
    fig.suptitle("Fig 6 — Parametric Sweep: Inter-Emitter Shielding vs. Pitch/Gap Ratio",
                 color=TXT, fontsize=12, fontweight='bold', y=0.98)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.36,
                           left=0.08, right=0.97, top=0.93, bottom=0.08)

    ratios  = [r['ratio']    for r in sweep_results]
    shields = [r['shielding'] for r in sweep_results]
    E_cor   = [r['E_corner']  for r in sweep_results]
    E_ctr   = [r['E_centre']  for r in sweep_results]
    E_edg   = [r.get('E_edge', np.nan) for r in sweep_results]
    E_max   = [r['E_max']     for r in sweep_results]
    E_min   = [r['E_min']     for r in sweep_results]

    # (a) Shielding %
    ax = fig.add_subplot(gs[0,0])
    ax.fill_between(ratios, shields, alpha=0.18, color=ACC)
    ax.plot(ratios, shields, 'o-', color=ACC, lw=2.2, ms=8)
    ax.axvline(P_base['pitch_x']/P_base['gap'], color=HOT, lw=1.2,
               ls='--', label=f"baseline r={P_base['pitch_x']/P_base['gap']:.1f}")
    ax.set_ylim(bottom=0)
    ax.legend(facecolor=PAN, edgecolor=BDR, labelcolor=TXT, fontsize=8.5)
    sax(ax, "Shielding (field non-uniformity) vs. pitch/gap",
        "pitch / gap", "shielding  [%]")

    # (b) Absolute tip E per class
    ax = fig.add_subplot(gs[0,1])
    ax.plot(ratios, E_cor, 'o-', color=HOT, lw=2.2, ms=8, label='Corner')
    ax.plot(ratios, E_edg, 's--',color=YLW, lw=1.8, ms=7, label='Edge')
    ax.plot(ratios, E_ctr, '^:', color=GRN, lw=1.8, ms=7, label='Centre')
    ax.axvline(P_base['pitch_x']/P_base['gap'], color=BDR, lw=1, ls='--', alpha=0.7)
    ax.legend(facecolor=PAN, edgecolor=BDR, labelcolor=TXT, fontsize=8.5)
    sax(ax, "Tip |E| by emitter class vs. pitch/gap",
        "pitch / gap", "|E|_tip  [MV/m]")

    # (c) Corner / centre ratio
    ratios_arr = np.array(ratios)
    cor_arr    = np.array(E_cor)
    ctr_arr    = np.array(E_ctr)
    ratio_EC   = cor_arr / np.where(ctr_arr>0, ctr_arr, np.nan)
    ax = fig.add_subplot(gs[1,0])
    ax.fill_between(ratios, ratio_EC, 1.0, alpha=0.15, color=MAG)
    ax.plot(ratios, ratio_EC, 'D-', color=MAG, lw=2.2, ms=8)
    ax.axhline(1.0, color=SUB, lw=0.8, ls='--', alpha=0.6)
    ax.axvline(P_base['pitch_x']/P_base['gap'], color=HOT, lw=1.2, ls='--')
    sax(ax, "Corner / centre E-field ratio  (coupling severity)",
        "pitch / gap", "E_corner / E_centre")

    # (d) Field uniformity band
    ax = fig.add_subplot(gs[1,1])
    ax.fill_between(ratios, E_min, E_max, alpha=0.2, color=ACC,
                    label='Min–max band')
    ax.plot(ratios, E_max, 'o-', color=HOT, lw=2, ms=7, label='Max (corner)')
    ax.plot(ratios, E_min, 's-', color=GRN, lw=2, ms=7, label='Min (centre)')
    ax.axvline(P_base['pitch_x']/P_base['gap'], color=BDR, lw=1, ls='--', alpha=0.7)
    ax.legend(facecolor=PAN, edgecolor=BDR, labelcolor=TXT, fontsize=8.5)
    sax(ax, "Field uniformity band across array", "pitch / gap", "|E|_tip  [MV/m]")

    fig.text(0.5, 0.005,
             "Larger pitch/gap → weaker coupling → more uniform field. "
             "Penalty: larger footprint.",
             ha='center', color=SUB, fontsize=8.5, style='italic')

    plt.savefig(f'{OUT}/06_parametric_sweep.png', dpi=150,
                bbox_inches='tight', facecolor=BG)
    plt.close(fig)
    print("  ✓ Fig 06 saved")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 07  E-field gradient dEz/dz — ion focusing
# ─────────────────────────────────────────────────────────────────────────────
def fig07_gradient(profiles, centres, cls_map, phi, Ez, x, y, z, P):
    L, gap = P['emitter_length'], P['gap']
    iy_c   = int(np.argmin(np.abs(y - 0.0)))

    fig = plt.figure(figsize=(17, 10), facecolor=BG)
    fig_bg(fig)
    fig.suptitle("Fig 7 — Axial E-Field Gradient dEz/dz  (Ion Acceleration & Focusing)",
                 color=TXT, fontsize=12, fontweight='bold', y=0.98)
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.44, wspace=0.36,
                           left=0.07, right=0.97, top=0.93, bottom=0.07)

    ax1 = fig.add_subplot(gs[0, :2])   # dEz/dz(z) by class
    ax2 = fig.add_subplot(gs[0,  2])   # dEz/dz at tip vs. x position
    ax3 = fig.add_subplot(gs[1, :])    # 2D map dEz/dz in XZ

    zm = z_arr * 1e3
    groups = {'corner':[], 'edge':[], 'centre':[]}
    for (cx,cy), pf in profiles.items():
        groups[pf['cls']].append(pf)

    for g, profs in groups.items():
        if not profs: continue
        col  = GCOL[g]
        grads = np.array([pf_['dEz_dz'] for pf_ in profs]) * 1e-9   # GV/m²
        mean  = grads.mean(0)
        ax1.plot(zm, mean, color=col, lw=2.2, label=g.capitalize())
        if len(profs) > 1:
            ax1.fill_between(zm, grads.min(0), grads.max(0), color=col, alpha=0.12)
    ax1.axvline(L*1e3,       color=ACC, lw=1.0, ls='--', alpha=0.7)
    ax1.axvline((L+gap)*1e3, color=SUB, lw=0.8, ls='--', alpha=0.6)
    ax1.axhline(0, color=BDR, lw=0.6)
    ax1.legend(facecolor=PAN, edgecolor=BDR, labelcolor=TXT, fontsize=9)
    sax(ax1, "Axial gradient dEz/dz — positive → accelerating, negative → decelerating",
        "z  [mm]", "dEz/dz  [GV/m²]")

    # dEz/dz at tip height along central row
    iz_L = int(np.argmin(np.abs(z_arr - L)))
    row  = sorted([(cx, pf) for (cx,cy), pf in profiles.items()
                   if np.isclose(cy, 0, atol=P['pitch_y']*0.45)],
                  key=lambda t: t[0])
    if row:
        xs_r = [t[0]*1e3 for t in row]
        gr_r = [t[1]['dEz_dz'][iz_L]*1e-9 for t in row]
        ax2.bar(xs_r, gr_r, width=P['pitch_x']*0.6*1e3, color=PRP, alpha=0.8,
                edgecolor=BDR)
        ax2.axhline(0, color=SUB, lw=0.7)
    sax(ax2, "dEz/dz at tip height (central row)", "x  [mm]", "dEz/dz  [GV/m²]")

    # 2D dEz/dz map in XZ
    ax3.set_facecolor(PAN)
    dEz_xz = np.gradient(Ez[:,iy_c,:], z, axis=1) * 1e-9   # GV/m²
    ext    = [z[0]*1e3, z[-1]*1e3, x[0]*1e3, x[-1]*1e3]
    vl     = np.percentile(np.abs(dEz_xz), 98)
    im     = ax3.imshow(dEz_xz, origin='lower', cmap='PiYG',
                        extent=ext, aspect='auto', vmin=-vl, vmax=vl)
    cbar(im, ax3, 'dEz/dz  [GV/m²]')
    ax3.axvline(L*1e3,       color='white', lw=0.9, ls='--', alpha=0.55)
    ax3.axvline((L+gap)*1e3, color=ACC,     lw=0.9, ls='--', alpha=0.7)
    # emitter positions
    for (cx,_) in [(c[0],c[1]) for c in centres if np.isclose(c[1],0,atol=P['pitch_y']*0.4)]:
        ax3.axhline(cx*1e3, color=BDR, lw=0.4, alpha=0.5)
    sax(ax3, "dEz/dz in XZ plane  (green=accelerating, pink=decelerating)",
        "z  [mm]", "x  [mm]", grid=False)

    plt.savefig(f'{OUT}/07_efield_gradient.png', dpi=150,
                bbox_inches='tight', facecolor=BG)
    plt.close(fig)
    print("  ✓ Fig 07 saved")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 08  Radial E-field profiles at tip height
# ─────────────────────────────────────────────────────────────────────────────
def fig08_radial(phi, Ex, Ey, Em, x, y, z, centres, cls_map, P):
    L, gap = P['emitter_length'], P['gap']
    iz_L   = int(np.argmin(np.abs(z - L)))
    iz_mid = int(np.argmin(np.abs(z - (L+gap*0.5))))

    # Radial profiles: from each emitter centre outward to pitch/2
    r_max  = min(P['pitch_x'], P['pitch_y']) * 0.48
    Nr     = 80
    r_arr  = np.linspace(0, r_max, Nr)

    fig = plt.figure(figsize=(17, 10), facecolor=BG)
    fig_bg(fig)
    fig.suptitle("Fig 8 — Radial E-Field Profiles at Tip & Mid-Gap Heights",
                 color=TXT, fontsize=12, fontweight='bold', y=0.98)
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.44, wspace=0.36,
                           left=0.07, right=0.97, top=0.93, bottom=0.07)

    ax1 = fig.add_subplot(gs[0, :2])  # Er(r) at tip height per class
    ax2 = fig.add_subplot(gs[0,  2])  # Ez(r) transverse variation at tip
    ax3 = fig.add_subplot(gs[1, :2])  # Er(r) at mid-gap
    ax4 = fig.add_subplot(gs[1,  2])  # contour |E| at tip height

    for (cx,cy), pf in profiles_glob.items():
        ix0 = pf['ix']; iy0 = pf['iy']
        col = GCOL[pf['cls']]
        lw  = 2.0 if pf['cls']=='corner' else 1.3

        Er_tip = [];  Ez_tip = []
        Er_mid = []
        for ri in range(Nr):
            # Radial angle: average over 4 cardinal directions
            vals_r_tip=[]; vals_z_tip=[]; vals_r_mid=[]
            for ang in np.linspace(0, 2*np.pi, 8, endpoint=False):
                dxr = r_arr[ri]*np.cos(ang)
                dyr = r_arr[ri]*np.sin(ang)
                ix_ = int(np.argmin(np.abs(x-(cx+dxr))))
                iy_ = int(np.argmin(np.abs(y-(cy+dyr))))
                ix_ = np.clip(ix_, 0, len(x)-1)
                iy_ = np.clip(iy_, 0, len(y)-1)
                Er_  = np.sqrt(Ex[ix_,iy_,iz_L]**2+Ey[ix_,iy_,iz_L]**2)
                vals_r_tip.append(Er_*1e-6)
                vals_z_tip.append(-Em[ix_,iy_,iz_L]*1e-6)  # total E
                vals_r_mid.append(np.sqrt(Ex[ix_,iy_,iz_mid]**2+Ey[ix_,iy_,iz_mid]**2)*1e-6)
            Er_tip.append(np.mean(vals_r_tip))
            Ez_tip.append(np.mean(vals_z_tip))
            Er_mid.append(np.mean(vals_r_mid))

        ax1.plot(r_arr*1e3, Er_tip, color=col, lw=lw, alpha=0.8, label=pf['label'] if lw>1.5 else None)
        ax2.plot(r_arr*1e3, Ez_tip, color=col, lw=lw, alpha=0.8)
        ax3.plot(r_arr*1e3, Er_mid, color=col, lw=lw, alpha=0.8)

    for ax_, ttl, yl_ in [
        (ax1, "Radial E_r(r) at tip height  z=L", "E_r  [MV/m]"),
        (ax2, "|E|(r) at tip height  z=L",         "|E|  [MV/m]"),
        (ax3, "Radial E_r(r) at mid-gap  z=L+gap/2","E_r  [MV/m]"),
    ]:
        ax_.axvline(P['pitch_x']*0.5e3, color=SUB, lw=0.7, ls='--', alpha=0.5)
        sax(ax_, ttl, "r from emitter axis  [mm]", yl_)

    ax1.legend(facecolor=PAN, edgecolor=BDR, labelcolor=TXT, fontsize=7.5, ncol=2)

    # |E| at tip height (XY map)
    ax4.set_facecolor(PAN)
    Em_tip = np.clip(Em[:,:,iz_L]*1e-6, 1e-2, None)
    ext    = [x[0]*1e3, x[-1]*1e3, y[0]*1e3, y[-1]*1e3]
    im     = ax4.imshow(Em_tip.T, origin='lower', cmap='inferno',
                        extent=ext, aspect='equal',
                        norm=LogNorm(vmin=max(Em_tip.min(),0.01), vmax=Em_tip.max()))
    cbar(im, ax4, '|E|  [MV/m]')
    for (cx,cy) in centres:
        ax4.plot(cx*1e3, cy*1e3, 'w+', ms=6, mew=1.5)
    sax(ax4, "|E| at tip height (XY)", "x  [mm]", "y  [mm]", grid=False)

    plt.savefig(f'{OUT}/08_radial_profiles.png', dpi=150,
                bbox_inches='tight', facecolor=BG)
    plt.close(fig)
    print("  ✓ Fig 08 saved")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 09  Classical equipotential field-line diagram
# ─────────────────────────────────────────────────────────────────────────────
def fig09_equipotential(phi, Ex, Ez, x, y, z, centres, P):
    L, gap = P['emitter_length'], P['gap']
    iy_c   = int(np.argmin(np.abs(y - 0.0)))

    phi_xz = phi[:, iy_c, :]
    Ex_xz  = Ex[ :, iy_c, :]
    Ez_xz  = Ez[ :, iy_c, :]

    fig, axes = plt.subplots(1, 2, figsize=(18, 9), facecolor=BG)
    fig_bg(fig)
    fig.suptitle("Fig 9 — Classical Equipotential & Field-Line Diagram  (XZ Plane)",
                 color=TXT, fontsize=12, fontweight='bold', y=0.98)

    ext  = [z[0]*1e3, z[-1]*1e3, x[0]*1e3, x[-1]*1e3]
    zz   = np.linspace(z[0]*1e3, z[-1]*1e3, phi_xz.shape[1])
    xx   = np.linspace(x[0]*1e3, x[-1]*1e3, phi_xz.shape[0])

    # (a) Dense equipotentials
    ax = axes[0]; ax.set_facecolor('#0a0f18')
    levels_phi = np.linspace(P['V_extractor']+50, P['V_emitter']-50, 35)
    cs = ax.contour(zz, xx, phi_xz, levels=levels_phi,
                    cmap='coolwarm', linewidths=0.7)
    ax.clabel(cs, inline=True, fontsize=5.5, colors=SUB, fmt='%g V',
              manual=False)
    ax.axvline(L*1e3,       color=GRN, lw=1.2, ls='-', alpha=0.6)
    ax.axvline((L+gap)*1e3, color=HOT, lw=1.2, ls='-', alpha=0.6)
    ax.text(L*1e3+0.02, xx[-1]*0.95, 'tip', color=GRN, fontsize=8, va='top')
    ax.text((L+gap)*1e3+0.02, xx[-1]*0.95, 'extr.', color=HOT, fontsize=8, va='top')
    sax(ax, "Equipotential Lines (dense)", "z  [mm]", "x  [mm]", grid=False)

    # (b) Field lines (streamplot) over potential shading
    ax = axes[1]; ax.set_facecolor('#0a0f18')
    im = ax.imshow(phi_xz, origin='lower', cmap='twilight_shifted',
                   extent=ext, aspect='auto', alpha=0.55,
                   vmin=P['V_extractor'], vmax=P['V_emitter'])
    cbar(im, ax, 'φ  [V]')
    # Dense field lines
    skz = max(1, len(z)//30); skx = max(1, len(x)//20)
    zs_ = z[::skz]*1e3; xs_ = x[::skx]*1e3
    Ezd = Ez_xz[::skx,::skz]; Exd = Ex_xz[::skx,::skz]
    spd = np.sqrt(Ezd**2+Exd**2)+1e-30
    try:
        strm = ax.streamplot(zs_, xs_, Ezd/spd, Exd/spd,
                             color='white', linewidth=0.7, density=1.4,
                             arrowsize=1.0, arrowstyle='->')
    except: pass
    # Equipotential contours (sparse)
    ax.contour(zz, xx, phi_xz, levels=levels_phi[::3],
               colors=ACC, linewidths=0.5, alpha=0.45)
    ax.axvline(L*1e3,       color=GRN, lw=1.2, ls='-', alpha=0.6)
    ax.axvline((L+gap)*1e3, color=HOT, lw=1.2, ls='-', alpha=0.6)
    sax(ax, "Field lines (white) + equipotentials (blue) over φ",
        "z  [mm]", "x  [mm]", grid=False)

    plt.tight_layout(rect=[0,0,1,0.95])
    plt.savefig(f'{OUT}/09_equipotential_fieldlines.png', dpi=150,
                bbox_inches='tight', facecolor=BG)
    plt.close(fig)
    print("  ✓ Fig 09 saved")


# ─────────────────────────────────────────────────────────────────────────────
# FIG 10  Convergence + parameter table
# ─────────────────────────────────────────────────────────────────────────────
def fig10_convergence(hist, P, profiles, E_tip_corr):
    fig = plt.figure(figsize=(17, 7.5), facecolor=BG)
    fig_bg(fig)
    fig.suptitle("Fig 10 — SOR Convergence History & Simulation Summary",
                 color=TXT, fontsize=12, fontweight='bold', y=0.98)
    gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.38,
                           left=0.06, right=0.97, top=0.90, bottom=0.10)

    # (a) Convergence
    ax = fig.add_subplot(gs[0, :2])
    iters = [h[0] for h in hist]
    ress  = [h[1] for h in hist]
    ax.fill_between(iters, ress, P['tol'], where=[r>P['tol'] for r in ress],
                    alpha=0.12, color=ACC)
    ax.semilogy(iters, ress, 'o-', color=ACC, lw=2.2, ms=5, label='Residual')
    ax.axhline(P['tol'], color=HOT, lw=1.2, ls='--',
               label=f"tol = {P['tol']:.0e}")
    conv_it = next((h[0] for h in hist if h[1] <= P['tol']), iters[-1])
    ax.axvline(conv_it, color=GRN, lw=1, ls=':', alpha=0.7,
               label=f"Converged @ iter {conv_it}")
    ax.legend(facecolor=PAN, edgecolor=BDR, labelcolor=TXT, fontsize=9)
    sax(ax, "SOR Relative Residual (interior free nodes only)",
        "Iteration", "Relative Residual  ||∇²φ||_∞ / (||φ||·D)")

    # (b) Parameter + results table
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.set_facecolor(PAN); ax2.axis('off')

    iz_L   = int(np.argmin(np.abs(z_arr - P['emitter_length'])))
    all_E  = [pf['Em'][iz_L]*1e-6 for pf in profiles.values()]
    E_cor_ = [pf['Em'][iz_L]*1e-6 for pf in profiles.values() if pf['cls']=='corner']
    E_ctr_ = [pf['Em'][iz_L]*1e-6 for pf in profiles.values() if pf['cls']=='centre']
    shld   = (max(all_E)-min(all_E))/max(all_E)*100 if len(all_E)>1 else 0

    Ne = P['n_x']*P['n_y']
    rows = [
        ("Parameter",              "Value"),
        ("Array",                  f"{P['n_x']}×{P['n_y']} = {Ne} emitters"),
        ("Pitch x/y",              f"{P['pitch_x']*1e3:.2f}/{P['pitch_y']*1e3:.2f} mm"),
        ("Emitter length",         f"{P['emitter_length']*1e3:.2f} mm"),
        ("Tip half-angle",         f"{P['tip_half_angle_deg']:.0f}°"),
        ("Tip radius",             f"{P['tip_radius']*1e6:.0f} µm"),
        ("V_emitter / V_extr.",    f"{P['V_emitter']:.0f} / {P['V_extractor']:.0f} V"),
        ("Gap (tip→extr.)",        f"{P['gap']*1e3:.2f} mm"),
        ("pitch / gap",            f"{P['pitch_x']/P['gap']:.2f}"),
        ("Flow/emitter",           f"{P['flow_per_emitter']*1e12:.0f} pL/s"),
        ("Grid N",                 f"{P['N_grid']}"),
        ("SOR ω",                  f"{P['omega']:.2f}"),
        ("Converged at iter",      f"{conv_it}"),
        ("─────────────────",      "─────────────"),
        ("Numerical E_tip (corner)", f"{np.mean(E_cor_):.3f} MV/m"),
        ("Numerical E_tip (centre)", f"{np.mean(E_ctr_):.3f} MV/m" if E_ctr_ else "N/A"),
        ("Analytical E_tip",       f"{E_tip_corr*1e-6:.3f} MV/m"),
        ("Shielding",              f"{shld:.2f}%"),
        ("Corner/centre ratio",    f"{np.mean(E_cor_)/np.mean(E_ctr_):.3f}" if E_ctr_ else "N/A"),
    ]
    tbl = ax2.table(cellText=[r for r in rows[1:]],
                    colLabels=rows[0],
                    cellLoc='left', loc='center', colWidths=[0.58, 0.42])
    tbl.auto_set_font_size(False); tbl.set_fontsize(8)
    for (r,c), cell in tbl.get_celld().items():
        cell.set_facecolor(PAN if r%2==0 else '#1a2130')
        cell.set_edgecolor(BDR)
        cell.set_text_props(color=TXT if r>0 else ACC)
    ax2.set_title("Parameters & Results", color=TXT, fontsize=10, pad=6)

    plt.savefig(f'{OUT}/10_convergence_summary.png', dpi=150,
                bbox_inches='tight', facecolor=BG)
    plt.close(fig)
    print("  ✓ Fig 10 saved")


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def run(P):
    global z_arr, profiles_glob

    os.makedirs(OUT, exist_ok=True)

    print("\n" + "="*72)
    print("  ELECTROSPRAY MULTI-EMITTER 3-D E-FIELD SOLVER")
    print(f"  Array: {P['n_x']}×{P['n_y']}  |  V={P['V_emitter']:.0f}V  |  "
          f"gap={P['gap']*1e3:.1f}mm  |  pitch={P['pitch_x']*1e3:.1f}mm")
    print("="*72)

    # 1. Grid
    print("\n── 1. GRID ──────────────────────────────────────────────────────────")
    x, y, z, dx, dy, dz = build_grid(P)
    z_arr = z   # expose globally for plotting

    # 2. Geometry
    print("\n── 2. GEOMETRY ──────────────────────────────────────────────────────")
    (X3,Y3,Z3,body,extr,centres,cls_map,
     cx_u,cy_u,iz_tip,iz_ext,iz_base) = build_geometry(x,y,z,P)

    # 3. Initial condition
    print("\n── 3. INITIAL CONDITION ─────────────────────────────────────────────")
    L_ref = P['emitter_length'] + P['gap']
    phi0  = P['V_emitter'] * np.maximum(0.0, 1.0 - z[np.newaxis,np.newaxis,:]/L_ref)
    phi0  = np.broadcast_to(phi0, (len(x),len(y),len(z))).copy()
    phi0[body] = P['V_emitter']
    phi0[extr] = P['V_extractor']
    fixed = body | extr
    print(f"  Fixed nodes: {fixed.sum():,} / {fixed.size:,}")

    # 4. Solve
    print("\n── 4. SOR SOLVER ────────────────────────────────────────────────────")
    phi, hist = sor_solve(phi0, fixed, dx, dy, dz, P, label='main')

    # 5. Electric field
    print("\n── 5. ELECTRIC FIELD ────────────────────────────────────────────────")
    Ex, Ey, Ez_, Em = compute_efield(phi, dx, dy, dz)
    print(f"  Peak |E| in domain: {Em.max()*1e-6:.4f} MV/m")

    # 6. Tip correction
    E_tip_corr = tip_field_corrected(P)
    print(f"  Analytical E_tip (hyperboloid): {E_tip_corr*1e-6:.4f} MV/m")

    # 7. Profiles
    print("\n── 6. AXIAL PROFILES ────────────────────────────────────────────────")
    profiles = extract_profiles(phi,Ex,Ey,Ez_,Em,x,y,z,centres,cls_map,P)
    profiles_glob = profiles  # expose for fig08

    # Coupling summary
    iz_L   = int(np.argmin(np.abs(z - P['emitter_length'])))
    all_E  = [pf['Em'][iz_L] for pf in profiles.values()]
    shld   = (max(all_E)-min(all_E))/max(all_E)*100 if len(all_E)>1 else 0
    E_cor_ = [pf['Em'][iz_L]*1e-6 for pf in profiles.values() if pf['cls']=='corner']
    E_ctr_ = [pf['Em'][iz_L]*1e-6 for pf in profiles.values() if pf['cls']=='centre']

    print(f"\n  ┌{'─'*52}┐")
    print(f"  │  COUPLING RESULTS                                  │")
    print(f"  │  Shielding (field non-uniformity): {shld:6.2f}%        │")
    print(f"  │  Corner emitters |E|_tip: {np.mean(E_cor_):.4f} MV/m        │")
    if E_ctr_:
        print(f"  │  Centre emitter  |E|_tip: {np.mean(E_ctr_):.4f} MV/m        │")
        print(f"  │  Corner/centre ratio:     {np.mean(E_cor_)/np.mean(E_ctr_):.4f}            │")
    print(f"  │  Analytical tip E:        {E_tip_corr*1e-6:.4f} MV/m        │")
    print(f"  └{'─'*52}┘")

    for (cx,cy), pf in sorted(profiles.items()):
        print(f"  [{pf['cls']:6s}]  {pf['label']:22s}  "
              f"|E|_tip = {pf['Em'][iz_L]*1e-6:.4f} MV/m")

    # 8. Parametric sweep
    print("\n── 7. PARAMETRIC COUPLING SWEEP ─────────────────────────────────────")
    sweep = parametric_sweep(P, P['sweep_ratios'])

    # 9. Plots
    print("\n── 8. GENERATING FIGURES ────────────────────────────────────────────")
    fig01_axial(profiles, centres, cls_map, cx_u, cy_u, E_tip_corr, P)
    fig02_xy_slices(phi,Ex,Ey,Ez_,Em,x,y,z,centres,P,iz_tip,iz_ext)
    fig03_xz(phi,Ex,Ey,Ez_,Em,x,y,z,centres,P)
    fig04_coupling(profiles,centres,cls_map,cx_u,cy_u,P)
    fig05_isosurface_stack(phi,Em,x,y,z,centres,P)
    fig06_parametric(sweep, P)
    fig07_gradient(profiles,centres,cls_map,phi,Ez_,x,y,z,P)
    fig08_radial(phi,Ex,Ey,Em,x,y,z,centres,cls_map,P)
    fig09_equipotential(phi,Ex,Ez_,x,y,z,centres,P)
    fig10_convergence(hist, P, profiles, E_tip_corr)

    print("\n── DONE ─────────────────────────────────────────────────────────────")
    print(f"  10 figures saved to {OUT}/")
    return phi, Ex, Ey, Ez_, Em, profiles, hist, sweep


# ══════════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    result = run(P)