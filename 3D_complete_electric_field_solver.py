#!/usr/bin/env python3
"""
Streamlined 3D Electric Field Solver for Multi-Emitter Arrays
==============================================================
Essential features only: Configure, Solve, Plot

Controls:
1. n_emitters_x  - Number of emitters in x-direction
2. n_emitters_y  - Number of emitters in y-direction  
3. spacing_x     - Spacing in x-direction (m)
4. spacing_y     - Spacing in y-direction (m)

Author: Multi-Emitter System
Date: 2026-02-12
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import time


class ElectricFieldSolver3D:
    """
    3D Electric Field Solver
    
    Parameters:
    -----------
    n_emitters_x : int
        Number of emitters in x-direction
    n_emitters_y : int
        Number of emitters in y-direction
    spacing_x : float
        Spacing between emitters in x (meters)
    spacing_y : float
        Spacing between emitters in y (meters)
    """
    
    def __init__(self, 
                 n_emitters_x=5, 
                 n_emitters_y=5,
                 spacing_x=400e-6,
                 spacing_y=400e-6,
                 V_emitter=0.0,
                 V_extractor=-2000.0,
                 z_extractor=1000e-6):
        
        # CONTROL PARAMETERS
        self.n_emitters_x = n_emitters_x
        self.n_emitters_y = n_emitters_y
        self.spacing_x = spacing_x
        self.spacing_y = spacing_y
        
        # Voltages
        self.V_emitter = V_emitter
        self.V_extractor = V_extractor
        self.z_extractor = z_extractor
        
        # Emitter geometry
        self.r_tip = 10e-6
        self.r_hole = 60e-6
        
        # Generate emitter positions
        self._generate_emitter_positions()
        
        # Create mesh
        self._create_mesh()
        
        print(f"\n{'='*70}")
        print(f"3D ELECTRIC FIELD SOLVER INITIALIZED")
        print(f"{'='*70}")
        print(f"  Array: {self.n_emitters_x} × {self.n_emitters_y} = {self.n_emitters} emitters")
        print(f"  Spacing: Δx={self.spacing_x*1e6:.0f} μm, Δy={self.spacing_y*1e6:.0f} μm")
        print(f"  Array size: {self.array_width*1e3:.2f} × {self.array_height*1e3:.2f} mm")
        print(f"  Mesh: {self.nx} × {self.ny} × {self.nz} = {self.nx*self.ny*self.nz:,} points")
    
    def _generate_emitter_positions(self):
        """Generate emitter positions in square grid"""
        x_pos = np.arange(self.n_emitters_x) * self.spacing_x
        y_pos = np.arange(self.n_emitters_y) * self.spacing_y
        
        # Center at origin
        x_pos -= np.mean(x_pos)
        y_pos -= np.mean(y_pos)
        
        X, Y = np.meshgrid(x_pos, y_pos)
        self.emitter_positions = np.column_stack([X.ravel(), Y.ravel()])
        self.n_emitters = len(self.emitter_positions)
        
        self.array_width = np.ptp(x_pos)
        self.array_height = np.ptp(y_pos)
    
    def _create_mesh(self):
        """Create computational mesh"""
        # Domain size (auto-sized to array with margin)
        margin = 1e-3
        x_size = self.array_width + 2*margin
        y_size = self.array_height + 2*margin
        z_size = 2.5e-3
        
        # Grid resolution (adaptive to array size)
        self.nx = 70
        self.ny = 70
        self.nz = 90
        
        # Create grid
        self.x = np.linspace(-x_size/2, x_size/2, self.nx)
        self.y = np.linspace(-y_size/2, y_size/2, self.ny)
        self.z = np.linspace(0, z_size, self.nz)
        
        self.dx = self.x[1] - self.x[0]
        self.dy = self.y[1] - self.y[0]
        self.dz = self.z[1] - self.z[0]
    
    def _idx(self, i, j, k):
        """Convert 3D indices to 1D"""
        return i * (self.ny * self.nz) + j * self.nz + k
    
    def _is_on_emitter(self, x, y, z):
        """Check if point is on emitter surface"""
        tolerance = max(self.dx, self.dy, self.dz) * 1.5
        
        for px, py in self.emitter_positions:
            r_dist = np.sqrt((x - px)**2 + (y - py)**2)
            
            if z < 100e-6:  # Cone height
                r_cone = z * np.tan(49.3 * np.pi / 180)
                if abs(r_dist - r_cone) < tolerance and r_dist < 3*self.r_tip:
                    return True
        return False
    
    def _is_on_extractor(self, x, y, z):
        """Check if point is on extractor plate"""
        if abs(z - self.z_extractor) > 1.5*self.dz:
            return False
        
        # Check if in hole
        for px, py in self.emitter_positions:
            if np.sqrt((x - px)**2 + (y - py)**2) < self.r_hole:
                return False
        
        return True
    
    def solve(self):
        """Solve Laplace equation"""
        print(f"\n{'='*70}")
        print("SOLVING 3D LAPLACE EQUATION")
        print(f"{'='*70}")
        
        t_start = time.time()
        
        # Build system
        print("\n  Building linear system...")
        N = self.nx * self.ny * self.nz
        A = sp.lil_matrix((N, N))
        b = np.zeros(N)
        
        dx2, dy2, dz2 = self.dx**2, self.dy**2, self.dz**2
        
        for i in range(self.nx):
            if i % 10 == 0:
                print(f"    Progress: {i/self.nx*100:.0f}%", end='\r')
            
            for j in range(self.ny):
                for k in range(self.nz):
                    n = self._idx(i, j, k)
                    x_val, y_val, z_val = self.x[i], self.y[j], self.z[k]
                    
                    # Boundary: Emitter
                    if self._is_on_emitter(x_val, y_val, z_val):
                        A[n, n] = 1.0
                        b[n] = self.V_emitter
                    
                    # Boundary: Extractor
                    elif self._is_on_extractor(x_val, y_val, z_val):
                        A[n, n] = 1.0
                        b[n] = self.V_extractor
                    
                    # Boundary: Domain edges (Neumann)
                    elif (i == 0 or i == self.nx-1 or j == 0 or 
                          j == self.ny-1 or k == self.nz-1):
                        A[n, n] = 1.0
                        if i == 0:
                            A[n, self._idx(1, j, k)] = -1.0
                        elif i == self.nx-1:
                            A[n, self._idx(self.nx-2, j, k)] = -1.0
                        elif j == 0:
                            A[n, self._idx(i, 1, k)] = -1.0
                        elif j == self.ny-1:
                            A[n, self._idx(i, self.ny-2, k)] = -1.0
                        elif k == self.nz-1:
                            A[n, self._idx(i, j, self.nz-2)] = -1.0
                    
                    # Interior: Laplace equation
                    else:
                        coeff = -2.0 * (1/dx2 + 1/dy2 + 1/dz2)
                        A[n, n] = coeff
                        A[n, self._idx(i-1, j, k)] = 1.0/dx2
                        A[n, self._idx(i+1, j, k)] = 1.0/dx2
                        A[n, self._idx(i, j-1, k)] = 1.0/dy2
                        A[n, self._idx(i, j+1, k)] = 1.0/dy2
                        A[n, self._idx(i, j, k-1)] = 1.0/dz2
                        A[n, self._idx(i, j, k+1)] = 1.0/dz2
        
        print(f"\n  System built: {N:,} equations")
        
        # Solve
        print("\n  Solving linear system (BiCGSTAB)...")
        A_csr = A.tocsr()
        
        # Preconditioner
        M = sp.diags(1.0 / A_csr.diagonal())
        
        # Solve with BiCGSTAB (different parameter names in different scipy versions)
        try:
            phi_vec, info = spla.bicgstab(A_csr, b, M=M, atol=1e-5, maxiter=5000)
        except TypeError:
            # Older scipy version uses 'tol' instead of 'atol'
            phi_vec, info = spla.bicgstab(A_csr, b, M=M, x0=None, tol=1e-5, maxiter=5000)
        except:
            # Fallback to basic call
            phi_vec, info = spla.bicgstab(A_csr, b, M=M, maxiter=5000)
        
        if info == 0:
            print("  ✓ Solution converged")
        else:
            print(f"  ⚠ Convergence issue (info={info})")
        
        # Reshape
        self.phi = phi_vec.reshape((self.nx, self.ny, self.nz))
        
        # Compute electric field
        print("\n  Computing electric field...")
        self._compute_electric_field()
        
        t_total = time.time() - t_start
        
        print(f"\n{'='*70}")
        print(f"  SOLUTION COMPLETE")
        print(f"  Potential range: [{np.min(self.phi):.1f}, {np.max(self.phi):.1f}] V")
        print(f"  Max field: {np.max(self.E_mag)/1e6:.2f} MV/m")
        print(f"  Total time: {t_total:.1f} seconds")
        print(f"{'='*70}")
    
    def _compute_electric_field(self):
        """Compute E = -∇φ"""
        self.Ex = np.zeros_like(self.phi)
        self.Ey = np.zeros_like(self.phi)
        self.Ez = np.zeros_like(self.phi)
        
        # Central differences
        self.Ex[1:-1, :, :] = -(self.phi[2:, :, :] - self.phi[:-2, :, :])/(2*self.dx)
        self.Ey[:, 1:-1, :] = -(self.phi[:, 2:, :] - self.phi[:, :-2, :])/(2*self.dy)
        self.Ez[:, :, 1:-1] = -(self.phi[:, :, 2:] - self.phi[:, :, :-2])/(2*self.dz)
        
        # Boundaries
        self.Ex[0, :, :] = -(self.phi[1, :, :] - self.phi[0, :, :])/self.dx
        self.Ex[-1, :, :] = -(self.phi[-1, :, :] - self.phi[-2, :, :])/self.dx
        self.Ey[:, 0, :] = -(self.phi[:, 1, :] - self.phi[:, 0, :])/self.dy
        self.Ey[:, -1, :] = -(self.phi[:, -1, :] - self.phi[:, -2, :])/self.dy
        self.Ez[:, :, 0] = -(self.phi[:, :, 1] - self.phi[:, :, 0])/self.dz
        self.Ez[:, :, -1] = -(self.phi[:, :, -1] - self.phi[:, :, -2])/self.dz
        
        self.E_mag = np.sqrt(self.Ex**2 + self.Ey**2 + self.Ez**2)
    
    def plot_all(self):
        """Generate all plots with error handling"""
        print(f"\n{'='*70}")
        print("GENERATING PLOTS")
        print(f"{'='*70}")
        
        # 1. 2D slices
        try:
            self._plot_2d_slices()
        except Exception as e:
            print(f"    ⚠ Warning: Could not generate 2D slices: {e}")
        
        # 2. 3D field magnitude
        try:
            self._plot_3d_field_surface()
        except Exception as e:
            print(f"    ⚠ Warning: Could not generate 3D surface: {e}")
        
        # 3. 3D field vectors
        try:
            self._plot_3d_field_vectors()
        except Exception as e:
            print(f"    ⚠ Warning: Could not generate 3D vectors: {e}")
        
        # 4. On-axis field
        try:
            self._plot_on_axis_field()
        except Exception as e:
            print(f"    ⚠ Warning: Could not generate on-axis plot: {e}")
        
        print(f"\n{'='*70}")
        print("PLOTTING COMPLETE")
        print(f"{'='*70}")
    
    def _plot_2d_slices(self):
        """Plot 2D slices (XY, XZ, YZ)"""
        print("\n  Creating 2D slice plots...")
        
        fig = plt.figure(figsize=(18, 5))
        
        # XY slice at emitter tips
        ax1 = fig.add_subplot(131)
        k = 5  # Near emitter tips
        X, Y = np.meshgrid(self.x*1e3, self.y*1e3, indexing='ij')
        
        cs1 = ax1.contourf(X, Y, self.E_mag[:, :, k]/1e6, levels=25, cmap='hot')
        ax1.contour(X, Y, self.phi[:, :, k], levels=15, colors='white', 
                    linewidths=0.5, alpha=0.3)
        
        # Mark emitters
        for px, py in self.emitter_positions:
            ax1.plot(px*1e3, py*1e3, 'go', markersize=8)
        
        plt.colorbar(cs1, ax=ax1, label='|E| (MV/m)')
        ax1.set_xlabel('x (mm)')
        ax1.set_ylabel('y (mm)')
        ax1.set_title(f'XY Slice at z={self.z[k]*1e3:.2f} mm')
        ax1.set_aspect('equal')
        
        # XZ slice through center
        ax2 = fig.add_subplot(132)
        j = self.ny // 2
        X, Z = np.meshgrid(self.x*1e3, self.z*1e3, indexing='ij')
        
        cs2 = ax2.contourf(X, Z, self.E_mag[:, j, :]/1e6, levels=25, cmap='hot')
        ax2.contour(X, Z, self.phi[:, j, :], levels=15, colors='white',
                    linewidths=0.5, alpha=0.3)
        
        ax2.axhline(self.z_extractor*1e3, color='cyan', linestyle='--', 
                    linewidth=2, label='Extractor')
        plt.colorbar(cs2, ax=ax2, label='|E| (MV/m)')
        ax2.set_xlabel('x (mm)')
        ax2.set_ylabel('z (mm)')
        ax2.set_title('XZ Slice (center)')
        ax2.legend()
        
        # YZ slice through center
        ax3 = fig.add_subplot(133)
        i = self.nx // 2
        Y, Z = np.meshgrid(self.y*1e3, self.z*1e3, indexing='ij')
        
        cs3 = ax3.contourf(Y, Z, self.E_mag[i, :, :]/1e6, levels=25, cmap='hot')
        ax3.contour(Y, Z, self.phi[i, :, :], levels=15, colors='white',
                    linewidths=0.5, alpha=0.3)
        
        ax3.axhline(self.z_extractor*1e3, color='cyan', linestyle='--',
                    linewidth=2, label='Extractor')
        plt.colorbar(cs3, ax=ax3, label='|E| (MV/m)')
        ax3.set_xlabel('y (mm)')
        ax3.set_ylabel('z (mm)')
        ax3.set_title('YZ Slice (center)')
        ax3.legend()
        
        plt.tight_layout()
        fig.savefig('field_2d_slices.png', dpi=200, bbox_inches='tight')
        print("    ✓ Saved: field_2d_slices.png")
        plt.close()
    
    def _plot_3d_field_surface(self):
        """Plot 3D field magnitude surface"""
        print("\n  Creating 3D field surface plot...")
        
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Slice at mid-height
        k = self.nz // 3
        X, Y = np.meshgrid(self.x*1e3, self.y*1e3, indexing='ij')
        E_slice = self.E_mag[:, :, k]
        
        surf = ax.plot_surface(X, Y, E_slice/1e6, cmap='hot', 
                              alpha=0.8, vmin=0, 
                              vmax=np.percentile(E_slice, 95)/1e6)
        
        # Mark emitters
        for px, py in self.emitter_positions:
            ax.plot([px*1e3], [py*1e3], [0], 'go', markersize=6)
        
        ax.set_xlabel('x (mm)')
        ax.set_ylabel('y (mm)')
        ax.set_zlabel('|E| (MV/m)')
        ax.set_title(f'3D Electric Field Magnitude (z={self.z[k]*1e3:.1f} mm)')
        
        fig.colorbar(surf, ax=ax, shrink=0.5, label='|E| (MV/m)')
        
        fig.savefig('field_3d_surface.png', dpi=200, bbox_inches='tight')
        print("    ✓ Saved: field_3d_surface.png")
        plt.close()
    
    def _plot_3d_field_vectors(self):
        """Plot 3D field vector visualization"""
        print("\n  Creating 3D field vector plot...")
        
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Subsample for quiver
        skip = 8
        X_sub = self.x[::skip]
        Y_sub = self.y[::skip]
        k = self.nz // 4  # Slice height
        Z_val = self.z[k]
        
        X_grid, Y_grid = np.meshgrid(X_sub, Y_sub, indexing='ij')
        Ex_sub = self.Ex[::skip, ::skip, k]
        Ey_sub = self.Ey[::skip, ::skip, k]
        Ez_sub = self.Ez[::skip, ::skip, k]
        E_mag_sub = self.E_mag[::skip, ::skip, k]
        
        # Normalize for visibility
        E_norm = np.sqrt(Ex_sub**2 + Ey_sub**2 + Ez_sub**2) + 1e-10
        Ex_norm = Ex_sub / E_norm
        Ey_norm = Ey_sub / E_norm
        Ez_norm = Ez_sub / E_norm
        
        # Flatten for quiver
        X_flat = X_grid.flatten()
        Y_flat = Y_grid.flatten()
        Z_flat = np.full_like(X_flat, Z_val)
        Ex_flat = Ex_norm.flatten()
        Ey_flat = Ey_norm.flatten()
        Ez_flat = Ez_norm.flatten()
        E_mag_flat = E_mag_sub.flatten()
        
        # Color by magnitude (single color value per arrow)
        colors = plt.cm.hot(E_mag_flat / (np.max(E_mag_flat) + 1e-10))
        
        # Plot arrows one by one with color
        for i in range(0, len(X_flat), 2):  # Skip some for clarity
            ax.quiver(X_flat[i]*1e3, Y_flat[i]*1e3, Z_flat[i]*1e3,
                     Ex_flat[i], Ey_flat[i], Ez_flat[i],
                     length=0.3, normalize=True,
                     color=colors[i], arrow_length_ratio=0.3)
        
        # Emitters
        for px, py in self.emitter_positions:
            ax.plot([px*1e3], [py*1e3], [0], 'go', markersize=8)
        
        # Extractor
        x_ext = np.linspace(self.x[0], self.x[-1], 20)
        y_ext = np.linspace(self.y[0], self.y[-1], 20)
        X_ext, Y_ext = np.meshgrid(x_ext, y_ext)
        Z_ext = np.ones_like(X_ext) * self.z_extractor
        
        ax.plot_surface(X_ext*1e3, Y_ext*1e3, Z_ext*1e3,
                       alpha=0.2, color='blue')
        
        ax.set_xlabel('x (mm)')
        ax.set_ylabel('y (mm)')
        ax.set_zlabel('z (mm)')
        ax.set_title('3D Electric Field Vectors')
        
        fig.savefig('field_3d_vectors.png', dpi=200, bbox_inches='tight')
        print("    ✓ Saved: field_3d_vectors.png")
        plt.close()
    
    def _plot_on_axis_field(self):
        """Plot field along emitter axes"""
        print("\n  Creating on-axis field plot...")
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # Center emitter
        center_idx = self.n_emitters // 2
        px, py = self.emitter_positions[center_idx]
        
        i = np.argmin(np.abs(self.x - px))
        j = np.argmin(np.abs(self.y - py))
        
        E_axis = self.E_mag[i, j, :]
        phi_axis = self.phi[i, j, :]
        
        ax1.plot(self.z*1e3, E_axis/1e6, 'b-', linewidth=2, label='Center emitter')
        
        # Corner emitter
        corner_idx = 0
        px_c, py_c = self.emitter_positions[corner_idx]
        i_c = np.argmin(np.abs(self.x - px_c))
        j_c = np.argmin(np.abs(self.y - py_c))
        E_corner = self.E_mag[i_c, j_c, :]
        
        ax1.plot(self.z*1e3, E_corner/1e6, 'r--', linewidth=2, label='Corner emitter')
        
        ax1.axvline(self.z_extractor*1e3, color='gray', linestyle='--', 
                   linewidth=2, alpha=0.5, label='Extractor')
        ax1.set_xlabel('z (mm)')
        ax1.set_ylabel('|E| (MV/m)')
        ax1.set_title('Electric Field Along Emitter Axis')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Potential
        ax2.plot(self.z*1e3, phi_axis, 'b-', linewidth=2, label='Center')
        ax2.plot(self.z*1e3, self.phi[i_c, j_c, :], 'r--', linewidth=2, label='Corner')
        ax2.axvline(self.z_extractor*1e3, color='gray', linestyle='--',
                   linewidth=2, alpha=0.5, label='Extractor')
        ax2.set_xlabel('z (mm)')
        ax2.set_ylabel('Potential (V)')
        ax2.set_title('Potential Along Emitter Axis')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        fig.savefig('field_on_axis.png', dpi=200, bbox_inches='tight')
        print("    ✓ Saved: field_on_axis.png")
        plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# MAIN EXECUTION
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    
    print("\n" + "="*80)
    print(" "*20 + "3D ELECTRIC FIELD SOLVER")
    print(" "*25 + "Streamlined Version")
    print("="*80)
    
    # ========================================================================
    # CONFIGURE: Set your 3 control parameters here
    # ========================================================================
    
    n_emitters_x = 5       # Number in x-direction
    n_emitters_y = 5       # Number in y-direction
    spacing_x = 400e-6     # Spacing in x (400 μm)
    spacing_y = 400e-6     # Spacing in y (400 μm)
    
    # ========================================================================
    # CREATE SOLVER
    # ========================================================================
    
    solver = ElectricFieldSolver3D(
        n_emitters_x=n_emitters_x,
        n_emitters_y=n_emitters_y,
        spacing_x=spacing_x,
        spacing_y=spacing_y,
        V_emitter=0.0,          # Ground
        V_extractor=-2000.0,    # -2 kV
        z_extractor=1000e-6     # 1 mm gap
    )
    
    # ========================================================================
    # SOLVE
    # ========================================================================
    
    solver.solve()
    
    # ========================================================================
    # PLOT ALL
    # ========================================================================
    
    solver.plot_all()
    
    # ========================================================================
    # SUMMARY
    # ========================================================================
    
    print("\n" + "="*80)
    print("COMPLETE!")
    print("="*80)
    print("\nGenerated files:")
    print("  1. field_2d_slices.png   - 2D slices (XY, XZ, YZ)")
    print("  2. field_3d_surface.png  - 3D field magnitude surface")
    print("  3. field_3d_vectors.png  - 3D field vectors")
    print("  4. field_on_axis.png     - Field along emitter axes")
    print("\n" + "="*80)