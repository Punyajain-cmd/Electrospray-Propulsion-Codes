import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from dataclasses import dataclass
from typing import Tuple, Optional

@dataclass
class ElectricFieldSolver2D:
    """
    2D Axisymmetric Electric Field Solver using Finite Element Method
    
    Solves Laplace's equation in cylindrical coordinates (r, z):
    1/r ∂/∂r(r ∂φ/∂r) + ∂²φ/∂z² = 0
    
    Boundary Conditions:
    - Emitter tip: φ = 0 (ground)
    - Extraction electrode: φ = V_ext
    - Far field: ∂φ/∂n = 0 (Neumann)
    """
    
    # Geometry parameters (meters)
    r_emitter: float = 10e-6      # Emitter tip radius
    z_emitter: float = 0.0         # Emitter tip z-position
    r_electrode: float = 500e-6    # Electrode radius
    z_electrode: float = 500e-6    # Electrode distance from tip
    
    # Applied voltage
    V_extraction: float = 2000.0   # Volts
    
    # Mesh parameters
    nr: int = 100                  # Number of radial points
    nz: int = 150                  # Number of axial points
    r_max: float = 1e-3            # Maximum radial extent
    z_max: float = 2e-3            # Maximum axial extent
    
    # Computed fields
    phi: Optional[np.ndarray] = None      # Electric potential (V)
    E_r: Optional[np.ndarray] = None      # Radial electric field (V/m)
    E_z: Optional[np.ndarray] = None      # Axial electric field (V/m)
    E_mag: Optional[np.ndarray] = None    # Field magnitude (V/m)
    
    # Mesh
    r: Optional[np.ndarray] = None
    z: Optional[np.ndarray] = None
    
    def __post_init__(self):
        """Initialize mesh after dataclass creation"""
        self.create_mesh()
    
    def create_mesh(self):
        """Create 2D structured mesh"""
        self.r = np.linspace(0, self.r_max, self.nr)
        self.z = np.linspace(self.z_emitter, self.z_max, self.nz)
        self.R, self.Z = np.meshgrid(self.r, self.z, indexing='ij')
        
        # Grid spacing
        self.dr = self.r[1] - self.r[0]
        self.dz = self.z[1] - self.z[0]
    
    def apply_boundary_conditions(self) -> Tuple[sp.lil_matrix, np.ndarray]:
        """
        Apply boundary conditions using finite difference method
        Returns system matrix A and RHS vector b for Ax = b
        """
        N = self.nr * self.nz
        A = sp.lil_matrix((N, N))
        b = np.zeros(N)
        
        def idx(i, j):
            """Convert 2D indices to 1D index"""
            return i * self.nz + j
        
        # Finite difference coefficients
        dr2 = self.dr ** 2
        dz2 = self.dz ** 2
        
        for i in range(self.nr):
            for j in range(self.nz):
                n = idx(i, j)
                r_val = self.r[i]
                
                # Boundary conditions
                # 1. Emitter tip (Dirichlet: φ = 0)
                if self.is_on_emitter(i, j):
                    A[n, n] = 1.0
                    b[n] = 0.0
                
                # 2. Extraction electrode (Dirichlet: φ = V_ext)
                elif self.is_on_electrode(i, j):
                    A[n, n] = 1.0
                    b[n] = self.V_extraction
                
                # 3. Axis boundary (r=0): Use L'Hospital's rule
                elif i == 0:
                    # At r=0: 2∂²φ/∂r² + ∂²φ/∂z² = 0
                    A[n, n] = -2.0 * (2.0/dr2 + 1.0/dz2)
                    if i + 1 < self.nr:
                        A[n, idx(i+1, j)] = 4.0 / dr2
                    if j > 0:
                        A[n, idx(i, j-1)] = 1.0 / dz2
                    if j + 1 < self.nz:
                        A[n, idx(i, j+1)] = 1.0 / dz2
                
                # 4. Outer boundaries (Neumann: ∂φ/∂n = 0)
                elif i == self.nr - 1:  # r = r_max
                    A[n, n] = 1.0
                    A[n, idx(i-1, j)] = -1.0
                    b[n] = 0.0
                
                elif j == self.nz - 1:  # z = z_max
                    A[n, n] = 1.0
                    A[n, idx(i, j-1)] = -1.0
                    b[n] = 0.0
                
                # 5. Interior points: Solve Laplace equation in cylindrical coords
                else:
                    # Standard 5-point stencil for axisymmetric Laplacian
                    # 1/r ∂/∂r(r ∂φ/∂r) + ∂²φ/∂z² = 0
                    
                    coeff_center = -2.0 * (1.0/dr2 + 1.0/dz2)
                    coeff_r_plus = (1.0/dr2) * (1.0 + self.dr/(2.0*r_val))
                    coeff_r_minus = (1.0/dr2) * (1.0 - self.dr/(2.0*r_val))
                    coeff_z = 1.0 / dz2
                    
                    A[n, n] = coeff_center
                    
                    if i > 0:
                        A[n, idx(i-1, j)] = coeff_r_minus
                    if i + 1 < self.nr:
                        A[n, idx(i+1, j)] = coeff_r_plus
                    if j > 0:
                        A[n, idx(i, j-1)] = coeff_z
                    if j + 1 < self.nz:
                        A[n, idx(i, j+1)] = coeff_z
        
        return A.tocsr(), b
    
    def is_on_emitter(self, i: int, j: int) -> bool:
        """Check if point is on emitter surface (sharp tip at origin)"""
        r_val = self.r[i]
        z_val = self.z[j]
        
        # Model emitter as cone with tip at origin
        # Cone half-angle ~ 49° (Taylor cone)
        cone_angle = 49.3 * np.pi / 180
        
        if z_val < self.z_emitter + 10 * self.r_emitter:
            # Near tip: define emitter surface
            r_cone = (z_val - self.z_emitter) * np.tan(cone_angle)
            if abs(r_val - r_cone) < 1.5 * self.dr and r_val <= self.r_emitter * 3:
                return True
        
        return False
    
    def is_on_electrode(self, i: int, j: int) -> bool:
        """Check if point is on extraction electrode"""
        r_val = self.r[i]
        z_val = self.z[j]
        
        # Electrode is a ring at z = z_electrode
        if abs(z_val - self.z_electrode) < 1.5 * self.dz:
            if r_val >= self.r_emitter * 5 and r_val <= self.r_electrode:
                return True
        
        return False
    
    def solve(self, solver='spsolve', verbose=True):
        """
        Solve the Laplace equation for electric potential
        
        Parameters:
        -----------
        solver : str
            'spsolve' (direct) or 'cg' (conjugate gradient)
        verbose : bool
            Print solver information
        """
        if verbose:
            print("Setting up FEM system...")
        
        A, b = self.apply_boundary_conditions()
        
        if verbose:
            print(f"System size: {A.shape[0]} equations")
            print(f"Solving using {solver}...")
        
        # Solve linear system
        if solver == 'spsolve':
            phi_vec = spla.spsolve(A, b)
        elif solver == 'cg':
            phi_vec, info = spla.cg(A, b, tol=1e-8, maxiter=5000)
            if info != 0:
                print(f"Warning: CG did not converge (info={info})")
        else:
            raise ValueError(f"Unknown solver: {solver}")
        
        # Reshape to 2D grid
        self.phi = phi_vec.reshape((self.nr, self.nz))
        
        if verbose:
            print("Computing electric field...")
        
        # Compute electric field E = -∇φ
        self.E_r = np.zeros_like(self.phi)
        self.E_z = np.zeros_like(self.phi)
        
        # Central differences for interior points
        for i in range(1, self.nr - 1):
            for j in range(1, self.nz - 1):
                self.E_r[i, j] = -(self.phi[i+1, j] - self.phi[i-1, j]) / (2 * self.dr)
                self.E_z[i, j] = -(self.phi[i, j+1] - self.phi[i, j-1]) / (2 * self.dz)
        
        # Boundaries (forward/backward differences)
        for j in range(self.nz):
            self.E_r[0, j] = -(self.phi[1, j] - self.phi[0, j]) / self.dr
            self.E_r[-1, j] = -(self.phi[-1, j] - self.phi[-2, j]) / self.dr
        
        for i in range(self.nr):
            self.E_z[i, 0] = -(self.phi[i, 1] - self.phi[i, 0]) / self.dz
            self.E_z[i, -1] = -(self.phi[i, -1] - self.phi[i, -2]) / self.dz
        
        # Field magnitude
        self.E_mag = np.sqrt(self.E_r**2 + self.E_z**2)
        
        if verbose:
            print(f"Max potential: {np.max(self.phi):.2f} V")
            print(f"Max field: {np.max(self.E_mag)/1e6:.2f} MV/m")
            print("Solution complete!")
    
    def get_field_at_point(self, r_query: float, z_query: float) -> Tuple[float, float, float]:
        """
        Interpolate electric field at arbitrary point
        
        Returns: (E_r, E_z, E_mag) in V/m
        """
        if self.phi is None:
            raise RuntimeError("Must call solve() first")
        
        interp_Er = RegularGridInterpolator((self.r, self.z), self.E_r, 
                                             bounds_error=False, fill_value=0)
        interp_Ez = RegularGridInterpolator((self.r, self.z), self.E_z,
                                             bounds_error=False, fill_value=0)
        
        Er = float(interp_Er([r_query, z_query]))
        Ez = float(interp_Ez([r_query, z_query]))
        
        return Er, Ez, np.sqrt(Er**2 + Ez**2)
    
    def plot_solution(self, figsize=(14, 10)):
        """Create comprehensive visualization of the solution"""
        if self.phi is None:
            raise RuntimeError("Must call solve() first")
        
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
        
        # Convert to mm for plotting
        R_mm = self.R * 1e3
        Z_mm = self.Z * 1e3
        
        # 1. Electric Potential
        ax1 = fig.add_subplot(gs[0, 0])
        levels = np.linspace(0, self.V_extraction, 25)
        cs1 = ax1.contourf(Z_mm, R_mm, self.phi, levels=levels, cmap='viridis')
        ax1.contour(Z_mm, R_mm, self.phi, levels=10, colors='white', linewidths=0.5, alpha=0.3)
        plt.colorbar(cs1, ax=ax1, label='Potential (V)')
        ax1.set_xlabel('z (mm)')
        ax1.set_ylabel('r (mm)')
        ax1.set_title('Electric Potential φ(r, z)')
        ax1.set_xlim(0, self.z_max * 1e3)
        ax1.set_ylim(0, self.r_max * 1e3)
        
        # Mark emitter and electrode
        ax1.plot([0], [0], 'ro', markersize=8, label='Emitter tip')
        ax1.axhline(y=self.r_electrode * 1e3, color='white', ls='--', lw=2, alpha=0.7)
        ax1.legend()
        
        # 2. Electric Field Magnitude
        ax2 = fig.add_subplot(gs[0, 1])
        cs2 = ax2.contourf(Z_mm, R_mm, self.E_mag / 1e6, levels=25, cmap='hot')
        plt.colorbar(cs2, ax=ax2, label='Field (MV/m)')
        ax2.set_xlabel('z (mm)')
        ax2.set_ylabel('r (mm)')
        ax2.set_title('Electric Field Magnitude |E|')
        ax2.set_xlim(0, self.z_max * 1e3)
        ax2.set_ylim(0, self.r_max * 1e3)
        
        # 3. Electric Field Vectors
        ax3 = fig.add_subplot(gs[1, 0])
        # Subsample for quiver plot
        skip = (slice(None, None, 5), slice(None, None, 5))
        ax3.quiver(Z_mm[skip], R_mm[skip], 
                   self.E_z[skip], self.E_r[skip],
                   self.E_mag[skip] / 1e6,
                   cmap='plasma', scale=5e8, width=0.003)
        ax3.contour(Z_mm, R_mm, self.phi, levels=15, colors='gray', linewidths=0.5, alpha=0.5)
        ax3.set_xlabel('z (mm)')
        ax3.set_ylabel('r (mm)')
        ax3.set_title('Electric Field Vectors (color = magnitude)')
        ax3.set_xlim(0, self.z_max * 1e3)
        ax3.set_ylim(0, self.r_max * 1e3)
        
        # 4. Field on axis (r=0)
        ax4 = fig.add_subplot(gs[1, 1])
        z_axis_mm = self.z * 1e3
        E_axis = self.E_mag[0, :] / 1e6
        phi_axis = self.phi[0, :]
        
        ax4_twin = ax4.twinx()
        l1 = ax4.plot(z_axis_mm, E_axis, 'r-', lw=2, label='|E| on axis')
        l2 = ax4_twin.plot(z_axis_mm, phi_axis, 'b-', lw=2, label='φ on axis')
        
        ax4.set_xlabel('z (mm)')
        ax4.set_ylabel('Electric Field (MV/m)', color='r')
        ax4_twin.set_ylabel('Potential (V)', color='b')
        ax4.tick_params(axis='y', labelcolor='r')
        ax4_twin.tick_params(axis='y', labelcolor='b')
        ax4.set_title('On-Axis Electric Field & Potential')
        ax4.grid(True, alpha=0.3)
        
        # Combined legend
        lns = l1 + l2
        labs = [l.get_label() for l in lns]
        ax4.legend(lns, labs, loc='upper right')
        
        plt.suptitle(f'2D Electric Field Solution (V_ext = {self.V_extraction} V)', 
                     fontsize=14, fontweight='bold', y=0.995)
        
        return fig


# ============================================================
# Example usage and validation
# ============================================================
if __name__ == "__main__":
    print("="*60)
    print("2D Electric Field Solver for Electrospray Thrusters")
    print("="*60)
    
    # Create solver instance
    solver = ElectricFieldSolver2D(
        r_emitter=10e-6,          # 10 μm tip radius
        z_electrode=500e-6,       # 500 μm gap
        V_extraction=2000.0,      # 2 kV
        nr=120,                   # Fine mesh
        nz=180,
        r_max=1e-3,              # 1 mm radial extent
        z_max=2e-3               # 2 mm axial extent
    )
    
    # Solve the field
    solver.solve(solver='spsolve', verbose=True)
    
    # Query field at specific points
    print("\nField at key locations:")
    print("-" * 40)
    
    # At emitter tip
    E_r, E_z, E_mag = solver.get_field_at_point(0, 5e-6)
    print(f"At emitter (r=0, z=5 μm):")
    print(f"  E_z = {E_z/1e6:.2f} MV/m")
    print(f"  |E| = {E_mag/1e6:.2f} MV/m")
    
    # Midpoint
    E_r, E_z, E_mag = solver.get_field_at_point(50e-6, 250e-6)
    print(f"\nAt midpoint (r=50 μm, z=250 μm):")
    print(f"  E_r = {E_r/1e6:.2f} MV/m")
    print(f"  E_z = {E_z/1e6:.2f} MV/m")
    print(f"  |E| = {E_mag/1e6:.2f} MV/m")
    
    # Near electrode
    E_r, E_z, E_mag = solver.get_field_at_point(100e-6, 495e-6)
    print(f"\nNear electrode (r=100 μm, z=495 μm):")
    print(f"  E_r = {E_r/1e6:.2f} MV/m")
    print(f"  E_z = {E_z/1e6:.2f} MV/m")
    print(f"  |E| = {E_mag/1e6:.2f} MV/m")
    
    # Plot solution
    print("\nGenerating visualization...")
    fig = solver.plot_solution()
    fig.savefig('electric_field_solution_2d.png', dpi=150, bbox_inches='tight')
    print("Saved: electric_field_solution_2d.png")
    
    plt.close()
    
    print("\n" + "="*60)
    print("Electric field solver complete!")
    print("="*60)