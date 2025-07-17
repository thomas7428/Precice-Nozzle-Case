from nutils import mesh, solver, export, function
from nutils.expression_v2 import Namespace
import numpy as np

class SolidSolver:
    """
    A class to encapsulate the Nutils-specific parts of the solid simulation.
    """
    def __init__(self, solid_params):
        """
        Initializes the Nutils domain, namespace, and variational formulation.
        """
        print("\n--- [SolidSolver] Initializing ---")
        self.domain, self.geom = mesh.gmsh("U_Channel_Solid_mesh.msh")
        self.ns = Namespace()
        self.ns.x = self.geom
        self.ns.define_for('x', gradient='d', normal='n', jacobians=('dV', 'dS'))
        self.ns.basis = self.domain.basis('std', degree=2)
        self.ns.t = self.ns.basis @ function.Argument('lhs', shape=(len(self.ns.basis),))
        self.ns.k = solid_params['thermal_conductivity']
        self.ns.Qhot = solid_params['hot_gas_flux']
        self.ns.h = solid_params['convection_coeff']
        self.ns.Tambient = solid_params['ambient_temp']
        # Store Tambient for access by the main script
        self.ambient_temp = self.ns.Tambient

        # --- Variational Formulation ---
        self.res = self.domain.integral('k d_j(basis_i) d_j(t) dV' @ self.ns, degree=4)
        self.res -= self.domain.boundary["U_Channel_Solid_Bottom"].integral('basis_i Qhot dS' @ self.ns, degree=4)
        self.res += self.domain.boundary["U_Channel_Solid_Top"].integral('basis_i h (t - Tambient) dS' @ self.ns, degree=4)
        self.cons = None

        # --- preCICE Coupling Setup ---
        self.coupling_boundaries = [
            "U_Channel_Solid_Channel_Top",
            "U_Channel_Solid_Channel_Right",
            "U_Channel_Solid_Channel_Bottom"
        ]
        self.bases = {name: self.domain.boundary[name].basis('std', degree=2) for name in self.coupling_boundaries}

        # Define a preCICE heat flux argument for each boundary
        self.ns.pqTop = function.Argument('pqTop', shape=(len(self.bases["U_Channel_Solid_Channel_Top"]),))
        self.ns.pqRight = function.Argument('pqRight', shape=(len(self.bases["U_Channel_Solid_Channel_Right"]),))
        self.ns.pqBottom = function.Argument('pqBottom', shape=(len(self.bases["U_Channel_Solid_Channel_Bottom"]),))

        # Pre-compute the interpolated heat flux fields.
        self.ns.interpFluxTop = self.bases["U_Channel_Solid_Channel_Top"] @ self.ns.pqTop
        self.ns.interpFluxRight = self.bases["U_Channel_Solid_Channel_Right"] @ self.ns.pqRight
        self.ns.interpFluxBottom = self.bases["U_Channel_Solid_Channel_Bottom"] @ self.ns.pqBottom

        # Add coupling terms to the residual
        self.res -= self.domain.boundary["U_Channel_Solid_Channel_Top"].integral('basis_i interpFluxTop dS' @ self.ns, degree=4)
        self.res -= self.domain.boundary["U_Channel_Solid_Channel_Right"].integral('basis_i interpFluxRight dS' @ self.ns, degree=4)
        self.res -= self.domain.boundary["U_Channel_Solid_Channel_Bottom"].integral('basis_i interpFluxBottom dS' @ self.ns, degree=4)
        
        print("--- [SolidSolver] Initialization complete ---")

    def get_coupling_mesh_coords(self):
        """Projects and returns the 3D coordinates of all coupling boundary nodes."""
        all_coords = []
        for name in self.coupling_boundaries:
            boundary_basis = self.bases[name]
            # CORRECTED: Project each coordinate component (x, y, z) individually
            # as the project function expects a scalar field.
            x_coords = self.domain.boundary[name].project(self.ns.x[0], onto=boundary_basis, geometry=self.ns.x, ischeme='gauss4')
            y_coords = self.domain.boundary[name].project(self.ns.x[1], onto=boundary_basis, geometry=self.ns.x, ischeme='gauss4')
            z_coords = self.domain.boundary[name].project(self.ns.x[2], onto=boundary_basis, geometry=self.ns.x, ischeme='gauss4')
            # Stack the individual coordinate arrays to form the final (N, 3) array.
            coords = np.stack([x_coords, y_coords, z_coords], axis=-1)
            all_coords.append(coords)
        return np.concatenate(all_coords)

    def solve_step(self, lhs, flux_data_split):
        """Solves the linear system for one time step."""
        arguments = {
            'pqTop': flux_data_split["U_Channel_Solid_Channel_Top"],
            'pqRight': flux_data_split["U_Channel_Solid_Channel_Right"],
            'pqBottom': flux_data_split["U_Channel_Solid_Channel_Bottom"]
        }
        # --- Final check before calling the solver ---
        print("--- DEBUG: Final pre-solve state ---")
        print(f"Constraints 'self.cons' are of type: {type(self.cons)}")
        print("------------------------------------")
        return solver.solve_linear('lhs', self.res, constrain=self.cons, arguments=arguments)

    def get_temperatures_for_precice(self, lhs):
        """Calculates the temperatures on the coupling boundaries to be sent to preCICE."""
        all_temps = []
        for name in self.coupling_boundaries:
            temps = self.domain.boundary[name].project(self.ns.t, onto=self.bases[name], geometry=self.ns.x, ischeme='gauss4', arguments=dict(lhs=lhs))
            all_temps.append(temps)
        return np.concatenate(all_temps)

    def export_vtk(self, lhs, results_folder, time_step_number):
        """Exports the current solution to a VTK file."""
        print(f"--- DEBUG: Exporting solid solution to VTK for step {time_step_number} ---")
        vtk_sample = self.domain.sample('vtk', degree=2)
        points = vtk_sample.eval(self.ns.x)
        print(f"--- DEBUG: points.shape = {points.shape}")
        print(f"--- DEBUG: lhs type: {type(lhs)}, shape: {getattr(lhs, 'shape', None)}")
        print(f"--- DEBUG: Number of basis functions: {len(self.ns.basis)}")
        # Print a small sample of lhs
        print(f"--- DEBUG: lhs sample: {lhs[:10] if hasattr(lhs, '__getitem__') else lhs}")
        try:
            temperature_values = vtk_sample.eval(self.ns.t, arguments={'lhs': lhs})
            print(f"--- DEBUG: temperature_values.shape = {temperature_values.shape}")
            heat_flux_vector = vtk_sample.eval('-k d_i(t)' @ self.ns, arguments={'lhs': lhs})
            print(f"--- DEBUG: heat_flux_vector.shape = {heat_flux_vector.shape}")
            full_path = f"{results_folder}/solid_solution_{time_step_number}"
            export.vtk(full_path, vtk_sample.tri, points, T=temperature_values, H=heat_flux_vector)
        except Exception as e:
            print(f"Failed to export VTK: {e}")
            import traceback
            traceback.print_exc()