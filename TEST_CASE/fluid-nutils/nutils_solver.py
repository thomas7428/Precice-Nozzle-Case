from nutils import mesh, solver, export, function
from nutils.expression_v2 import Namespace
import numpy as np

class NutilsSolver:
    """
    A class to encapsulate the Nutils-specific parts of the fluid simulation.
    """
    def __init__(self, fluid_params, inlet_boundary_name, outlet_boundary_name):
        """
        Initializes the Nutils domain, namespace, and variational formulation.
        """
        print("\n--- [NutilsSolver] Initializing ---")
        # --- Nutils setup ---
        self.domain, self.geom = mesh.gmsh("U_Channel_Fluid_mesh.msh")
        self.ns = Namespace()
        self.ns.x = self.geom
        self.ns.define_for('x', gradient='d', normal='n', jacobians=('dV', 'dS'))

        # --- Fluid & Thermal Properties ---
        self.ns.rho = fluid_params['density']
        self.ns.mu = fluid_params['dynamic_viscosity']
        self.ns.k = fluid_params['thermal_conductivity']
        self.ns.cp = fluid_params['specific_heat']
        self.ns.Tinlet = fluid_params['inlet_temp']
        self.ns.velocity = np.array([
            fluid_params['inlet_velocity_x'],
            fluid_params['inlet_velocity_y'],
            fluid_params['inlet_velocity_z']
        ])

        # --- Field Definitions ---
        self.ns.t = self.domain.field('t', btype='std', degree=2)
        self.ns.v = self.domain.field('v', btype='std', degree=2)

        # --- Variational Formulation ---
        self.res = self.domain.integral('(rho cp velocity_j d_j(t) v + k d_j(t) d_j(v)) dV' @ self.ns, degree=4)
        self.res += self.domain.boundary[outlet_boundary_name].integral('-k d_j(t) n_j v dS' @ self.ns, degree=4)

        # --- Boundary Conditions ---
        sqr = self.domain.boundary[inlet_boundary_name].integral('(t - Tinlet)^2 dS' @ self.ns, degree=4)
        self.cons = solver.System(sqr, trial='t').solve_constraints(droptol=1e-15)

        # --- preCICE Coupling Setup ---
        self.coupling_boundaries = ["U_Channel_Fluid_Channel_Top", "U_Channel_Fluid_Channel_Right", "U_Channel_Fluid_Channel_Bottom"]
        self.bases = {name: self.domain.boundary[name].basis('std', degree=2) for name in self.coupling_boundaries}
        
        # Define arguments and interpolated fields for coupling
        for name in self.coupling_boundaries:
            # Define unique names for the argument, basis, and interpolated field
            arg_name = f'pt{name.split("_")[-1]}'
            basis_name = f'basis{name.split("_")[-1]}'
            interp_name = f'interpTemp{name.split("_")[-1]}'
            
            # Add the function.Argument and basis to the namespace
            setattr(self.ns, arg_name, function.Argument(arg_name, shape=(len(self.bases[name]),)))
            setattr(self.ns, basis_name, self.bases[name])
            
            # Pre-calculate the interpolated temperature function and add it to the namespace
            interp_temp_func = getattr(self.ns, arg_name) @ getattr(self.ns, basis_name)
            setattr(self.ns, interp_name, interp_temp_func)
            
            # Add the coupling boundary condition to the residual using the simplified name
            self.res += self.domain.boundary[name].integral(f'v (t - {interp_name}) dS' @ self.ns, degree=4)
        
        print("--- [NutilsSolver] Initialization complete ---")

    def get_coupling_mesh_coords(self):
        """
        Projects and returns the 3D coordinates of all coupling boundary nodes.
        """
        all_coords = []
        for name in self.coupling_boundaries:
            boundary_basis = self.bases[name]
            coords = np.stack([
                self.domain.boundary[name].project(self.ns.x[i], onto=boundary_basis, geometry=self.ns.x, ischeme='gauss4')
                for i in range(3)
            ], axis=-1)
            all_coords.append(coords)
        return np.concatenate(all_coords)
    
    def solve_step(self, lhs, temp_data_split):
        """
        Solves the system for one time step with the given boundary temperatures.
        
        :param lhs: The current solution vector (initial guess). Can be a dict on first run.
        :param temp_data_split: A dictionary mapping boundary names to temperature arrays.
        :return: The new solution vector as a numpy array.
        """
        print("\n--- DEBUG: Entering solve_step ---")
        print(f"Type of incoming 'lhs': {type(lhs)}")

        # This logic handles both the initial dictionary and subsequent vector objects.
        if isinstance(lhs, dict):
            print("--- DEBUG: 'lhs' is a dictionary (initial run). Using as is. ---")
            initial_guess = lhs
        else:
            print(f"--- DEBUG: 'lhs' is a vector object of shape {lhs.shape}. Sanitizing to numpy array. ---")
            # This is the crucial sanitization step for iterations > 0.
            # It prevents the 'unhashable type' error.
            initial_guess = np.array(lhs)
        
        print(f"Type of 'initial_guess' being passed to solver: {type(initial_guess)}")
        print("---")

        # --- Build the arguments dictionary ---
        # The initial guess for the primary unknown ('t') is passed with an empty string key.
        # This is the correct pattern for this solver setup.
        arguments = {f'pt{name.split("_")[-1]}': temps for name, temps in temp_data_split.items()}
        arguments[''] = initial_guess
        
        # --- Final check before calling the solver ---
        print("--- DEBUG: Final pre-solve state ---")
        print(f"Solver will be called with initial guess of type: {type(arguments[''])}")
        print(f"Constraints 'self.cons' are of type: {type(self.cons)}")
        # CORRECTED: Dictionaries do not have a .shape attribute. Use len() to get the number of constraints.
        print(f"Number of constraints in 'self.cons': {len(self.cons)}")
        print("------------------------------------")

        # --- Solve the system ---
        solution_dict = solver.System(self.res, trial='t', test='v').solve(constrain=self.cons, arguments=arguments)
        
        # --- Post-solve information ---
        print("--- DEBUG: System solved successfully ---")
        print(f"Type of returned solution: {type(solution_dict['t'])}")
        
        # Return the solution vector for the next iteration
        return solution_dict['t']

    def _get_complete_eval_args(self, lhs):
        """Helper function to build a complete arguments dictionary for evaluation."""
        # Ensure lhs is a numpy array, not a dict.
        if isinstance(lhs, dict):
            lhs = lhs.get('t', np.zeros_like(next(iter(lhs.values()))))

        eval_args = {'t': lhs}
        for name in self.coupling_boundaries:
            arg_name = f'pt{name.split("_")[-1]}'
            eval_args[arg_name] = np.zeros(len(self.bases[name]))
        return eval_args
    
    def calculate_flux(self, lhs):
        """
        Calculates and returns the heat flux on all coupling boundaries.
        """
        eval_args = self._get_complete_eval_args(lhs)
        flux_expr = self.ns.k * self.ns.d(self.ns.t) @ self.ns.n
        all_fluxes = []
        for name in self.coupling_boundaries:
            flux = self.domain.boundary[name].project(flux_expr, onto=self.bases[name], geometry=self.ns.x, ischeme='gauss4', arguments=eval_args)
            all_fluxes.append(flux)
        return np.concatenate(all_fluxes)
    
    def export_vtk(self, lhs, results_folder, time_step_number):
        """
        Exports the current solution to a VTK file.
        """
        print(f"--- DEBUG: Exporting fluid solution to VTK for step {time_step_number} ---")
        vtk_sample = self.domain.sample('vtk', degree=2)
        eval_args = self._get_complete_eval_args(lhs)

        points = vtk_sample.eval(self.ns.x)
        temperature = vtk_sample.eval(self.ns.t, arguments=eval_args)
        heat_flux_vector = -self.ns.k * self.ns.d(self.ns.t)
        heat_flux = vtk_sample.eval(heat_flux_vector, arguments=eval_args)
        
        tri = vtk_sample.tri
        if tri.shape[1] == 4:
            tri = tri[:, :3]
            
        full_path = f"{results_folder}/fluid_solution_{time_step_number}"
        try:
            export.vtk(full_path, tri, points, T=temperature, H=heat_flux)
        except Exception as e:
            print(f"--- ERROR: Failed to export VTK file: {e} ---")
