from nutils import mesh, solver, export, function, cli, topology
from nutils.expression_v2 import Namespace
import numpy as np
import json
import os
import precice

# --- Load Configuration ---
with open('../config.json', 'r') as f:
    config = json.load(f)
fluid_params = config['fluid']
sim_params = config['simulation']

# --- preCICE setup ---
participant_name = "Fluid"
config_file_name = "../precice-config.xml"
interface = precice.Participant(participant_name, config_file_name, 0, 1)
mesh_name = "Fluid-Mesh"
read_data_name = "Temperature"
write_data_name = "Heat-Flux"

# --- Folder name for results ---
timestamp = np.datetime64('now').astype(str).replace(':', '-').replace(' ', '_')
results_folder = f"{sim_params['results_folder_prefix']}_fluid_{timestamp}"
if not os.path.exists(results_folder):
    os.makedirs(results_folder)
                                                            
# --- Nutils setup ---
print("\n--- DEBUG: Attempting to load mesh with nutils.mesh.gmsh... ---")
domain, geom = mesh.gmsh("U_Channel_Fluid_mesh.msh")
print("--- DEBUG: Mesh loaded successfully. ---")

ns = Namespace()
ns.x = geom
ns.define_for('x', gradient='d', normal='n', jacobians=('dV', 'dS'))

# --- Fluid & Thermal Properties (from config.json) ---
ns.rho = fluid_params['density']
ns.mu = fluid_params['dynamic_viscosity']
ns.k = fluid_params['thermal_conductivity']
ns.cp = fluid_params['specific_heat']
ns.Tinlet = fluid_params['inlet_temp'] # New inlet temperature parameter
ns.velocity = np.array([
    fluid_params['inlet_velocity_x'],
    fluid_params['inlet_velocity_y'],
    fluid_params['inlet_velocity_z']
])

# --- Dynamic Boundary Assignment for Reversible Flow ---
# Determine flow direction (assuming y-axis is dominant)
if fluid_params['inlet_velocity_y'] >= 0:
    inlet_boundary_name = 'U_Channel_Fluid_Entry'
    outlet_boundary_name = 'U_Channel_Fluid_Outlet'
else:
    # If flow is reversed, swap the boundary definitions
    inlet_boundary_name = 'U_Channel_Fluid_Outlet'
    outlet_boundary_name = 'U_Channel_Fluid_Entry'
print(f"--- DEBUG: Flow direction detected. Inlet: '{inlet_boundary_name}', Outlet: '{outlet_boundary_name}' ---")

# --- Field Definitions ---
# Create the field objects and store them in dedicated Python variables.
temp_field = domain.field('t', btype='std', degree=2)
test_field = domain.field('v', btype='std', degree=2)

# Add the fields to the namespace so the expression parser can find them by name.
ns.t = temp_field
ns.v = test_field

# --- Variational Formulation (Non-Linear Style for 3D Mesh) ---
try:
    print("--- DEBUG: Setting up non-linear style variational formulation. ---")
    
    # PDE term using 'v' as the test function. No .var('t') needed.
    res = domain.integral('(rho cp velocity_j d_j(t) v + k d_j(t) d_j(v)) dV' @ ns, degree=4)

    # Outlet boundary condition using 'v' and the dynamic outlet name
    res += domain.boundary[outlet_boundary_name].integral('-k d_j(t) n_j v dS' @ ns, degree=4)

     # --- Boundary Conditions (Non-Linear Style) ---
    # Inlet: Fixed temperature, enforced by squaring the difference
    print("--- DEBUG: Setting up Dirichlet BC for non-linear solver. ---")
    # Use the dynamic inlet name and the new temperature parameter
    sqr = domain.boundary[inlet_boundary_name].integral('(t - Tinlet)^2 dS' @ ns, degree=4)
    # Use the older, explicit System object to generate constraints
    cons = solver.System(sqr, trial='t').solve_constraints(droptol=1e-15)

    # --- preCICE Coupling Setup for separate boundaries ---
    coupling_boundaries = ["U_Channel_Fluid_Channel_Top", "U_Channel_Fluid_Channel_Right", "U_Channel_Fluid_Channel_Bottom"]
    bases = {name: domain.boundary[name].basis('std', degree=2) for name in coupling_boundaries}
    
    # Add the bases to the namespace so the expression parser can find them
    ns.basesTop = bases["U_Channel_Fluid_Channel_Top"]
    ns.basesRight = bases["U_Channel_Fluid_Channel_Right"]
    ns.basesBottom = bases["U_Channel_Fluid_Channel_Bottom"]
    
    # Define arguments for the temperatures received from preCICE
    ns.ptTop = function.Argument('ptTop', shape=(len(ns.basesTop),))
    ns.ptRight = function.Argument('ptRight', shape=(len(ns.basesRight),))
    ns.ptBottom = function.Argument('ptBottom', shape=(len(ns.basesBottom),))

    # Pre-compute the interpolated temperature fields. This simplifies the expression for the parser.
    ns.interpTempTop = ns.ptTop @ ns.basesTop
    ns.interpTempRight = ns.ptRight @ ns.basesRight
    ns.interpTempBottom = ns.ptBottom @ ns.basesBottom

    # Add coupling conditions to the residual using 'v'
    res += domain.boundary["U_Channel_Fluid_Channel_Top"].integral('v (t - interpTempTop) dS' @ ns, degree=4)
    res += domain.boundary["U_Channel_Fluid_Channel_Right"].integral('v (t - interpTempRight) dS' @ ns, degree=4)
    res += domain.boundary["U_Channel_Fluid_Channel_Bottom"].integral('v (t - interpTempBottom) dS' @ ns, degree=4)
    print("--- DEBUG: Variational formulation and coupling setup complete. ---")

except KeyError as e:
    print(f"\n\nFATAL ERROR: A required boundary is missing from the mesh file: {e}")
    exit(1)
    
# --- Time Stepping Loop ---
print("--- DEBUG: Initializing solution vector lhs ---")
lhs = cons.copy()

# Get vertex IDs and split points for preCICE data mapping
print("--- DEBUG: Projecting boundary coordinates for preCICE mesh vertices ---")
all_coords = []
for name in coupling_boundaries:
    print(f"--- DEBUG: Projecting coordinates for boundary: {name} ---")
    boundary_basis = bases[name]
    # Project each component of the geometry 'x' onto the basis to find the node coordinates.
    x_coords = domain.boundary[name].project(ns.x[0], onto=boundary_basis, geometry=ns.x, ischeme='gauss4')
    y_coords = domain.boundary[name].project(ns.x[1], onto=boundary_basis, geometry=ns.x, ischeme='gauss4')
    z_coords = domain.boundary[name].project(ns.x[2], onto=boundary_basis, geometry=ns.x, ischeme='gauss4')
    # Stack the components into an (N, 3) array.
    coords = np.stack([x_coords, y_coords, z_coords], axis=-1)
    print(f"--- DEBUG: Projected {coords.shape[0]} coordinates for {name} ---")
    all_coords.append(coords)
all_coords_np = np.concatenate(all_coords)
print(f"--- DEBUG: Total projected coordinates for all coupling boundaries: {all_coords_np.shape} ---")

# Set all vertices for the mesh in a single call
print(f"--- DEBUG: Registering {all_coords_np.shape[0]} mesh vertices with preCICE ---")
vertex_ids = interface.set_mesh_vertices(mesh_name, all_coords_np)

# Data structures for splitting the preCICE data array
split_indices = [len(bases[name]) for name in coupling_boundaries]
split_points = np.cumsum(split_indices)[:-1]
print(f"--- DEBUG: Split indices for boundaries: {split_indices} ---")
print(f"--- DEBUG: Split points for np.split: {split_points} ---")

# Initialize preCICE before the main loop
print("--- DEBUG: Initializing preCICE interface ---")
interface.initialize()
print("--- DEBUG: preCICE interface initialized. ---")

time_step_number = 0
# Pre-allocate the array for reading data
num_vertices = len(vertex_ids)
precice_temp_all = np.zeros(num_vertices)

while interface.is_coupling_ongoing():
    print("--- DEBUG: Starting new preCICE coupling iteration ---")
    # Get the time step size for the current window
    dt = interface.get_max_time_step_size()
    print(f"--- DEBUG: preCICE max time step size: {dt}")

    # Verify for a checkpoint action
    if interface.requires_writing_checkpoint():
        print("--- DEBUG: preCICE action write iteration checkpoint required ---")
        saved_lhs = lhs.copy()
    
    # Read temperature data from preCICE
    precice_temp_all = interface.read_data(mesh_name, read_data_name, vertex_ids, 0.0 * dt)
    temps_top, temps_right, temps_bottom = np.split(precice_temp_all, split_points)
    print(f"--- DEBUG: Received temperature data shape: {precice_temp_all.shape}")
    
    # Split the array into parts for each boundary
    temps_top, temps_right, temps_bottom = np.split(precice_temp_all, split_points)
    print(f"--- DEBUG: Split temperature data: top={temps_top.shape}, right={temps_right.shape}, bottom={temps_bottom.shape}")
    
    # --- Update the temperature fields in the namespace
    print("--- DEBUG: Solving non-linear system with updated boundary temperatures ---")
    # Set the runtime arguments from preCICE
    arguments = dict(ptTop=temps_top, ptRight=temps_right, ptBottom=temps_bottom)
    # Add the initial guess for the temperature field 't' to the arguments.
    # The key for the solution vector is an empty string ''.
    arguments[''] = lhs
    # The solve method returns a dictionary: {'t': array([...])}
    solution_dict = solver.System(res, trial='t', test='v').solve(constrain=cons, arguments=arguments)
    # We must extract the numerical array from the dictionary.
    lhs = solution_dict['t']
    print("--- DEBUG: Linear system solved ---")
    
    # --- FLUX CALCULATION (String-Based Method) ---
    # 1. Define the normal heat flux as a string expression.
    #    Nutils uses spaces for multiplication. 'd_i(t)' is the gradient.    
    flux_expr = '-k d_i(t) n_i' @ ns
    
    # 2. Create the arguments dictionary to substitute the solution 'lhs' for the field 't'.
    #    The key must be the name of the field, 't'.
    eval_args = {'t': lhs}

    # 3. Project the expression. The 'arguments' parameter tells project how to evaluate 't'.
    flux_top = domain.boundary["U_Channel_Fluid_Channel_Top"].project(flux_expr, onto=bases["U_Channel_Fluid_Channel_Top"], geometry=ns.x, ischeme='gauss4', arguments=eval_args)
    flux_right = domain.boundary["U_Channel_Fluid_Channel_Right"].project(flux_expr, onto=bases["U_Channel_Fluid_Channel_Right"], geometry=ns.x, ischeme='gauss4', arguments=eval_args)
    flux_bottom = domain.boundary["U_Channel_Fluid_Channel_Bottom"].project(flux_expr, onto=bases["U_Channel_Fluid_Channel_Bottom"], geometry=ns.x, ischeme='gauss4', arguments=eval_args)
    print(f"--- DEBUG: Heat flux shapes: top={flux_top.shape}, right={flux_right.shape}, bottom={flux_bottom.shape}")
    
    # Concatenate fluxes and write back to preCICE
    heat_flux_all = np.concatenate([flux_top, flux_right, flux_bottom])
    interface.write_data(mesh_name, write_data_name, vertex_ids, heat_flux_all)
    print(f"--- DEBUG: Writing heat flux data to preCICE, shape: {heat_flux_all.shape} ---")
    
    print("--- DEBUG: Advancing preCICE interface ---")
    interface.advance(dt)
    
    # Verify if reading a checkpoint is required
    if interface.requires_reading_checkpoint():
        print("--- DEBUG: preCICE action read iteration checkpoint required ---")
        # Restore the saved lhs if a checkpoint was written
        lhs = saved_lhs.copy()
    
    # --- VTK Export ---
    print("--- DEBUG: Exporting fluid solution to VTK ---")
    vtk_sample = domain.sample('vtk', 1.0)
    points = vtk_sample.eval(ns.x)
    # Evaluate the temperature and heat flux vector using the same argument dictionary.
    temperature_values = vtk_sample.eval(ns.t, arguments=eval_args)
    heat_flux_values = vtk_sample.eval('-k d_i(t)' @ ns, arguments=eval_args)
    filename = f"fluid_solution_{time_step_number}"
    full_path = f"{results_folder}/{filename}"
    time_step_number += 1
    export.vtk(full_path, vtk_sample.tri, points, T=temperature_values, H=heat_flux_values)


# --- Finalize ---
print("--- DEBUG: Finalizing preCICE interface ---")
interface.finalize()
print("Fluid solver finished.")