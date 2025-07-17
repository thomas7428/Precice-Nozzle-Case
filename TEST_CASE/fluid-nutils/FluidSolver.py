from nutils import solver
import numpy as np
import json
import os
import precice
import time
from nutils_solver import NutilsSolver # Import the new solver class

# --- Load Configuration ---
with open('../config.json', 'r') as f:
    config = json.load(f)
fluid_params = config['fluid']
sim_params = config['simulation']
export_vtk = fluid_params['export_vtk']

# --- preCICE setup ---
participant_name = "Fluid"
config_file_name = "../precice-config.xml"
interface = precice.Participant(participant_name, config_file_name, 0, 1)
mesh_name = "Fluid-Mesh"
read_data_name = "Temperature"
write_data_name = "Heat-Flux"

# --- Folder name for results ---
timestamp = time.strftime("%M_%H_%d_%m")
results_folder = f"../{sim_params['results_folder_prefix']}_{timestamp}"
if not os.path.exists(results_folder):
    os.makedirs(results_folder)

# --- Dynamic Boundary Assignment for Reversible Flow ---
if fluid_params['inlet_velocity_y'] >= 0:
    inlet_boundary_name = 'U_Channel_Fluid_Entry'
    outlet_boundary_name = 'U_Channel_Fluid_Outlet'
else:
    inlet_boundary_name = 'U_Channel_Fluid_Outlet'
    outlet_boundary_name = 'U_Channel_Fluid_Entry'
print(f"--- DEBUG: Flow direction detected. Inlet: '{inlet_boundary_name}', Outlet: '{outlet_boundary_name}' ---")

# --- Instantiate the Nutils Solver ---
try:
    nutils_solver = NutilsSolver(fluid_params, inlet_boundary_name, outlet_boundary_name)
except KeyError as e:
    print(f"\n\nFATAL ERROR: A required boundary is missing from the mesh file: {e}")
    exit(1)

# --- preCICE Mesh and Data Mapping Setup ---
print("--- DEBUG: Projecting boundary coordinates for preCICE mesh vertices ---")
all_coords_np = nutils_solver.get_coupling_mesh_coords()
vertex_ids = interface.set_mesh_vertices(mesh_name, all_coords_np)
print(f"--- DEBUG: Registered {len(vertex_ids)} mesh vertices with preCICE ---")

# Create data structures for splitting the preCICE data array
split_indices = [len(nutils_solver.bases[name]) for name in nutils_solver.coupling_boundaries]
split_points = np.cumsum(split_indices)[:-1]

# --- Time Stepping Loop ---
lhs = nutils_solver.cons.copy()
time_step_number = 0

interface.initialize()
print("--- DEBUG: preCICE interface initialized. ---")

while interface.is_coupling_ongoing():
    print(f"\n--- Coupling Iteration: {time_step_number} ---")
    dt = interface.get_max_time_step_size()

    if interface.requires_writing_checkpoint():
        print("--- DEBUG: preCICE action: write iteration checkpoint ---")
        saved_lhs = lhs.copy()
    
    # Read temperature data from preCICE
    precice_temp_all = interface.read_data(mesh_name, read_data_name, vertex_ids, 0.0 * dt)
    
    # Split the flat array into a list of arrays for each boundary
    temp_arrays = np.split(precice_temp_all, split_points)
    # Create a dictionary mapping boundary names to their temperature data
    temp_data_split = dict(zip(nutils_solver.coupling_boundaries, temp_arrays))
    
    # Solve the system for the current time step
    lhs = nutils_solver.solve_step(lhs, temp_data_split)
    print("--- DEBUG: Nutils system solved ---")
    
    # Calculate heat flux on the coupling boundaries
    heat_flux_all = nutils_solver.calculate_flux(lhs)
    
    # Write heat flux data to preCICE
    interface.write_data(mesh_name, write_data_name, vertex_ids, heat_flux_all)
    print(f"--- DEBUG: Wrote heat flux data to preCICE, shape: {heat_flux_all.shape} ---")
    
    interface.advance(dt)
    
    if interface.requires_reading_checkpoint():
        print("--- DEBUG: preCICE action: read iteration checkpoint ---")
        lhs = saved_lhs.copy()
    
    if export_vtk:
        nutils_solver.export_vtk(lhs, results_folder, time_step_number)

    time_step_number += 1
    
# --- Finalize ---
print("--- DEBUG: Finalizing preCICE interface ---")
interface.finalize()
print("Fluid solver finished.")