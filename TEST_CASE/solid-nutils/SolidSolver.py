import numpy as np
import json
import os
import precice
from solid_solver import SolidSolver
import time

# --- Load Configuration ---
with open('../config.json', 'r') as f:
    config = json.load(f)
solid_params = config['solid']
sim_params = config['simulation']
export_vtk = solid_params.get('export_vtk', False)

# --- preCICE setup ---
participant_name = "Solid"
config_file_name = "../precice-config.xml"
interface = precice.Participant(participant_name, config_file_name, 0, 1)
mesh_name = "Solid-Mesh"
read_data_name = "Heat-Flux"
write_data_name = "Temperature"

# --- Folder name for results ---
timestamp = time.strftime("%M_%H_%d_%m")
results_folder = f"../{sim_params['results_folder_prefix']}_{timestamp}"
if not os.path.exists(results_folder):
    os.makedirs(results_folder)

# --- Instantiate the Nutils Solver ---
try:
    solid_solver = SolidSolver(solid_params)
except KeyError as e:
    print(f"\n\nFATAL ERROR: A required boundary is missing from the mesh file: {e}")
    exit(1)

# --- preCICE Mesh and Data Mapping Setup ---
all_coords_np = solid_solver.get_coupling_mesh_coords()
vertex_ids = interface.set_mesh_vertices(mesh_name, all_coords_np)

split_indices = [len(solid_solver.bases[name]) for name in solid_solver.coupling_boundaries]
split_points = np.cumsum(split_indices)[:-1]

# --- Time Stepping Loop ---
print(f"--- DEBUG: Initializing solid temperature field to ambient: {solid_solver.ambient_temp} ---")
lhs = np.full(len(solid_solver.ns.basis), solid_solver.ambient_temp)
time_step_number = 0
saved_lhs = lhs.copy()

interface.initialize()
print("--- DEBUG: preCICE interface initialized. ---")

while interface.is_coupling_ongoing():

    print(f"\n--- Coupling Iteration: {time_step_number} ---\n")

    dt = interface.get_max_time_step_size()

    if interface.requires_writing_checkpoint():
        print("--- DEBUG: preCICE action: write iteration checkpoint ---")
        saved_lhs = lhs.copy()
    
    # CORRECTED: The read_data function requires the timestep 'dt' as a fourth argument.
    precice_flux_all = interface.read_data(mesh_name, read_data_name, vertex_ids, dt)
    
    flux_arrays = np.split(precice_flux_all, split_points)
    flux_data_split = dict(zip(solid_solver.coupling_boundaries, flux_arrays))
    
    # Solve the system
    solution = solid_solver.solve_step(lhs, flux_data_split)
    
    # CRITICAL FIX: Ensure the solution is a numerical numpy array of floats.
    # The solver sometimes returns symbolic Array<> objects, which must be converted.
    lhs = np.asarray(solution, dtype=float)

    print("--- DEBUG: Nutils system solved ---")
    
    interface_temps_all = solid_solver.get_temperatures_for_precice(lhs)

    interface.write_data(mesh_name, write_data_name, vertex_ids, interface_temps_all)
    print(f"--- DEBUG: Wrote heat flux data to preCICE, shape: {precice_flux_all.shape} ---")

    interface.advance(dt)
    
    if interface.requires_reading_checkpoint():
        print("--- DEBUG: preCICE action: read iteration checkpoint ---")
        lhs = saved_lhs.copy()
    
    if export_vtk:
        solid_solver.export_vtk(lhs, results_folder, time_step_number)

    time_step_number += 1
    
interface.finalize()
print("Solid solver finished.")
