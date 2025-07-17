from nutils import mesh, solver, export, function, cli, topology
from nutils.expression_v2 import Namespace
import numpy as np
import json
import os
import precice

# --- Load Configuration ---
with open('../config.json', 'r') as f:
    config = json.load(f)
solid_params = config['solid']
sim_params = config['simulation']

# --- preCICE setup ---
participant_name = "Solid"
config_file_name = "../precice-config.xml"
interface = precice.Participant(participant_name, config_file_name, 0, 1)

mesh_name = "Solid-Mesh"
# NOTE: Read/Write data names are swapped compared to the fluid solver
read_data_name = "Heat-Flux"
write_data_name = "Temperature"

# --- Folder name for results ---
timestamp = np.datetime64('now').astype(str).replace(':', '-').replace(' ', '_')
results_folder = f"{sim_params['results_folder_prefix']}_solid_{timestamp}"
if not os.path.exists(results_folder):
    os.makedirs(results_folder)
    
# --- Nutils setup ---
print("\n--- DEBUG: Attempting to load mesh with nutils.mesh.gmsh... ---")
domain, geom = mesh.gmsh("U_Channel_Solid_mesh.msh")
print("--- DEBUG: Mesh loaded successfully. ---")
ns = Namespace()
ns.x = geom
ns.define_for('x', gradient='d', normal='n', jacobians=('dV', 'dS'))
ns.basis = domain.basis('std', degree=2)
ns.t = ns.basis @ function.Argument('lhs', shape=(len(ns.basis),))
ns.k = solid_params['thermal_conductivity']
ns.Qhot = solid_params['hot_gas_flux']
ns.h = solid_params['convection_coeff']
ns.Tambient = solid_params['ambient_temp']
print("--- DEBUG: Namespace initialized with parameters. ---")
print(f"Loaded parameters: k={ns.k}, Qhot={ns.Qhot}, h={ns.h}, Tambient={ns.Tambient}")

# --- Variational Formulation ---
try:
    print("--- DEBUG: Setting up variational formulation. ---")
    res = domain.integral('k d_j(basis_i) d_j(t) dV' @ ns, degree=4)

    # --- Boundary Conditions ---
    # 1. Hot Gas Contact (Neumann BC on U_Channel_Solid_Bottom)
    res -= domain.boundary["U_Channel_Solid_Bottom"].integral('basis_i Qhot dS' @ ns, degree=4)
    
    # 2. Exterior Air Contact (Robin BC on U_Channel_Solid_Top)
    res += domain.boundary["U_Channel_Solid_Top"].integral('basis_i h (t - Tambient) dS' @ ns, degree=4)

    # No constraints are needed if there are no Dirichlet BCs
    cons = None

    # --- preCICE Coupling Setup ---
    # 3. Cooling Channel Contact (Coupled via preCICE)
    coupling_boundaries = [
        "U_Channel_Solid_Channel_Top",
        "U_Channel_Solid_Channel_Right",
        "U_Channel_Solid_Channel_Bottom"
    ]

    bases = {name: domain.boundary[name].basis('std', degree=2) for name in coupling_boundaries}
    # Define a preCICE heat flux argument for each boundary
    ns.pqTop = bases["U_Channel_Solid_Channel_Top"] @ function.Argument('pqTop', shape=(len(bases["U_Channel_Solid_Channel_Top"]),))
    ns.pqRight = bases["U_Channel_Solid_Channel_Right"] @ function.Argument('pqRight', shape=(len(bases["U_Channel_Solid_Channel_Right"]),))
    ns.pqBottom = bases["U_Channel_Solid_Channel_Bottom"] @ function.Argument('pqBottom', shape=(len(bases["U_Channel_Solid_Channel_Bottom"]),))

    # Add coupling condition (Neumann-type) to residual for EACH boundary
    # The minus sign indicates flux entering the domain.
    res -= domain.boundary["U_Channel_Solid_Channel_Top"].integral('basis_i pqTop dS' @ ns, degree=4)
    res -= domain.boundary["U_Channel_Solid_Channel_Right"].integral('basis_i pqRight dS' @ ns, degree=4)
    res -= domain.boundary["U_Channel_Solid_Channel_Bottom"].integral('basis_i pqBottom dS' @ ns, degree=4)
    print("--- DEBUG: Variational formulation and coupling setup complete. ---")

except KeyError as e:
    print(f"\n\nFATAL ERROR: A required boundary is missing from the mesh file: {e}")
    exit(1)

# --- Time Stepping Loop ---
lhs = np.zeros(len(ns.basis))

print("--- DEBUG: Projecting coupling node coordinates for preCICE mesh registration. ---")
# Get coupling node coordinates using the same robust projection method
all_coords = []
for name in coupling_boundaries:
    print(f"--- DEBUG: Projecting coordinates for boundary '{name}'. ---")
    boundary_basis = bases[name]
    x_coords = domain.boundary[name].project(ns.x[0], onto=boundary_basis, geometry=ns.x, ischeme='gauss4')
    y_coords = domain.boundary[name].project(ns.x[1], onto=boundary_basis, geometry=ns.x, ischeme='gauss4')
    z_coords = domain.boundary[name].project(ns.x[2], onto=boundary_basis, geometry=ns.x, ischeme='gauss4')
    coords = np.stack([x_coords, y_coords, z_coords], axis=-1)
    print(f"--- DEBUG: {len(coords)} nodes projected for '{name}'. ---")
    all_coords.append(coords)
all_coords_np = np.concatenate(all_coords)
print(f"--- DEBUG: Total coupling nodes for preCICE: {len(all_coords_np)} ---")

print("--- DEBUG: Registering mesh vertices with preCICE. ---")
vertex_ids = interface.set_mesh_vertices(mesh_name, all_coords_np)
print("--- DEBUG: Mesh vertices registered. ---")

# Data structures for splitting the preCICE data array
split_indices = [len(bases[name]) for name in coupling_boundaries]
split_points = np.cumsum(split_indices)[:-1]
print(f"--- DEBUG: Split points for coupling boundaries: {split_points} ---")

print("--- DEBUG: Initializing preCICE interface. ---")
interface.initialize()
print("--- DEBUG: preCICE interface initialized. ---")

# Pre-allocate the array for reading data
num_vertices = len(vertex_ids)
precice_flux_all = np.zeros(num_vertices)

time_step_number = 0

while interface.is_coupling_ongoing():
    dt = interface.get_max_time_step_size()
    print(f"--- DEBUG: New time step with dt = {dt} ---")
    
    #Verify if a checkpoint action is required
    if interface.requires_writing_checkpoint():
        print("--- DEBUG: preCICE action write iteration checkpoint required ---")
        saved_lhs = lhs.copy()
    
    print("--- DEBUG: Reading heat flux data from preCICE. ---")
    relative_read_time = 0.0 * dt
    precice_flux_all = interface.read_data(mesh_name, read_data_name, vertex_ids, relative_read_time)
    print(f"--- DEBUG: Received heat flux data: {precice_flux_all} ---")
    
    # Split the array into parts for each boundary
    flux_top, flux_right, flux_bottom = np.split(precice_flux_all, split_points)
    print(f"--- DEBUG: Split heat flux into boundaries: top({len(flux_top)}), right({len(flux_right)}), bottom({len(flux_bottom)}) ---")
    
    # Solve the system
    arguments = dict(pqTop=flux_top, pqRight=flux_right, pqBottom=flux_bottom)
    print("--- DEBUG: Solving linear system for temperature. ---")
    lhs = solver.solve_linear('lhs', res, constrain=cons, arguments=arguments)
    print(f"--- DEBUG: Solution vector (lhs) computed. ---")
    
    # Evaluate temperature at the coupling nodes and write back to preCICE
    all_temps = []
    for name in coupling_boundaries:
        boundary_basis = bases[name]
        print(f"--- DEBUG: Projecting temperature on boundary '{name}'. ---")
        temps = domain.boundary[name].project(ns.t, onto=boundary_basis, geometry=ns.x, ischeme='gauss4', arguments=dict(lhs=lhs))
        print(f"--- DEBUG: Projected temperatures for '{name}': {temps} ---")
        all_temps.append(temps)
    
    interface_temps_all = np.concatenate(all_temps)
    print(f"--- DEBUG: Writing temperature data to preCICE: {interface_temps_all} ---")
    interface.write_data(mesh_name, write_data_name, vertex_ids, interface_temps_all)
    
    print("--- DEBUG: Advancing preCICE interface. ---")
    interface.advance(dt)
    
    # Verify if reading a checkpoint is required
    if interface.requires_reading_checkpoint():
        print("--- DEBUG: preCICE action read iteration checkpoint required ---")
        # Restore the saved lhs if a checkpoint was written
        lhs = saved_lhs.copy()

    # --- VTK Export ---
    print("--- DEBUG: Exporting solid solution to VTK. ---")
    # Create a sample for VTK export and pass it to the exporter
    vtk_sample = domain.sample('vtk', 1.0)
    temperature_values = vtk_sample.eval(ns.t, lhs=lhs)
    heat_flux_values = vtk_sample.eval(ns.k * ns.d(ns.t), lhs=lhs)
    filename = f"solid_solution_{time_step_number}"
    full_path = f"{results_folder}/{filename}"
    points = vtk_sample.eval(ns.x)
    time_step_number += 1
    export.vtk(full_path, vtk_sample.tri, points, T=temperature_values, H=heat_flux_values)

# --- Finalize ---
print("--- DEBUG: Finalizing preCICE interface. ---")
interface.finalize()
print("Solid solver finished.")