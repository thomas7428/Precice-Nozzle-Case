
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
participant_name = "CombustionSolver"
config_file_name = "../precice-config.xml"
interface = precice.Participant(participant_name, config_file_name, 0, 1)
mesh_name = "Inner_Fluid_Mesh"
read_data_name = "Temperature"
write_data_name = "Heat-Flux"

# --- Folder name for results ---
timestamp = np.datetime64('now').astype(str).replace(':', '-').replace(' ', '_')
results_folder = f"{sim_params['results_folder_prefix']}_combustion_{timestamp}"
if not os.path.exists(results_folder):
    os.makedirs(results_folder)

# --- Nutils setup ---
print("\n--- DEBUG: Attempting to load mesh with nutils.mesh.gmsh... ---")
domain, geom = mesh.gmsh("Inner_Fluid_mesh.msh")
print("--- DEBUG: Mesh loaded successfully. ---")

ns = Namespace()
ns.x = geom
ns.define_for('x', gradient='d', normal='n', jacobians=('dV', 'dS'))

# --- Fluid & Combustion Properties (from config.json) ---
ns.rho = fluid_params['density']
ns.mu = fluid_params['dynamic_viscosity']
ns.k = fluid_params['thermal_conductivity']
ns.cp = fluid_params['specific_heat']
ns.Tinlet = fluid_params['inlet_temp']
ns.velocity = np.array([
    fluid_params['inlet_velocity_x'],
    fluid_params['inlet_velocity_y'],
    fluid_params['inlet_velocity_z']
])

# --- Boundary Names (adapt to your mesh) ---
inlet_boundary_name = 'Inlet'
outlet_boundary_name = 'Outlet'
combustion_wall_boundaries = ["Wall"]
coupling_boundaries = ["Wall"]  # For PreCICE coupling, adapt as needed


# --- Field Definitions: velocity (u,v), pressure (p), temperature (t), fuel mass fraction (y) ---
u_field = domain.field('u', btype='std', degree=2)
v_field = domain.field('v', btype='std', degree=2)
p_field = domain.field('p', btype='std', degree=1)
t_field = domain.field('t', btype='std', degree=2)
y_field = domain.field('y', btype='std', degree=2)

ns.u = u_field
ns.v = v_field
ns.p = p_field
ns.t = t_field
ns.y = y_field

test_u = domain.field('tu', btype='std', degree=2)
test_v = domain.field('tv', btype='std', degree=2)
test_p = domain.field('tp', btype='std', degree=1)
test_t = domain.field('tt', btype='std', degree=2)
test_y = domain.field('ty', btype='std', degree=2)


# --- Variational Formulation (with combustion source term placeholder) ---
try:
    print("--- DEBUG: Setting up variational formulation for reacting flow. ---")
    # Arrhenius parameters (example values, adjust for real rocket chemistry)
    A = 1e6  # Pre-exponential factor [1/s]
    E = 8e4  # Activation energy [J/mol]
    R = 8.314  # Gas constant [J/(mol*K)]
    deltaH = 2e7  # Heat of reaction [J/kg]
    Y_fuel_inlet = 1.0  # Inlet fuel mass fraction

    # Reaction rate (Arrhenius law, single-step, irreversible)
    ns.omega = A * ns.y * function.exp(-E/(R*(ns.t+1e-8)))
    ns.Q_comb = -ns.omega * deltaH

    # Navier-Stokes + energy + species (steady, low Mach, 2D)
    # Momentum x
    res_u = domain.integral('rho u d_j(u) tu + d_j(p) tu - mu d_j(u) d_j(tu) dV' @ ns, degree=4)
    # Momentum y
    res_v = domain.integral('rho v d_j(v) tv + d_j(p) tv - mu d_j(v) d_j(tv) dV' @ ns, degree=4)
    # Continuity
    res_p = domain.integral('(d_j(u) + d_j(v)) tp dV' @ ns, degree=4)
    # Energy
    res_t = domain.integral('(rho cp (u d_0(t) + v d_1(t)) tt + k d_j(t) d_j(tt) - Q_comb tt) dV' @ ns, degree=4)
    # Species
    D = 1e-5  # Diffusivity [m^2/s]
    res_y = domain.integral('(rho (u d_0(y) + v d_1(y)) ty + D d_j(y) d_j(ty) + omega ty) dV' @ ns, degree=4)

    res = res_u + res_v + res_p + res_t + res_y

    # Boundary conditions
    print("--- DEBUG: Setting up Dirichlet BCs for reacting flow. ---")
    sqr_t = domain.boundary[inlet_boundary_name].integral('(t - Tinlet)^2 dS' @ ns, degree=4)
    sqr_y = domain.boundary[inlet_boundary_name].integral('(y - Y_fuel_inlet)^2 dS' @ ns, degree=4)
    # Inlet velocity (u,v)
    sqr_u = domain.boundary[inlet_boundary_name].integral('(u - velocity_0)^2 dS' @ ns, degree=4)
    sqr_v = domain.boundary[inlet_boundary_name].integral('(v - velocity_1)^2 dS' @ ns, degree=4)
    # Outlet: zero pressure
    sqr_p = domain.boundary[outlet_boundary_name].integral('p^2 dS' @ ns, degree=4)
    cons = solver.System(sqr_t + sqr_y + sqr_u + sqr_v + sqr_p, trial=['t','y','u','v','p']).solve_constraints(droptol=1e-15)

    # --- preCICE Coupling Setup (temperature only) ---
    bases = {name: domain.boundary[name].basis('std', degree=2) for name in coupling_boundaries}
    ns.basesWall = bases["Wall"]
    ns.ptWall = function.Argument('ptWall', shape=(len(ns.basesWall),))
    ns.interpTempWall = ns.ptWall @ ns.basesWall
    # Add coupling condition (temperature)
    res += domain.boundary["Wall"].integral('tt (t - interpTempWall) dS' @ ns, degree=4)
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
    x_coords = domain.boundary[name].project(ns.x[0], onto=boundary_basis, geometry=ns.x, ischeme='gauss4')
    y_coords = domain.boundary[name].project(ns.x[1], onto=boundary_basis, geometry=ns.x, ischeme='gauss4')
    z_coords = np.zeros_like(x_coords)  # 2D mesh, z=0
    coords = np.stack([x_coords, y_coords, z_coords], axis=-1)
    print(f"--- DEBUG: Projected {coords.shape[0]} coordinates for {name} ---")
    all_coords.append(coords)
all_coords_np = np.concatenate(all_coords)
print(f"--- DEBUG: Total projected coordinates for all coupling boundaries: {all_coords_np.shape} ---")

print(f"--- DEBUG: Registering {all_coords_np.shape[0]} mesh vertices with preCICE ---")
vertex_ids = interface.set_mesh_vertices(mesh_name, all_coords_np)

split_indices = [len(bases[name]) for name in coupling_boundaries]
split_points = np.cumsum(split_indices)[:-1]
print(f"--- DEBUG: Split indices for boundaries: {split_indices} ---")
print(f"--- DEBUG: Split points for np.split: {split_points} ---")

print("--- DEBUG: Initializing preCICE interface ---")
interface.initialize()
print("--- DEBUG: preCICE interface initialized. ---")

time_step_number = 0
num_vertices = len(vertex_ids)
precice_temp_all = np.zeros(num_vertices)

while interface.is_coupling_ongoing():
    print("--- DEBUG: Starting new preCICE coupling iteration ---")
    dt = interface.get_max_time_step_size()
    print(f"--- DEBUG: preCICE max time step size: {dt}")

    if interface.requires_writing_checkpoint():
        print("--- DEBUG: preCICE action write iteration checkpoint required ---")
        saved_lhs = lhs.copy()

    precice_temp_all = interface.read_data(mesh_name, read_data_name, vertex_ids, 0.0 * dt)
    temps_wall, = np.split(precice_temp_all, split_points) if split_points.size else (precice_temp_all,)
    print(f"--- DEBUG: Received temperature data shape: {precice_temp_all.shape}")
    print(f"--- DEBUG: Split temperature data: wall={temps_wall.shape}")

    arguments = dict(ptWall=temps_wall)
    arguments[''] = lhs
    solution_dict = solver.System(res, trial='t', test='v').solve(constrain=cons, arguments=arguments)
    lhs = solution_dict['t']
    print("--- DEBUG: Linear system solved ---")

    # --- FLUX CALCULATION ---
    flux_expr = '-k d_i(t) n_i' @ ns
    eval_args = {'t': lhs}
    flux_wall = domain.boundary["Wall"].project(flux_expr, onto=bases["Wall"], geometry=ns.x, ischeme='gauss4', arguments=eval_args)
    print(f"--- DEBUG: Heat flux shape: wall={flux_wall.shape}")

    heat_flux_all = flux_wall
    interface.write_data(mesh_name, write_data_name, vertex_ids, heat_flux_all)
    print(f"--- DEBUG: Writing heat flux data to preCICE, shape: {heat_flux_all.shape} ---")

    print("--- DEBUG: Advancing preCICE interface ---")
    interface.advance(dt)

    if interface.requires_reading_checkpoint():
        print("--- DEBUG: preCICE action read iteration checkpoint required ---")
        lhs = saved_lhs.copy()

    # --- VTK Export ---
    print("--- DEBUG: Exporting combustion solution to VTK ---")
    vtk_sample = domain.sample('vtk', 1.0)
    points = vtk_sample.eval(ns.x)
    temperature_values = vtk_sample.eval(ns.t, arguments=eval_args)
    heat_flux_values = vtk_sample.eval('-k d_i(t)' @ ns, arguments=eval_args)
    filename = f"combustion_solution_{time_step_number}"
    full_path = f"{results_folder}/{filename}"
    time_step_number += 1
    export.vtk(full_path, vtk_sample.tri, points, T=temperature_values, H=heat_flux_values)

print("--- DEBUG: Finalizing preCICE interface ---")
interface.finalize()
print("Combustion solver finished.")
