import numpy as np
from heat_conduction_cylindrical import solve_heat_conduction_steady
import precice

def main():
    meshName = "Fluid-1D-Mesh"
    temperatureName = "Temperature"
    heatfluxName = "Heat-Flux"
    # Example: 20 points along a 1D tube of length 1.0, padded to 3D for preCICE
    n_points = 20
    length = 1.0
    x = np.linspace(0, length, n_points)
    interface_coords = np.zeros((n_points, 3))
    interface_coords[:,0] = x  # x varies, y=z=0
    interface = precice.Participant("Fluid", "../precice-config.xml", 0, 1)
    vertex_ids = interface.set_mesh_vertices(meshName, interface_coords)
    vertex_ids = np.array(vertex_ids, dtype=int)
    temperature = np.ones(len(vertex_ids), dtype=np.float64) * 300.0
    if interface.requires_initial_data():
        interface.write_data(meshName, temperatureName, vertex_ids, temperature)
    interface.initialize()
    checkpoint_temperature = temperature.copy()
    checkpoint_heatflux = np.zeros(len(vertex_ids), dtype=np.float64)
    with open("fluid_temperature_history.txt", "w") as temp_log, open("fluid_heatflux_history.txt", "w") as flux_log:
        temp_log.write("# time " + " ".join([f"T{i}" for i in range(len(vertex_ids))]) + "\n")
        flux_log.write("# time " + " ".join([f"q{i}" for i in range(len(vertex_ids))]) + "\n")
        step = 0
        time = 0.0
        while interface.is_coupling_ongoing():
            dt = interface.get_max_time_step_size()
            print(f"[Fluid] Coupling step {step}, time={time}, dt={dt}")
            print(f"[Fluid] temperature: {temperature}")
            step += 1
            # Read heat flux from solid (for all vertices)
            heatflux = np.zeros(len(vertex_ids), dtype=np.float64)
            for i, vid in enumerate(vertex_ids):
                flux_val = np.zeros(1, dtype=np.float64)
                interface.read_data(meshName, heatfluxName, np.array([vid], dtype=int), flux_val)
                heatflux[i] = flux_val[0]
            print(f"[Fluid] heatflux: {heatflux}")
            # Use received heat flux as Neumann BC if nonzero, else Dirichlet
            if np.any(np.abs(heatflux) > 1e-12):
                interface_bc = heatflux[-1]  # Use last node's heat flux as Neumann BC
                print(f"[Fluid] Using Neumann BC at interface: heatflux={interface_bc}")
                x_nodes_solver, temperature_new = solve_heat_conduction_steady(x, interface_bc, bc_type='neumann')
            else:
                interface_bc = temperature[-1]
                print(f"[Fluid] Using Dirichlet BC at interface: temperature={interface_bc}")
                x_nodes_solver, temperature_new = solve_heat_conduction_steady(x, interface_bc, bc_type='dirichlet')
            temperature_full = np.ones(len(vertex_ids), dtype=np.float64) * 300.0
            for i, xval in enumerate(x_nodes_solver):
                idx = np.argmin(np.abs(x - xval))
                temperature_full[idx] = temperature_new[i]
            temperature = temperature_full
            current_time = interface.get_time() if hasattr(interface, "get_time") else time
            temp_log.write(f"{current_time} " + " ".join(map(str, temperature)) + "\n")
            flux_log.write(f"{current_time} " + " ".join(map(str, heatflux)) + "\n")
            if interface.requires_writing_checkpoint():
                checkpoint_temperature = temperature.copy()
                checkpoint_heatflux = heatflux.copy()
            interface.write_data(meshName, temperatureName, vertex_ids, temperature)
            interface.advance(dt)
            time += dt
            if interface.requires_reading_checkpoint():
                temperature = checkpoint_temperature.copy()
                heatflux = checkpoint_heatflux.copy()
    interface.finalize()

if __name__ == "__main__":
    main()
