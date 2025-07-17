# Fluid Solver (fluid-nutils)

This folder contains the fluid-side implementation for the coupled heat transfer simulation using [Nutils](https://nutils.org/) and [preCICE](https://www.precice.org/).

## Contents

- `fluid_solver.py` — Main fluid solver logic, including coupling with preCICE. Not physcially accurate for the moment.
- `FluidSolver.py` — Example or alternative fluid solver script.
- `nutils_solver.py` — Additional solver logic or utilities.
- `create_fluid_mesh.geo` — GMSH script to generate the fluid mesh from the BREP geometry.
- `U_Channel_Fluid_mesh.brep` — Simple BREP geometry for the fluid domain.
- `Cooling_Channels_mesh.brep` — More complex BREP geometry for advanced cases.
- `U_Channel_Fluid_mesh.msh` — Generated mesh file (output of GMSH).
- `run.sh` — Example script to run the fluid solver.

## How it Works

- The `.geo` file loads the `.brep` geometry and defines physical groups (named boundaries) for coupling and boundary conditions.
- The mesh is generated with GMSH and used by the Nutils-based solver.
- The solver exchanges temperature and heat flux data with the solid solver via preCICE at the defined coupling boundaries.

## Modifying the Fluid Solver

- **Geometry**: Edit the `.brep` file and update the `.geo` script as needed. Regenerate the mesh with GMSH.
- **Coupling Boundaries**: Ensure physical group names in the `.geo`, mesh, code, and preCICE config all match.
- **Physics**: Adjust the variational formulation in `fluid_solver.py` as needed.

## Tips

- Use GMSH to inspect and verify boundary tags and physical group assignments.
- Keep boundary names consistent across all files.
