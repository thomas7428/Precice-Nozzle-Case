# Solid Solver (solid-nutils)

This folder contains the solid-side implementation for the coupled heat transfer simulation using [Nutils](https://nutils.org/) and [preCICE](https://www.precice.org/).

## Contents

- `solid_solver.py` — Main solid solver logic, including coupling with preCICE. It's more a dummy than a real solver.
- `SolidSolver.py` — Example or alternative solid solver script.
- `create_solid_mesh.geo` — GMSH script to generate the solid mesh from the BREP geometry.
- `U_Channel_Solid_mesh.brep` — Simple BREP geometry for the solid domain.
- `Solid_mesh.brep` — More complex BREP geometry for advanced cases.
- `U_Channel_Solid_mesh.msh` — Generated mesh file (output of GMSH).
- `run.sh` — Example script to run the solid solver.

## How it Works

- The `.geo` file loads the `.brep` geometry and defines physical groups (named boundaries) for coupling and boundary conditions.
- The mesh is generated with GMSH and used by the Nutils-based solver.
- The solver exchanges heat flux and temperature data with the fluid solver via preCICE at the defined coupling boundaries.

## Modifying the Solid Solver

- **Geometry**: Edit the `.brep` file and update the `.geo` script as needed. Regenerate the mesh with GMSH.
- **Coupling Boundaries**: Ensure physical group names in the `.geo`, mesh, code, and preCICE config all match.
- **Physics**: Adjust the variational formulation in `solid_solver.py` as needed.

## Tips

- Use GMSH to inspect and verify boundary tags and physical group assignments.
- Keep boundary names consistent across all files.
