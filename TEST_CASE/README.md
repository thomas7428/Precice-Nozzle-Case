
# TEST_CASE: Template for 1D Fluid-Solid Coupling with preCICE

This folder provides a ready-to-use template for coupled heat transfer simulations between a 1D fluid and a solid domain using [preCICE](https://www.precice.org/) and [Nutils](https://nutils.org/). Both the fluid and solid solvers are implemented in Python and are located in the `fluid-nutils` and `solid-nutils` folders, respectively.

## Overview

The template demonstrates a partitioned simulation where the fluid and solid domains exchange heat fluxes and temperatures at their interface boundaries via preCICE. The geometry and mesh for each domain are defined using `.brep` and `.geo` files, and meshing is performed with GMSH.

## Folder Structure

- `config.json` — Simulation parameters for both solvers (material properties, time settings, etc.).
- `precice-config.xml` — preCICE configuration file (participants, meshes, data mapping, coupling scheme).
- `solid-nutils/` — Solid solver code, mesh generation script (`create_solid_mesh.geo`), and geometry (`U_Channel_Solid_mesh.brep` for simple case and `Solid_mesh.brep` for a more concrete case).
- `fluid-nutils/` — Fluid solver code, mesh generation script (`create_fluid_mesh.geo`), and geometry (`U_Channel_Fluid_mesh.brep` for simple case and `Cooling_Channels_mesh.brep` for a more concrete case).
- `U_Channel_Solid_mesh.msh` / `U_Channel_Fluid_mesh.msh` — GMSH mesh files generated from the `.geo` scripts.

## Mesh Import and Geometry Definition

**Mesh importation is done via the `.geo` files in each solver folder**, which load the corresponding `.brep` geometry file. The `.geo` scripts define mesh sizes and assign physical groups (named boundaries) that are referenced in the code and preCICE configuration.

To modify the geometry, edit the `.brep` file and update the `.geo` script as needed. To visualize and identify face numbers or tags, open the `.brep` file in GMSH and use the GUI to display entity numbers. This helps ensure the physical group assignments in the `.geo` file match your intended coupling boundaries.

### Example: Physical Group Naming

In both the solid and fluid `.geo` scripts, physical surfaces are defined for coupling and boundary conditions. For example:

- `U_Channel_Solid_Channel_Top`, `U_Channel_Solid_Channel_Right`, `U_Channel_Solid_Channel_Bottom`
- `U_Channel_Fluid_Channel_Top`, `U_Channel_Fluid_Channel_Right`, `U_Channel_Fluid_Channel_Bottom`

These names must match between the mesh, the solver code, and the preCICE configuration.

## Explanation of Faces

- **Coupling Faces**: The `*_Channel_Top`, `*_Channel_Right`, and `*_Channel_Bottom` surfaces are used for coupling between the fluid and solid domains. Heat fluxes and temperatures are exchanged here via preCICE.
- **Other Faces**: Additional surfaces (e.g., `Entry`, `Exit`, `Left_Top`, `Left_Bottom`) are used for boundary conditions such as inlets, outlets, insulation, or convection.

## Simulation Logic

1. **Initialization**:
   - Each solver loads its mesh (generated from the `.geo` and `.brep` files) and sets up the variational problem using Nutils.
   - preCICE is initialized and coupling boundaries are registered in both solvers.

2. **Coupling**:
   - At each time step, the fluid and solid solvers exchange data (heat fluxes and temperatures) at the coupling boundaries using preCICE.
   - The solid solver projects received heat fluxes onto its boundary basis and solves for the new temperature field.
   - The fluid solver receives updated temperatures and advances its own solution accordingly.

3. **Time Stepping**:
   - The simulation advances in time, repeating the coupling and solution steps until the end time is reached.

4. **Output**:
   - Results (e.g., temperature fields, heat fluxes) are exported in VTK format for visualization by both solvers.

## How to Modify the Template

- **Change the Geometry**: Edit the `.brep` file and update the `.geo` script in the relevant solver folder. Regenerate the mesh with GMSH.
- **Update Coupling Boundaries**: Ensure physical group names in the `.geo` scripts, mesh, code, and preCICE config all match.
- **Modify Physics**: Adjust the variational formulation in the solver code as needed.
- **Change Simulation Parameters**: Edit `config.json` for time step size, end time, material properties, etc.
- **Update preCICE Configuration**: Edit `precice-config.xml` to change coupling data, mapping, or participant settings.

## Tips

- Always check that boundary names are consistent across mesh, code, and preCICE config.
- Use GMSH to visually inspect and verify face numerotation and physical group assignments.
- The template is designed for easy extension to more complex scenarios (e.g., 2D/3D, different physics).

## References

- [preCICE Documentation](https://precice.org/docs/)
- [Nutils Documentation](https://docs.nutils.org/)
- [GMSH Documentation](http://gmsh.info/doc/texinfo/gmsh.html)

---

For further questions or to contribute improvements, please open an issue or pull request in the main repository.
