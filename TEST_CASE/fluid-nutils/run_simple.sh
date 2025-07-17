#!/bin/bash
#
# Run script for the Fluid participant (nutils)

echo "Generating fluid mesh..."
# Generate the mesh, ensuring the output filename is correct
gmsh create_fluid_mesh_simple.geo -3 -o U_Channel_Fluid_mesh.msh

echo "Starting Fluid solver..."
python FluidSolver.py