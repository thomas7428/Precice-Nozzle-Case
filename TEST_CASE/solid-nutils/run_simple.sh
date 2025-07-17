#!/bin/bash
#
# Run script for the Solid participant (dummy)

echo "Generating solid mesh..."
# Generate the mesh, ensuring the output filename is correct
gmsh create_solid_mesh_simple.geo -3 -o U_Channel_Solid_mesh.msh

echo "Starting Solid solver..."
python SolidSolver.py