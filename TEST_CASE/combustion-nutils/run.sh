#!/bin/bash
# Run the Nutils-based combustion solver with PreCICE coupling
set -e

# Activate your Python environment if needed
# source ~/venv/bin/activate

# Run the combustion solver
python3 CombustionSolver.py --precice-config ../precice-config.xml --mesh create_inner_fluid_mesh.geo
