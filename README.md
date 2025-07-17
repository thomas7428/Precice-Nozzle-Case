# Precice-Nozzle-Case

This repository contains simulation cases and scripts for coupled fluid-solid simulations using preCICE and nutils. The structure is organized to separate fluid and solid domains, as well as mesh generation and solver scripts.

## Structure

- `TEST_CASE/`
  - `config.json` — Configuration for the test case.
  - `precice-config.xml` — preCICE configuration file.
  - `fluid-nutils/` — Fluid domain mesh files, mesh generation scripts, and solver scripts.
  - `solid-nutils/` — Solid domain mesh files, mesh generation scripts, and solver scripts.
  - `combustion-nutils/` — (if present) Scripts and meshes for combustion simulations.

## Getting Started

1. **Install Dependencies**
   - Python (recommended: 3.8+)
   - nutils
   - preCICE
   - Any other dependencies as required by the solver scripts.

2. **Generate Meshes**
   - Use the `.geo` scripts in `fluid-nutils/` and `solid-nutils/` to generate the required mesh files (e.g., with Gmsh).

3. **Run Simulations**
   - Use the provided `run.sh` scripts in each domain folder to start the respective solvers.
   - Make sure to configure `precice-config.xml` as needed for your coupling setup.

## Notes

- See the `README.md` files in each subfolder for more details on running specific solvers or generating meshes.
- This repository is intended for research and educational purposes.
