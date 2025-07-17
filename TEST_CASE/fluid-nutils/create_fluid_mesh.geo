// This script loads a BREP file for the fluid domain and defines physical groups.
// It is designed to work with the Gmsh mesh generator.

// Load the geometry from your .brep file
Merge "Cooling_Channels_mesh.brep";

// Define a mesh size parameter
lc = 0.005; // Adjusted for finer mesh, change as needed
Mesh.CharacteristicLengthMax = lc;

// --- Name the Volume ---
Physical Volume("U_Channel_Fluid_Volume") = {1};

// --- Name ALL Boundaries ---
// You must find the correct surface numbers (tags) in the Gmsh GUI.
// Replace the placeholder numbers with the real tags from your geometry.

// Coupling Surfaces
Physical Surface("U_Channel_Fluid_Channel_Top") = {1};    // This seems to be tag 2 in your file
Physical Surface("U_Channel_Fluid_Channel_Right") = {2};  // This seems to be tag 3 in your file
Physical Surface("U_Channel_Fluid_Channel_Bottom") = {6}; // This seems to be tag 4 in your file

// Other Boundaries (THESE ARE MISSING)
Physical Surface("U_Channel_Fluid_Entry") = {5}; // Example: replace 5 with the correct tag
Physical Surface("U_Channel_Fluid_Outlet") = {3}; // Example: replace 6 with the correct tag

// Generate the 3D mesh
Mesh 3;

// Save the final mesh file.
Save "U_Channel_Fluid_mesh.msh";