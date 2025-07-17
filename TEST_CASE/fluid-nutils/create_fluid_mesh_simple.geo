// This script loads a BREP file for the fluid domain and defines physical groups.
// It is designed to work with the Gmsh mesh generator.

// Load the geometry from your .brep file
Merge "U_Channel_Fluid_mesh.brep";

// Define a mesh size parameter
lc = 0.001; // Adjusted for finer mesh, change as needed
lm = 0.0001; // Minimum mesh size for finer details
Mesh.CharacteristicLengthMin = lm;
Mesh.CharacteristicLengthMax = lc;

// --- Name the Volume ---
Physical Volume("U_Channel_Fluid_Volume") = {1};

// --- Name ALL Boundaries ---
// You must find the correct surface numbers (tags) in the Gmsh GUI.
// Replace the placeholder numbers with the real tags from your geometry.

// Coupling Surfaces
Physical Surface("U_Channel_Fluid_Channel_Top") = {3};    // This seems to be tag 2 in your file
Physical Surface("U_Channel_Fluid_Channel_Right") = {2};  // This seems to be tag 3 in your file
Physical Surface("U_Channel_Fluid_Channel_Bottom") = {1}; // This seems to be tag 4 in your file

// Other Boundaries (THESE ARE MISSING)
Physical Surface("U_Channel_Fluid_Entry") = {6}; // Example: replace 5 with the correct tag
Physical Surface("U_Channel_Fluid_Outlet") = {5}; // Example: replace 6 with the correct tag

// Generate the 3D mesh
Mesh 3;

// Save the final mesh file.
Mesh.MshFileVersion = 4.1;
Save "U_Channel_Fluid_mesh.msh";