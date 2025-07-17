// This script loads a BREP file for the solid domain and defines physical groups
// that match the provided U_Channel_Solid_mesh.msh file.

// Merge the BREP file which contains the geometry
Merge "Solid_mesh.brep";

// Define a mesh size parameter
lc = 0.005; // Adjusted for finer mesh, change as needed
Mesh.CharacteristicLengthMax = lc;

// --- Define Physical Groups ---
// The names and the second argument (e.g., 2, 3, 4...) are the physical tags.
// They are taken directly from your U_Channel_Solid_mesh.msh file.
// The numbers in curly braces {} are the geometric entity tags from the BREP file.
// These may need to be verified by inspecting the BREP file in the Gmsh GUI.

// Volume
Physical Volume("U_Channel_Solid_Volume", 4) = {1};

// Surfaces
Physical Surface("U_Channel_Solid_Entry", 2) = {2}; // Entry surface of the solid
Physical Surface("U_Channel_Solid_Exit", 3) = {3}; // Exit surface of the solid
Physical Surface("U_Channel_Solid_Bottom", 5) = {10}; // Bottom surface of the solid in contact with the hot gaz
Physical Surface("U_Channel_Solid_Right", 6) = {5}; // Right surface of the solid
Physical Surface("U_Channel_Solid_Top", 7) = {1}; // Top surface of the solid in contact with the exterior
Physical Surface("U_Channel_Solid_Left_Top", 8) = {4}; // Left face on top of the cooling channel
Physical Surface("U_Channel_Solid_Channel_Top", 9) = {9}; // Top surface of the channel
Physical Surface("U_Channel_Solid_Channel_Right", 10) = {8}; // Right wall of the channel
Physical Surface("U_Channel_Solid_Channel_Bottom", 11) = {7}; // Bottom of the channel
Physical Surface("U_Channel_Solid_Left_Bottom", 12) = {6}; // Left face under the cooling channel 

// Generate 3D mesh
Mesh 3;

// Save the mesh in the correct format (Version 4.1)
Mesh.MshFileVersion = 4.1;
Save "U_Channel_Solid_mesh.msh";