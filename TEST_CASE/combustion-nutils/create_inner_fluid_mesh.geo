// create_inner_fluid_mesh.geo
// Gmsh script to generate the inner fluid mesh for the rocket nozzle combustion domain

// Parameters
lc = 0.0025; // Mesh size, adjust as needed

// Geometry definition (placeholder: replace with actual nozzle geometry)
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 0.2, 0, lc};
Point(4) = {0, 0.2, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// Physical groups
Physical Surface("Inner_Fluid_Domain", 1) = {6};
Physical Line("Inlet", 2) = {1};
Physical Line("Outlet", 3) = {3};
Physical Line("Wall", 4) = {2, 4};

// Mesh generation
Mesh 2;
Mesh.MshFileVersion = 4.1;
Save "Inner_Fluid_mesh.msh";
