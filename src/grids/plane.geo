// Gmsh project created on Thu 26 May 2011
// Infinite plane (Gouy-Chapman model) for iPBS

// define characteristic length (mesh quality)
lc = 0.1;

// Set geometry
box_size = 20;

// Setup mesh
Point(1) = {-box_size/2, 0, 0};
Point(2) = {box_size/2, 0, 0};
Point(3) = {-box_size/2, box_size, 0};
Point(4) = {box_size/2, box_size, 0};
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {1, 3};
Line(4) = {3, 4};
Line Loop(7) = {4, -2, -1, 3};
Plane Surface(8) = {7};

// Define boundaries
Physical Line(0) = {4};
Physical Line(1) = {3, 2};
Physical Line(2) = {1};
Physical Surface(12) = {8};
