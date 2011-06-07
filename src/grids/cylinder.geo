// Gmsh project created on Mon 23 May 2011
// Infinite cylinder for iPBS

// define characteristic length (mesh quality)
lc = 0.1;

// Set geometry
outer_radius = 10;
length = 20;
inner_radius = 0.5;

// Setup mesh
Point(1) = {-length/2, inner_radius, 0};
Point(2) = {length/2, inner_radius, 0};
Point(3) = {-length/2, outer_radius, 0};
Point(4) = {length/2, outer_radius, 0};
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {1, 3};
Line(4) = {3, 4};
Line Loop(7) = {4, -2, -1, 3};
Plane Surface(8) = {7};

// Define boundaries
Physical Line(1) = {3, 4, 2};
Physical Line(2) = {1};
Physical Surface(11) = {8};
