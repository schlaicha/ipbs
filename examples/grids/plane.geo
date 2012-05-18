// Gmsh project created on Mon 23 May 2011
// Infinite plane (Gouy-Chapman model) for iPBS

// Set geometry
box_length = 100;
box_width = 12;
wall_refinement = .1;
outer_refinement = 1;

// Setup mesh
Point(1) = {0, 0, 0, wall_refinement};
Point(2) = {box_width, 0, 0, outer_refinement};
Point(3) = {0, box_length, 0, wall_refinement};
Point(4) = {box_width, box_length, 0, outer_refinement};
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {1, 3};
Line(4) = {3, 4};
Line Loop(7) = {4, -2, -1, 3};
Plane Surface(8) = {7};

// Define boundaries
Physical Line(0) = {1,2,4};
Physical Line(1) = {3};
Physical Surface(12) = {8};
