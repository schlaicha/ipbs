// Gmsh project created on Tue 13 Mar 2012

// Set geometry
sqrt_dist = 6;
radius = 1;
circle_refinement = .1;
corner_refinement = .1;
faraway_refinement = 10;
outer_refinement = 1;
delta = .3;
box_height = 110.;
box_width = 110;

// Setup mesh
Point(0) = {0, 0, 0, corner_refinement};
//Point(1) = {delta, delta, 0, corner_refinement};
Point(2) = {0, box_height, 0, outer_refinement};
Point(3) = {box_width, 0, 0, outer_refinement};
Point(4) = {box_width, box_height, 0, faraway_refinement};
Point(5) = {sqrt_dist, sqrt_dist, 0};
Point(6) = {sqrt_dist-radius, sqrt_dist, 0, circle_refinement};
Point(7) = {sqrt_dist+radius, sqrt_dist, 0, circle_refinement};
Point(8) = {sqrt_dist, sqrt_dist-radius, 0, circle_refinement};
Point(9) = {sqrt_dist, sqrt_dist+radius, 0, circle_refinement};
//Point(10) = {delta, 0, 0, corner_refinement};
//Point(11) = {0, delta, 0, corner_refinement};
Line(1) = {0, 2};
Line(2) = {2, 4};
Line(3) = {0, 3};
Line(4) = {3, 4};
Circle(5) = {7, 5, 9};
Circle(6) = {9, 5, 6};
Circle(7) = {6, 5, 8};
Circle(8) = {8, 5, 7};
//Circle(9) = {11, 1, 10};
Line Loop(10) = {2, -4, -3, 1};
//Line Loop(10) = {2, -4, -3, -9, 1};
Line Loop(11) = {5, 6, 7, 8};
Plane Surface(12) = {10, 11};

// Define boundaries
Physical Line(1) = {2, 4};
Physical Line(2) = {1, 3};
Physical Line(3) = {7, 6, 5, 8};
Physical Surface(13) = {12};
