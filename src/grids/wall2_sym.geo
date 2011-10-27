// Gmsh project created on Tue Aug 23 17:13:33 2011
// Creates 2 charged walls for iPBS using spherical symmetry
// refinement is done via the characteristic length

box_length = 15;
box_width = 12;
wall_width = 1;
wall_height = 7;
wall_refinement = .1;
wall_distance = 1.5;
lower_left_refinement =1;
upper_left_refinement = 1;
lower_right_refinement = 1.;
upper_right_refinement = 1;
center_refinement = .1;
middle_refinement = 1;

// Define the geometry for the wall (iPBS)

Point(0) = {0, 0, 0, center_refinement};
Point(1) = {-box_width/2, 0, 0, lower_left_refinement};
Point(2) = {-box_width/2, box_length, 0, upper_left_refinement};
Point(3) = {-box_width/2, box_length/2, 0, middle_refinement};
Point(8) = {0, box_length, 0, middle_refinement};
Point(9) = {-wall_distance-wall_width/2, 0, 0, wall_refinement};
Point(10) = {-wall_distance+wall_width/2, 0, 0, wall_refinement};
Point(11) = {-wall_distance-wall_width/2, wall_height, 0, wall_refinement};
Point(12) = {-wall_distance+wall_width/2, wall_height, 0, wall_refinement};
Point(13) = {-wall_distance, wall_height, 0, wall_refinement};
Point(14) = {-wall_distance, wall_height+wall_width/2, 0, wall_refinement};
Line(1) = {1, 3};
Line(2) = {3, 2};
Line(3) = {2, 8};
Line(4) = {8, 0};
Line(5) = {0, 10};
Line(6) = {10, 12};
Line(7) = {9, 1};
Line(8) = {9, 11};
Circle(9) = {11, 13, 14};
Circle(10) = {14, 13, 12};
Line Loop(11) = {4, 5, 6, -10, -9, -8, 7, 1, 2, 3};
Plane Surface(12) = {11};

// Define physical groups for assigning B.C.
// 0 is for Dirichlet boundary elements
// 1 for Neumann
// 2 for iPBS
Physical Line(1) = {4, 7, 1, 2, 3, 5};
Physical Line(2) = {6, 8, 9, 10};

// At least one physical surface (or physical volume
// is needed by DUNE gmshreader (otherwise all elements
// are boundaries)
Physical Surface(15) = {12};
