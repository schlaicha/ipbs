// Gmsh project created on Tue Aug 23 10:26:33 2011
// Creates a charged wall for iPBS using spherical symmetry
// refinement is done via the characteristic length

box_length = 150;
box_width = 30;
wall_position = 0;
wall_width = 2;
wall_height = 140;
wall_refinement = .1;
lower_left_refinement =4.;
upper_left_refinement = 10.;
lower_right_refinement = 8.;
upper_right_refinement = 10.;
center_refinement = 0.1;

// Define the geometry for the wall (iPBS)

Point(1) = {-box_width/2, 0, 0, lower_left_refinement};
Point(2) = {-box_width/2, box_length, 0, upper_left_refinement};
Point(3) = {-box_width/2, box_length/2, 0, 1};
Point(4) = {box_width/2, 0, 0, lower_right_refinement};
Point(5) = {box_width/2, box_length, 0, upper_right_refinement};
Point(7) = {box_width/2, box_length/2, 0, 1};
Point(8) = {0, box_length, 0, center_refinement};
Point(9) = {wall_position-wall_width/2, 0, 0, wall_refinement};
Point(10) = {wall_position+wall_width/2, 0, 0, wall_refinement};
Point(11) = {wall_position-wall_width/2, wall_height, 0, wall_refinement};
Point(12) = {wall_position+wall_width/2, wall_height, 0, wall_refinement};
Point(13) = {wall_position, wall_height, 0, wall_refinement};
Point(14) = {wall_position, wall_height+wall_width/2, 0, wall_refinement};
Line(1) = {1, 3};
Line(2) = {2, 3};
Line(3) = {2, 8};
Line(4) = {8, 5};
Line(5) = {5, 7};
Line(6) = {7, 4};
Line(7) = {4, 10};
Line(8) = {10, 12};
Line(9) = {1, 9};
Line(10) = {9, 11};
Circle(11) = {11, 13, 14};
Circle(12) = {14, 13, 12};
Line Loop(13) = {5, 6, 7, 8, -12, -11, -10, -9, 1, -2, 3, 4};
Plane Surface(14) = {13};

// Define physical groups for assigning B.C.
// 0 is for Dirichlet boundary elements
// 1 for Neumann
// 2 for iPBS
Physical Line(1) = {9, 1, 2, 3, 4, 5, 6, 7};
Physical Line(2) = {8};
Physical Line(3) = {10};
Physical Line(4) = {11, 12};

// At least one physical surface (or physical volume
// is needed by DUNE gmshreader (otherwise all elements
// are boundaries)
Physical Surface(19) = {14};
