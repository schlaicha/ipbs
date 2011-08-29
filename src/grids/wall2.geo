// Gmsh project created on Tue Aug 23 17:13:33 2011
// Creates 2 charged walls for iPBS using spherical symmetry
// refinement is done via the characteristic length

box_length = 60;
box_width = 30;
wall_width = 1;
wall_height = 50;
wall_refinement = .1;
wall_distance = 1.5;
lower_left_refinement =1;
upper_left_refinement = 3;
lower_right_refinement = 1.;
upper_right_refinement = 3;
center_refinement = .2;
middle_refinement = 2;

// Define the geometry for the wall (iPBS)

Point(0) = {0, 0, 0, center_refinement};
Point(1) = {-box_width/2, 0, 0, lower_left_refinement};
Point(2) = {-box_width/2, box_length, 0, upper_left_refinement};
Point(3) = {-box_width/2, box_length/2, 0, middle_refinement};
Point(4) = {box_width/2, 0, 0, lower_right_refinement};
Point(5) = {box_width/2, box_length, 0, upper_right_refinement};
Point(7) = {box_width/2, box_length/2, 0, middle_refinement};
Point(8) = {0, box_length, 0, middle_refinement};
Point(9) = {-wall_distance-wall_width/2, 0, 0, wall_refinement};
Point(10) = {-wall_distance+wall_width/2, 0, 0, wall_refinement};
Point(11) = {-wall_distance-wall_width/2, wall_height, 0, wall_refinement};
Point(12) = {-wall_distance+wall_width/2, wall_height, 0, wall_refinement};
Point(13) = {-wall_distance, wall_height, 0, wall_refinement};
Point(14) = {-wall_distance, wall_height+wall_width/2, 0, wall_refinement};
Point(15) = {wall_distance-wall_width/2, 0, 0, wall_refinement};
Point(16) = {wall_distance+wall_width/2, 0, 0, wall_refinement};
Point(17) = {wall_distance-wall_width/2, wall_height, 0, wall_refinement};
Point(18) = {wall_distance+wall_width/2, wall_height, 0, wall_refinement};
Point(19) = {wall_distance, wall_height, 0, wall_refinement};
Point(20) = {wall_distance, wall_height+wall_width/2, 0, wall_refinement};
Line(2) = {1, 3};
Line(3) = {3, 2};
Line(4) = {2, 8};
Line(5) = {8, 5};
Line(6) = {5, 7};
Line(7) = {7, 4};
Line(8) = {4, 16};
Line(9) = {16, 18};
Line(10) = {17, 15};
Line(11) = {15, 0};
Line(21) = {0, 10};
Line(12) = {10, 12};
Line(13) = {11, 9};
Line(14) = {9, 1};
Circle(15) = {18, 19, 20};
Circle(16) = {20, 19, 17};
Circle(17) = {12, 13, 14};
Circle(18) = {14, 13, 11};
Line Loop(22) = {6, 7, 8, 9, 15, 16, 10, 11, 21, 12, 17, 18, 13, 14, 2, 3, 4, 5};
Plane Surface(23) = {22};

// Define physical groups for assigning B.C.
// 0 is for Dirichlet boundary elements
// 1 for Neumann
// 2 for iPBS
Physical Line(2) = {9, 10, 15, 16};
Physical Line(3) = {12, 17, 18, 13};
Physical Line(1) = {14, 8, 2, 3, 4, 5, 6, 7, 11, 21};

// At least one physical surface (or physical volume
// is needed by DUNE gmshreader (otherwise all elements
// are boundaries)
Physical Surface(25) = {23};
