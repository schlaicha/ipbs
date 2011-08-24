// Gmsh project created on Tue Aug 23 17:13:33 2011
// Creates 2 charged walls for iPBS using spherical symmetry
// refinement is done via the characteristic length

box_length = 100;
box_width = 300;
wall_distance = 50;
wall_width = 50;
wall_height = 01;
wall_refinement = .4;
lower_left_refinement =4.;
upper_left_refinement = 10.;
lower_right_refinement = 8.;
upper_right_refinement = 10.;
center_refinement = 5.0;

// Define the geometry for the wall (iPBS)

Point(1) = {-box_width/2, 0, 0, lower_left_refinement};
Point(2) = {-box_width/2, box_length, 0, upper_left_refinement};
Point(3) = {-box_width/2, box_length/2, 0, center_refinement};
Point(4) = {box_width/2, 0, 0, lower_right_refinement};
Point(5) = {box_width/2, box_length, 0, upper_right_refinement};
Point(7) = {box_width/2, box_length/2, 0, center_refinement};
Point(8) = {0, box_length, 0, center_refinement};
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
Line(1) = {1, 1};
Line(2) = {1, 3};
Line(3) = {3, 2};
Line(4) = {2, 8};
Line(5) = {8, 5};
Line(6) = {5, 7};
Line(7) = {7, 4};
Line(8) = {4, 16};
Line(9) = {16, 18};
Line(10) = {17, 15};
Line(11) = {15, 10};
Line(12) = {10, 12};
Line(13) = {11, 9};
Line(14) = {9, 1};
Circle(15) = {18, 19, 20};
Circle(16) = {20, 19, 17};
Circle(17) = {12, 13, 14};
Circle(18) = {14, 13, 11};
Line Loop(19) = {5, 6, 7, 8, 9, 15, 16, 10, 11, 12, 17, 18, 13, 14, 2, 3, 4};
Plane Surface(20) = {19};

// Define physical groups for assigning B.C.
// 0 is for Dirichlet boundary elements
// 1 for Neumann
// 2 for iPBS
Physical Line(2) = {9, 10, 15, 16};
Physical Line(3) = {12, 17, 18, 13};
Physical Line(1) = {14, 11, 8};
Physical Line(0) = {7, 6, 5, 4, 3, 2};

// At least one physical surface (or physical volume
// is needed by DUNE gmshreader (otherwise all elements
// are boundaries)

