// Gmsh project created on Mon May  30 14:41:33 2011
// Creates the iPBS single sphere with a symmetric mesh
// refinement is done via the characteristic length

radius = 3;
position = 0;  // Center position of the sphere
box_size = 20;
outer_refinement =2.;
middle_refinement =1.;
inner_refinement = 0.1;

// Define the geometry for 2d-sphere (iPBS)

Point(1) = {position, 0, 0, 1.0};
Point(2) = {position, radius, 0, inner_refinement};
Point(3) = {position-radius, 0, 0, inner_refinement};
Point(4) = {position+radius, 0, 0, inner_refinement};
Point(5) = {-box_size, 0, 0, outer_refinement};
Point(6) = {box_size, 0, 0, outer_refinement};
Point(7) = {0, box_size, 0, outer_refinement};
Point(8) = {position+radius+box_size/2, 0, 0, middle_refinement};
Point(9) = {position-radius-box_size/2, 0, 0, middle_refinement};
Point(10) = {0, position+radius+box_size/2, 0, middle_refinement};
Circle(1) = {4, 1, 2};
Circle(2) = {2, 1, 3};
Line(4) = {3, 9};
Line(5) = {4, 8};
Line(10) = {9, 5};
Line(11) = {8, 6};
Circle(6) = {5, 1, 7};
Circle(7) = {7, 1, 6};
Circle(8) = {8, 1, 10};
Circle(9) = {9, 1, 10};

Line Loop(12) = {5, 8, -9, -4, -2, -1};
Plane Surface(13) = {12};
Line Loop(14) = {8, -9, 10, 6, 7, -11};
Plane Surface(15) = {14};

// Define physical groups for assigning B.C.
// 0 is for Dirichlet boundary elements
// 1 for Neumann
// 2 for iPBS
Physical Line(0) = {2, 1};
Physical Line(1) = {5, 11, 4, 10};
Physical Line(2) = {6, 7};


// At least one physical surface (or physical volume
// is needed by DUNE gmshreader (otherwise all elements
// are boundaries)
Physical Surface(19) = {13, 15};
