// Gmsh project created on Mon May  30 14:41:33 2011
// Creates the iPBS single sphere with a symmetric mesh
// refinement is done via the characteristic length

a = 12;
b = 2;
position = 0;  // Center position of the sphere
box_size = 50;
outer_refinement = 5;
inner_refinement = 0.5;

// Define the geometry for 2d-sphere (iPBS)

Point(1) = {position, 0, 0, 1.0};
Point(2) = {position, b/2, 0, inner_refinement};
Point(3) = {position-a/2, 0, 0, inner_refinement};
Point(4) = {position+a/2, 0, 0, inner_refinement};
Point(5) = {-box_size, 0, 0, outer_refinement};
Point(6) = {box_size, 0, 0, outer_refinement};
Point(7) = {-box_size, box_size, 0, outer_refinement};
Point(8) = {box_size, box_size, 0, outer_refinement};
Line(3) = {4, 6};
Line(4) = {3, 5};
Line(5) = {5, 7};
Line(6) = {7, 8};
Line(7) = {8, 6};
Ellipse(12) = {2, 1, 2, 4};
Ellipse(13) = {2, 1, 2, 3};
Line Loop(14) = {12, 3, -7, -6, -5, -4, -13};
Plane Surface(15) = {14};

// Define physical groups for assigning B.C.
// 0 is for Dirichlet boundary elements
// 1 for Neumann
// 2 for iPBS
Physical Line(2) = {12, 13};
Physical Line(1) = {4, 3};
Physical Line(0) = {7, 5, 6};

// At least one physical surface (or physical volume
// is needed by DUNE gmshreader (otherwise all elements
// are boundaries)
Physical Surface(19) = {15};
