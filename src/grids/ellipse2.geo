// Gmsh project created on Mon May  30 14:41:33 2011
// Creates the iPBS single sphere with a symmetric mesh
// refinement is done via the characteristic length

a = 3;
b = 25;
position = 0;  // Center position of the sphere
box_size = 50;
outer_refinement = 5;
inner_refinement = 0.2;

// Define the geometry for 2d-sphere (iPBS)

Point(1) = {position, 0, 0, 1.0};
Point(2) = {position, b, 0, inner_refinement/3.};
Point(3) = {position-a, 0, 0, inner_refinement};
Point(4) = {position+a, 0, 0, inner_refinement};
Point(5) = {-box_size, 0, 0, outer_refinement};
Point(6) = {box_size, 0, 0, outer_refinement};
Point(7) = {-box_size, box_size, 0, outer_refinement};
Point(8) = {box_size, box_size, 0, outer_refinement};
Line(3) = {4, 6};
Line(4) = {3, 5};
Line(5) = {5, 7};
Line(6) = {7, 8};
Line(7) = {8, 6};
Ellipse(8) = {4, 1, 2, 2};
Ellipse(9) = {3, 1, 2, 2};
Line Loop(10) = {6, 7, -3, 8, -9, 4, 5};
Plane Surface(11) = {10};

// Define physical groups for assigning B.C.
// 0 is for Dirichlet boundary elements
// 1 for Neumann
// 2 for iPBS
Physical Line(2) = {8, 9};
Physical Line(1) = {3, 4};
Physical Line(0) = {5, 6, 7};

// At least one physical surface (or physical volume
// is needed by DUNE gmshreader (otherwise all elements
// are boundaries)
Physical Surface(15) = {11};
