// Gmsh project created on Mon May  30 14:41:33 2011
// Creates the iPBS single sphere with a symmetric mesh
// refinement is done via the characteristic length

radius = 5;
position = 0;  // Center position of the sphere
box_size = 30;
outer_refinement = 2;
inner_refinement = 0.1;

// Define the geometry for 2d-sphere (iPBS)

Point(1) = {0, position, 0, 1.0};
Point(2) = {0, radius, 0, inner_refinement};
Point(3) = {-radius, 0, 0, inner_refinement};
Point(4) = {radius, 0, 0, inner_refinement};
Point(5) = {-box_size, 0, 0, outer_refinement};
Point(6) = {box_size, 0, 0, outer_refinement};
Point(7) = {-box_size, box_size, 0, outer_refinement};
Point(8) = {box_size, box_size, 0, outer_refinement};
Circle(1) = {4, 1, 2};
Circle(2) = {2, 1, 3};
Line(3) = {4, 6};
Line(4) = {3, 5};
Line(5) = {5, 7};
Line(6) = {7, 8};
Line(7) = {8, 6};
Line Loop(8) = {6, 7, -3, 1, 2, 4, 5};
Plane Surface(9) = {8};

// Define physical groups for assigning B.C.
// 0 is for Dirichlet boundary elements
// 1 for Neumann
// 2 for iPBS
Physical Line(0) = {6,7,5};
Physical Line(1) = {4};
Physical Line(2) = {2, 1};


// At least one physical surface (or physical volume
// is needed by DUNE gmshreader (otherwise all elements
// are boundaries)
Physical Surface(13) = {9};
