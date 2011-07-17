// Gmsh project created on Tue Jun 07 12:07:33 2011
// Creates the iPBS single sphere with a symmetric mesh
// refinement is done via the characteristic length

inner_radius = 5;
position = 0;  // Center position of the sphere
outer_radius = 50;
outer_refinement = 3;
inner_refinement = 0.3;

// Define the geometry for 2d-sphere (iPBS)

Point(0) = {0, 0, 0, 1.0};
Point(1) = {position, 0, 0, 1.0};
Point(2) = {position, inner_radius, 0, inner_refinement};
Point(3) = {-inner_radius+position, 0, 0, inner_refinement};
Point(4) = {inner_radius+position, 0, 0, inner_refinement};
Point(5) = {-outer_radius, 0, 0, outer_refinement};
Point(6) = {outer_radius, 0, 0, outer_refinement};
Point(7) = {0, outer_radius, 0, outer_refinement};
Circle(8) = {2,1,3};
Circle(9) = {2,1,4};
Circle(10) = {6, 0, 7};
Circle(11) = {7, 0, 5};
Line(12) = {4, 6};
Line(13) = {3, 5};
Line Loop(14) = {10, 11, -13, -8, 9, 12};
Plane Surface(15) = {14};

// Define physical groups for assigning B.C.
// 0 is for Dirichlet boundary elements
// 1 for Neumann
// 2 for iPBS

Physical Line(2) = {9, 8};
Physical Line(1) = {13, 12};
// Dirichlet on outer bpundary
Physical Line(0) = {10,11};

// At least one physical surface (or physical volume
// is needed by DUNE gmshreader (otherwise all elements
// are boundaries)
Physical Surface(16) = {15};
