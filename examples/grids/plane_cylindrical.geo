// Gmsh project created on Tue Aug 23 10:26:33 2011
// Creates a charged wall for iPBS using spherical symmetry
// refinement is done via the characteristic length

box_length = 10;
box_width = 100;
left_refinement = .05;
right_refinement = 2;
center_refinement = 1;

// Define the geometry for the wall (iPBS)

Point(0) = {box_width/2, box_length/2, 0, center_refinement};
Point(1) = {0, 0, 0, left_refinement};
Point(2) = {0, box_length, 0, left_refinement};
Point(3) = {box_length, 0, 0, right_refinement};
Point(4) = {box_length, box_length, 0, right_refinement};
Line(1) = {1, 2};
Line(2) = {3, 4};
Line(3) = {1, 3};
Line(4) = {2, 4};
Line Loop(5) = {4, -2, -3, 1};
Plane Surface(6) = {5};

// Define physical groups for assigning B.C.
// 0 is for Dirichlet boundary elements
// 1 for Neumann
// 2 for iPBS

// upper und lower wall are zero flux due to symmetry
Physical Line(1) = {4, 3};
Physical Line(2) = {1};
Physical Line(0) = {2};

// At least one physical surface (or physical volume
// is needed by DUNE gmshreader (otherwise all elements
// are boundaries)

Physical Surface(7) = {6};
