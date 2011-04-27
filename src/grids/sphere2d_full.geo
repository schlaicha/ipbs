// Gmsh project created on Tue Mar  1 19:33:33 2011

radius = 5;
box_size = 50;
refine_length = 0;   // this is the length on symmetry axis
                    // where we want better refinement

// Define the geometry for 2d-sphere (iPBS)

Point(1) = {0, 0, 0};
Point(2) = {radius, 0, 0};
Point(3) = {-radius, 0, 0};
Point(4) = {-box_size, -box_size, 0};
Point(5) = {box_size, -box_size, 0};
Point(6) = {0, radius, 0};
Point(10) = {0, -radius, 0, 1.0};
Point(7) = {-box_size, box_size, 0};
Point(8) = {box_size, box_size, 0};

Line(1) = {4, 7};
Line(2) = {7, 8};
Line(3) = {8, 5};
Line(4) = {5, 4};

Circle(5) = {6, 1, 2};
Circle(6) = {2, 1, 10};
Circle(7) = {10, 1, 3};
Circle(8) = {3, 1, 6};

Line Loop(9) = {3, 4, 1, 2};
Line Loop(10) = {6, 7, 8, 5};
Plane Surface(11) = {9, 10};

// Refinement
Transfinite Line {6, 5, 8, 7} = 32 Using Progression 1;

// Define physical groups for assigning B.C.
// 0 is for Dirichlet boundary elements
// 1 for Neumann
// 2 for iPBS
Physical Line(0) = {4, 3, 2, 1};
Physical Line(2) = {8, 5, 6, 7};

// At least one physical surface (or physical volume
// is needed by DUNE gmshreader (otherwise all elements
// are boundaries)
Physical Surface(12) = {11};
