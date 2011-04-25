// Gmsh project created on Tue Mar  1 19:33:33 2011

radius = 5;
box_size = 50;
refine_length = 2;   // this is the length on symmetry axis
                    // where we want better refinement

// Define the geometry for 2d-sphere (iPBS)

Point(1) = {0, 0, 0};
Point(2) = {radius, 0, 0};
Point(3) = {-radius, 0, 0};
Point(4) = {-box_size, 0, 0};
Point(5) = {box_size, 0, 0};
Point(6) = {0, radius, 0};
Point(7) = {-box_size, box_size, 0};
Point(8) = {box_size, box_size, 0};
Line(1) = {4, 7};
Line(2) = {7, 8};
Line(3) = {8, 5};
Point(9) = {radius+refine_length, 0, 0};
Point(10) = {-radius-refine_length, 0, 0};
Line(4) = {2, 9};
Line(5) = {9, 5};
Line(6) = {3, 10};
Line(7) = {10, 4};
Circle(8) = {3, 1, 6};
Circle(9) = {6, 1, 2};
Transfinite Line {9, 8} = 64 Using Progression 1;
Transfinite Line {4, 6} = 32 Using Progression 1;
Transfinite Line {7, 5} = 16 Using Progression 1;
Line Loop(10) = {2, 3, -5, -4, -9, -8, 6, 7, 1};
Plane Surface(11) = {10};

// Define physical groups for assigning B.C.
// 0 is for Dirichlet boundary elements
// 1 for Neumann
// 2 for iPBS

Physical Line(0) = {1, 2, 3};
Physical Line(1) = {5, 7, 6, 4};
Physical Line(2) = {9, 8};

// At least one physical surface (or physical volume
// is needed by DUNE gmshreader (otherwise all elements
// are boundaries)

Physical Surface(12) = {11};
