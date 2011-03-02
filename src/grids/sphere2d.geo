// Gmsh project created on Tue Mar  1 19:33:33 2011

// Define the geometry for 2d-sphere (iPBS)

Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {-1, 0, 0};
Point(4) = {-5, 0, 0};
Point(5) = {5, 0, 0};
Point(6) = {0, 1, 0};
Point(7) = {-5, 5, 0};
Point(8) = {5, 5, 0};
Line(1) = {4, 7};
Line(2) = {7, 8};
Line(3) = {8, 5};
Point(9) = {1.2, 0, 0};
Point(10) = {-1.2, 0, 0};
Line(4) = {2, 9};
Line(5) = {9, 5};
Line(6) = {3, 10};
Line(7) = {10, 4};
Circle(8) = {3, 1, 6};
Circle(9) = {6, 1, 2};
Transfinite Line {9, 8} = 32 Using Progression 1;
Transfinite Line {4, 6} = 4 Using Progression 1;
Line Loop(10) = {2, 3, -5, -4, -9, -8, 6, 7, 1};
Plane Surface(11) = {10};

// Define physical groups for assigning B.C.
// 0 is for Dirichlet boundary elements
// 1 for Neumann
// 2 for iPBS

Physical Point(0) = {10, 4, 9, 5};
Physical Point(1) = {2, 6, 3};
Physical Point(2) = {4, 7, 8, 5};
Physical Line(0) = {1, 2, 3};
Physical Line(1) = {5, 7, 6, 4};
Physical Line(2) = {9, 8};
