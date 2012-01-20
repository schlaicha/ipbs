// Gmsh project created on Tue Mar  1 19:33:33 2011

radius = 5;
box_size = 50;
refine_sphere = 0.2;
refine_box = 5;
refine_length = 0;   // this is the length on symmetry axis
                    // where we want better refinement

// Define the geometry for 2d-sphere (iPBS)

// origin
Point(0) = {0, 0, 0};
Point(1) = {radius, 0, 0, refine_sphere};
Point(2) = {-radius, 0, 0, refine_sphere};
Point(3) = {0, radius, 0, refine_sphere};
Point(4) = {0, -radius, 0 , refine_sphere};
Point(5) = {0, 0, radius , refine_sphere};
Point(6) = {0, 0, -radius , refine_sphere};
// Box
Point(7) = {-box_size/2, -box_size/2, -box_size/2, refine_box};
Point(8) = {-box_size/2, box_size/2, -box_size/2, refine_box};
Point(9) = {box_size/2, -box_size/2, -box_size/2, refine_box};
Point(10) = {box_size/2, box_size/2, -box_size/2, refine_box};
Point(11) = {-box_size/2, -box_size/2, box_size/2, refine_box};
Point(12) = {-box_size/2, box_size/2, box_size/2, refine_box};
Point(13) = {box_size/2, -box_size/2, box_size/2, refine_box};
Point(14) = {box_size/2, box_size/2, box_size/2, refine_box};

// Inner sphere
Circle(1) = {3, 0, 6};
Circle(2) = {6, 0, 4};
Circle(3) = {4, 0, 5};
Circle(4) = {5, 0, 3};
Circle(5) = {3, 0, 1};
Circle(6) = {1, 0, 4};
Circle(7) = {4, 0, 2};
Circle(8) = {2, 0, 3};
Circle(9) = {6, 0, 2};
Circle(10) = {2, 0, 5};
Circle(11) = {5, 0, 1};
Circle(12) = {1, 0, 6};

Line(13) = {8, 12};
Line(14) = {12, 14};
Line(15) = {14, 10};
Line(16) = {10, 8};
Line(17) = {7, 11};
Line(18) = {11, 13};
Line(19) = {13, 9};
Line(20) = {9, 7};
Line(21) = {7, 8};
Line(22) = {11, 12};
Line(23) = {13, 14};
Line(24) = {9, 10};

// Box surface
Line Loop(25) = {14, 15, 16, 13};
Plane Surface(26) = {25};
Plane Surface(27) = {25};
Line Loop(28) = {22, -13, -21, 17};
Plane Surface(29) = {28};
Line Loop(30) = {18, 23, -14, -22};
Plane Surface(31) = {30};
Line Loop(32) = {23, 15, -24, -19};
Plane Surface(33) = {32};
Line Loop(34) = {24, 16, -21, -20};
Plane Surface(35) = {34};
Line Loop(36) = {19, 20, 17, 18};
Plane Surface(37) = {36};

// Sphere Surface
Line Loop(38) = {5, -11, 4};
Ruled Surface(39) = {38};
Line Loop(40) = {1, -12, -5};
Ruled Surface(41) = {40};
Line Loop(42) = {1, 9, 8};
Ruled Surface(43) = {42};
Line Loop(44) = {8, -4, -10};
Ruled Surface(45) = {44};
Line Loop(46) = {10, -3, 7};
Ruled Surface(47) = {46};
Line Loop(48) = {7, -9, 2};
Ruled Surface(49) = {48};
Line Loop(50) = {2, -6, 12};
Ruled Surface(51) = {50};
Line Loop(52) = {3, 11, 6};
Ruled Surface(53) = {52};

// Create the volume
Surface Loop(54) = {33, 31, 37, 35, 29, 26};
Surface Loop(55) = {39, 41, 43, 49, 47, 45, 53, 51};
Volume(56) = {54, 55};

// Define physical groups for assigning B.C.
// 0 is for Dirichlet boundary elements
// 1 for Neumann
// 2 for iPBS
// Sphere
Physical Surface(2) = {45, 39, 51, 49, 47, 53, 41, 43};
// Box
Physical Surface(1) = {31, 29, 26, 37, 35, 33};

// At least one physical surface (or physical volume
// is needed by DUNE gmshreader (otherwise all elements
// are boundaries)
Physical Volume(59) = {56};
