// Define setup variables
radius = 1;
distance = 2;
height = 10;
width = 9 + distance;
lc_sphere = .01;
//lc_center = @lc_center@;
lc_center = distance/radius*0.1;
lc_boarder = 1;

// ========== Setup Geometry ===========

// Define points
Point(0) = {0, 0, 0, lc_center};
Point(1) = {-distance+radius, 0, 0, lc_sphere};
Point(2) = {-distance, radius, 0, lc_sphere};
Point(3) = {-distance-radius, 0, 0, lc_sphere};
Point(4) = {-width, 0, 0, lc_boarder};
Point(5) = {-width, height, 0, lc_boarder};
Point(6) = {0, height, 0, lc_boarder};
Point(7) = {-distance, 0, 0, lc_boarder};

// Define lines
Line(1) = {0, 1};
Line(2) = {3, 4};
Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {6, 0};
Circle(6) = {1, 7, 2};
Circle(7) = {2, 7, 3};

// Define surface
Line Loop(8) = {5, 1, 6, 7, 2, 3, 4};
Plane Surface(9) = {8};

// ========== Physical Groups (Boundary Conditions) ===========

Physical Line(1) = {5, 4, 3, 2, 1};
Physical Line(2) = {7, 6};
Physical Surface(10) = {9};
