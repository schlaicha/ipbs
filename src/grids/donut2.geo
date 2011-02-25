// Origin
Point(1) = {0, 0, 0, 0};
// Inner Circle Points
Point(2) = {1, 0, 0, 0};
Point(3) = {0, 1, 0, 0};
Point(4) = {-1, 0, 0, 0};
Point(5) = {0, -1, 0, 0};
// Outer Circle Points
Point(6) = {5, 0, 0, 0};
Point(7) = {-5, 0, 0, 0};
Point(8) = {0, 5, 0, 0};
Point(9) = {0, -5, 0, 0};

// Define Outer Circle
Circle(1) = {8, 1, 6};
Circle(2) = {6, 1, 9};
Circle(3) = {9, 1, 7};
Circle(4) = {7, 1, 8};
// Define Inner Circle
Circle(5) = {3, 1, 2};
Circle(6) = {2, 1, 5};
Circle(7) = {5, 1, 4};
Circle(8) = {4, 1, 3};

// Define Surface
Line Loop(9) = {1, 2, 3, 4};
Line Loop(10) = {5, 6, 7, 8};
Plane Surface(11) = {9, 10};
Transfinite Line {5, 8, 7, 6} = 32 Using Progression 1;
Transfinite Line {1, 2, 3, 4} = 16 Using Progression 1;
