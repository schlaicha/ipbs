lc = 0.1;
Point(1) = {0, 0, 0, 5.0};
Point(2) = {5, 0, 0, 5.0};
Point(3) = {5, 5, 0, 5};
Point(4) = {0, 5, 0, 5};
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};
Line Loop(5) = {3, 4, 1, 2};
Plane Surface(6) = {5};

