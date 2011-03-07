// Define setup variables
radius = 1;
distance = 3;
world_size = 8;

// ========== Setuo Geometry ===========

// Define points
Point(0) = {0, 0, 0, 0.2}; // origin
Point(1) = {distance, 0, 0, 1.0}; // center of right sphere
Point(2) = {-distance, 0, 0, 1.0}; // center of left sphere
Point(3) = {-world_size, 0, 0, 0.3}; // lower left world coordinate
Point(4) = {world_size, 0, 0, 0.3}; // lower right world coordinate
Point(5) = {-world_size, world_size, 0, 1.0}; // upper right world coordinate
Point(6) = {world_size, world_size, 0, 1.0}; // upper right world coordinate
Point(7) = {distance+radius, 0, 0, 0.2}; // right boarder of right sphere
Point(8) = {distance-radius, 0, 0, 0.2}; // left boarder of right sphere
Point(9) = {distance, radius, 0, 0.2}; // top boarder of right sphere
Point(10) = {-distance+radius, 0, 0, 0.2}; // right boarder of left sphere
Point(11) = {-distance-radius, 0, 0, 0.2}; // left boarder of left sphere
Point(12) = {-distance, radius, 0, 0.2}; // top boarder of left sphere

// Define lines
Line(1) = {3, 5}; // left boerder
Line(2) = {5, 6}; // upper boarder
Line(3) = {6, 4}; // left boerder
Line(4) = {4, 7}; // lower right boarder
Line(5) = {3, 11}; // lower left boarder
Line(6) = {8, 0}; //right middle lower boarder
Line(7) = {0, 10}; // left middle lower boarder
Circle(8) = {7, 1, 9}; // right sphere (right part)
Circle(9) = {9, 1, 8}; // right sphere (left part)
Circle(10) = {10, 2, 12}; // left sphere (right part)
Circle(11) = {12, 2, 11}; // left sphere (left part)

// Define surface
Line Loop(12) = {2, 3, 4, 8, 9, 6, 7, 10, 11, -5, 1};
Plane Surface(13) = {12};


// ========== Refinement ===========

// Refine the line interconnecting the spheres
Transfinite Line {7, 6} = 12 Using Progression 1;
// Refine the spheres
Transfinite Line {9, 8, 10, 11} = 18 Using Progression 1;
// Refine lines from spheres to outer boundary
Transfinite Line {4, 5} = 12 Using Progression 1;


// ========== Physical Groups (Boundary Conditions) ===========

// Outer boundaries, which are set to Dirichlet B.C.
Physical Line(0) = {1, 2, 3};
// Lower boundaries, wich are set to Neumann B.C.
Physical Line(1) = {4, 6, 7, 5};
// Spheres which are set to IPBS B.C.
Physical Line(2) = {8, 9, 10, 11};

// Finally we also need a physical surface
Physical Surface(14) = {13};
