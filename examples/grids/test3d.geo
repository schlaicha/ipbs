// Gmsh project created on Tue Mar  1 19:33:33 2011

radius = 5;
box_x = 100;
box_y = 100;
box_z = 100;
refine_sphere = 2;
refine_box = 10;

// reference sphere
pos0_x = 0;
pos0_y = 0;
pos0_z = 0;

pos1_x = +15;
pos1_y = 0;
pos1_z = 0;

pos2_x = -15;
pos2_y = 0;
pos2_z = 0;

pos3_x = -30;
pos3_y = 0;
pos3_z = 0;

pos4_x = 30;
pos4_y = 0;
pos4_z = 0;

pos5_x = 0;
pos5_y = 15;
pos5_z = 0;

pos6_x = 0;
pos6_y = 30;
pos6_z = 0;

pos7_x = 0;
pos7_y = -15;
pos7_z = 0;

pos8_x = 0;
pos8_y = -30;
pos8_z = 0;

pos9_x = 0;
pos9_y = 0;
pos9_z = 15;

pos10_x = 0;
pos10_y = 0;
pos10_z = 30;

pos11_x = 0;
pos11_y = 0;
pos11_z = 30;

pos12_x = 0;
pos12_y = 0;
pos12_z = -15;

pos13_x = 0;
pos13_y = 0;
pos13_z = -30;

// Box
Point(7) = {- box_x/2,    -  box_y/2,   -box_z/2, refine_box};
Point(8) = {- box_x/2,       box_y/2,   -box_z/2, refine_box};
Point(9) = {  box_x/2,     -  box_y/2,  -box_z/2, refine_box};
Point(10) = { box_x/2,       box_y/2,   -box_z/2, refine_box};
Point(11) = {-box_x/2,   -  box_y/2,     box_z/2, refine_box};
Point(12) = {-box_x/2,      box_y/2,     box_z/2, refine_box};
Point(13) = { box_x/2,    -  box_y/2,    box_z/2, refine_box};
Point(14) = { box_x/2,       box_y/2,    box_z/2, refine_box};
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

// sphere 
Point(0) = {pos0_x,           pos0_y,           pos0_z};
Point(1) = {pos0_x + radius,  pos0_y,           pos0_z,           refine_sphere};
Point(2) = {pos0_x - radius,  pos0_y,           pos0_z,           refine_sphere};
Point(3) = {pos0_x,           pos0_y + radius,  pos0_z,           refine_sphere};
Point(4) = {pos0_x,           pos0_y + -radius, pos0_z,           refine_sphere};
Point(5) = {pos0_x,           pos0_y,           pos0_z + radius,  refine_sphere};
Point(6) = {pos0_x,           pos0_y,           pos0_z - radius,  refine_sphere};
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

// Sphere 1 Surface
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



// now create other spheres
part_1[] = Translate {pos1_x - pos0_x, pos1_y - pos0_y, pos1_z - pos0_z} {
  Duplicata { Surface{39, 53, 47, 45, 43, 51, 49, 41}; }
};
part_2[] = Translate {pos2_x - pos0_x, pos2_y - pos0_y, pos2_z - pos0_z} {
  Duplicata { Surface{39, 53, 47, 45, 43, 51, 49, 41}; }
};
part_3[] = Translate {pos3_x - pos0_x, pos3_y - pos0_y, pos3_z - pos0_z} {
  Duplicata { Surface{39, 53, 47, 45, 43, 51, 49, 41}; }
};
part_4[] = Translate {pos4_x - pos0_x, pos4_y - pos0_y, pos4_z - pos0_z} {
  Duplicata { Surface{39, 53, 47, 45, 43, 51, 49, 41}; }
};
part_5[] = Translate {pos5_x - pos0_x, pos5_y - pos0_y, pos5_z - pos0_z} {
  Duplicata { Surface{39, 53, 47, 45, 43, 51, 49, 41}; }
};
part_6[] = Translate {pos6_x - pos0_x, pos6_y - pos0_y, pos6_z - pos0_z} {
  Duplicata { Surface{39, 53, 47, 45, 43, 51, 49, 41}; }
};
part_7[] = Translate {pos7_x - pos0_x, pos7_y - pos0_y, pos7_z - pos0_z} {
  Duplicata { Surface{39, 53, 47, 45, 43, 51, 49, 41}; }
};
part_8[] = Translate {pos8_x - pos0_x, pos8_y - pos0_y, pos8_z - pos0_z} {
  Duplicata { Surface{39, 53, 47, 45, 43, 51, 49, 41}; }
};
part_9[] = Translate {pos9_x - pos0_x, pos9_y - pos0_y, pos9_z - pos0_z} {
  Duplicata { Surface{39, 53, 47, 45, 43, 51, 49, 41}; }
};
part_10[] = Translate {pos10_x - pos0_x, pos10_y - pos0_y, pos10_z - pos0_z} {
  Duplicata { Surface{39, 53, 47, 45, 43, 51, 49, 41}; }
};
part_11[] = Translate {pos11_x - pos0_x, pos11_y - pos0_y, pos11_z - pos0_z} {
  Duplicata { Surface{39, 53, 47, 45, 43, 51, 49, 41}; }
};
part_12[] = Translate {pos12_x - pos0_x, pos12_y - pos0_y, pos12_z - pos0_z} {
  Duplicata { Surface{39, 53, 47, 45, 43, 51, 49, 41}; }
};
part_13[] = Translate {pos13_x - pos0_x, pos13_y - pos0_y, pos13_z - pos0_z} {
  Duplicata { Surface{39, 53, 47, 45, 43, 51, 49, 41}; }
};

// Create the volume
Surface Loop(54) = {33, 31, 37, 35, 29, 26};
Surface Loop(55) = {39, 41, 43, 49, 47, 45, 53, 51};
Surface Loop(56) = part_1[];
Surface Loop(57) = part_2[];
Surface Loop(58) = part_3[];
Surface Loop(59) = part_4[];
Surface Loop(60) = part_5[];
Surface Loop(61) = part_6[];
Surface Loop(62) = part_7[];
Surface Loop(63) = part_8[];
Surface Loop(64) = part_9[];
Surface Loop(65) = part_10[];
Surface Loop(66) = part_11[];
Surface Loop(67) = part_12[];
Surface Loop(68) = part_13[];
Volume(100) = {54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 68};

// Define physical groups for assigning B.C.

// Box
Physical Surface(0) = {31, 29, 26, 37, 35, 33};

// Sphere
//Physical Surface(1) = {45, 39, 51, 49, 47, 53, 41, 43};
Physical Surface(1) = {45, 39, 51, 49, 47, 53, 41, 43, part_1[], part_2[], part_3[], part_4[], part_5[], part_6[], part_7[], part_8[], part_9[], part_10[], part_11[], part_12[], part_13[]};
//Physical Surface(1) = {45, 39, 51, 49, 47, 53, 41, 43, part_1[], part_2[], part_3[], part_4[]};

// At least one physical surface (or physical volume
// is needed by DUNE gmshreader (otherwise all elements
// are boundaries)
Physical Volume(59) = {100};
