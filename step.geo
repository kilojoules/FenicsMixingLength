C = 1;
h3 = .001 / C;
h2 = .01 / C;
h1 = .1 / C;
//+
Point(1) = {0, 0, 0, h2};
//+
Point(2) = {0, 1, 0, h2};
//+
Point(3) = {-0.2, 1, 0, h2};
//+
Point(4) = {-0.2, 2.5, 0, h1};
//+
Point(5) = {10, 2.5, 0, h1};
//+
Point(6) = {10, 0, 0, h1};
Point(7) = {1, 0, 0, h2};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
Line(7) = {7, 1};
//+
Line Loop(1) = {4, 5, 6, 7, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Field[1] = Box;
//+
Field[1].VIn = h3;
//+
Field[1].VOut = h1;
//+
Field[1].XMax = 0.02;
//+
Field[1].XMin = -3.2;
//+
Field[1].YMax = 1.03;
//+
Field[1].YMin = -0.2;
//+
Background Field = 1;
