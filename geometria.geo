// Gmsh project created on Tue Oct 06 10:52:26 2020
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {8, 0, 0, 1.0};
//+
Point(3) = {8, 2, 0, 1.0};
//+
Point(4) = {0, 2, 0, 1.0};
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Circle(5) = {4, 1, 0, 0.5, 0, 2*Pi};
//+
Line Loop(1) = {4, 1, 2, 3};
//+
Line Loop(2) = {5};
//+
Plane Surface(1) = {1, 2};
