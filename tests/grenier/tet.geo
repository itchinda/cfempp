// Single Tetrahedron Debug Mesh
SetFactory("OpenCASCADE");

h=2;
// 1. Define the 4 vertices (a simple unit tet)
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {0, 1, 0, h};
Point(4) = {0, 0, 1, h};

// 2. Build the geometry
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line(4) = {1, 4};
Line(5) = {2, 4};
Line(6) = {3, 4};

Line Loop(1) = {1, 2, 3};      // Base
Line Loop(2) = {1, 5, -4};     // Side 1
Line Loop(3) = {2, 6, -5};     // Side 2
Line Loop(4) = {3, 4, -6};     // Side 3

Surface(1) = {1}; Surface(2) = {2}; Surface(3) = {3}; Surface(4) = {4};
Surface Loop(1) = {1, 2, 3, 4};
Volume(1) = {1};

// 3. Define Physical Groups for your solver
Physical Volume("Universe", -1) = {1};
Physical Surface("Boundary", -2) = {1, 2, 3, 4};

// 4. Mesh Order
// Change to 2 for Tet10, 1 for Tet4
Mesh.ElementOrder = 2;
