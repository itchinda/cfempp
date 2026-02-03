h=1;

// 1. Points (Lower face)
Point(1) = {0, 0, 0, h};
Point(2) = {1, 0, 0, h};
Point(3) = {1, 1, 0, h};
Point(4) = {0, 1, 0, h};

// 2. Points (Upper face)
Point(5) = {0, 0, 1, h};
Point(6) = {1, 0, 1, h};
Point(7) = {1, 1, 1, h};
Point(8) = {0, 1, 1, h};

// 3. Lines (Bottom face)
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// 4. Lines (Top face)
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

// 5. Lines (Verticals)
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// 6. Surfaces
Curve Loop(1) = {1, 2, 3, 4};        Plane Surface(1) = {1}; // Bottom
Curve Loop(2) = {5, 6, 7, 8};        Plane Surface(2) = {2}; // Top
Curve Loop(3) = {1, 10, -5, -9};     Plane Surface(3) = {3}; // Front
Curve Loop(4) = {2, 11, -6, -10};    Plane Surface(4) = {4}; // Right
Curve Loop(5) = {3, 12, -7, -11};    Plane Surface(5) = {5}; // Back
Curve Loop(6) = {4, 9, -8, -12};     Plane Surface(6) = {6}; // Left

// 7. Volume
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Mesh.ElementOrder = 2;

// 8. Force Hexahedral Meshing
Transfinite Curve {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12} = 2; 
Transfinite Surface {1, 2, 3, 4, 5, 6};
Transfinite Volume{1};
//Recombine Surface {1, 2, 3, 4, 5, 6};

Physical Surface("Bottom", 1)={1};
Physical Surface("Top", 2)={2};
Physical Surface("Front", 3)={3};
Physical Surface("Right", 4)={4};
Physical Surface("Back", 5)={5};
Physical Surface("Left", 6)={6};

Physical Volume("Universe", 1)= {1};

Physical Surface("Bord", 7) = {1,2};
