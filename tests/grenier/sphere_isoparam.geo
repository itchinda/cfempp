// Set the geometry kernel to OpenCASCADE
SetFactory("OpenCASCADE");

// 1. Define Parameters
R = 1.0;        // Radius
mesh_size = 0.5; // Characteristic length

// 2. Create the Sphere Geometry
// Sphere(tag) = {centerX, centerY, centerZ, radius};
Sphere(1) = {0, 0, 0, R};

// 3. Set Mesh Options for Isoparametric Elements
// ElementOrder 2 = Quadratic (8-node quads / 10-node tets)
// ElementOrder 3 = Cubic, etc.
Mesh.ElementOrder = 2;

// High-order optimization (helps prevent inverted elements on curved surfaces)
Mesh.HighOrderOptimize = 1;

// 4. Define Physical Groups for your Solver
// Domain (Internal Volume)
Physical Volume("Universe", -1) = {1};

// Boundary (The Shell)
// In OpenCASCADE, the surface of Sphere(1) is typically Surface(1)
Physical Surface("Boundary", 1) = {1};

// 5. Mesh Size Control
// This sets the global mesh refinement
Mesh.CharacteristicLengthMin = mesh_size;
Mesh.CharacteristicLengthMax = mesh_size;
