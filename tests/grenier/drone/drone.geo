// 1. Tell Gmsh to use the OpenCASCADE CAD kernel
SetFactory("OpenCASCADE");

// 2. Load your STEP file
// Replace with your actual filename
v() = ShapeFromFile("drone.STEP"); 

// 3. Optional: Force a fine mesh size to push CFEM++ to its limits
// Tweak these numbers depending on the physical scale of your CAD model

Mesh.MeshSizeMin = 0.25; 
Mesh.MeshSizeMax = 5.0;

// 4. (Crucial for FEM) Define a Physical Volume so CFEM++ can identify the domain
// This assigns the ID 1 to the generated volume
Physical Surface("Universe") = {v()};
