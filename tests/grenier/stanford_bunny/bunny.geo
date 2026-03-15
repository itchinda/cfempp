Merge "Bunny.stl"; // Importe la géométrie

// Désactive l'influence du maillage STL importé sur l'intérieur
Mesh.CharacteristicLengthExtendFromBoundary = 0;

// Définit la taille de maille souhaitée (ex: 2.0 unités)
Mesh.CharacteristicLengthMax = 2.0;
Mesh.CharacteristicLengthMin = 0.5;


Surface Loop(1) = {1}; // Définit l'enveloppe (vérifiez l'ID de la surface)
Volume(1) = {1}; // Crée le volume interne pour le maillage 3D
Physical Volume("Universe") = {1};