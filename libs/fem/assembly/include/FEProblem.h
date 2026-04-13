//************************************************************************
// --- CFEM++ Library
// ---
// --- Copyright 2024-2026 Ismaël Tchinda Ngueyong et al.
// --- ALL RIGHTS RESERVED
// ---
// --- This software is protected by international copyright laws.
// --- Use, distribution, or modification of this software in any form,
// --- source or binary, for personal, academic, or non-commercial 
// --- purposes is permitted free of charge, provided that this 
// --- copyright notice and this permission notice appear in all 
// --- copies and supporting documentation.
// ---
// --- The authors and contributors provide this code "as is" without 
// --- any express or implied warranty. In no event shall they be 
// --- held liable for any damages arising from the use of this software.
//************************************************************************

// class FEProblem {
// public:
//     // Ajoute une inconnue au problème
//     void addUnknown(FEField& field) {
//         m_fields.push_back(&field);
//     }

//     // Calcule les offsets et prépare PETSc
//     void finalize() {
//         CfemInt currentOffset = 0;
//         for (auto* field : m_fields) {
//             m_offsets.push_back(currentOffset);
//             currentOffset += field->getSpace().getGlobalNumDofs();
//         }
//         m_globalSize = currentOffset;

//         // Initialisation de la matrice et des vecteurs avec tes wrappers RAII
//         m_A = std::make_unique<PetscMatrix>(PETSC_DECIDE, m_globalSize);
//         m_rhs = std::make_unique<PetscVector>(PETSC_DECIDE, m_globalSize);
//         m_sol = std::make_unique<PetscVector>(PETSC_DECIDE, m_globalSize);
//     }

//     // Utilitaire pour l'assemblage : convertit l'index local d'un champ en index global
//     CfemInt getGlobalIndex(int fieldIdx, CfemInt localDofIdx) const {
//         return m_offsets[fieldIdx] + localDofIdx;
//     }

//     // Après le solve, on redistribue les valeurs vers les FEFields
//     void distributeSolution() {
//         // Logique pour copier les segments de m_sol vers chaque FEField::getData()
//         // Puis appeler FEField::updateLocalBuffer()
//     }

// private:
//     std::vector<FEField*> m_fields;
//     std::vector<CfemInt> m_offsets;
//     CfemInt m_globalSize = 0;

//     std::unique_ptr<PetscMatrix> m_A;
//     std::unique_ptr<PetscVector> m_rhs;
//     std::unique_ptr<PetscVector> m_sol;
// };