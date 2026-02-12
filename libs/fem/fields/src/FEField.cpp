#include "FEField.h"
#include <stdexcept>
#include <iostream>
#include "ReferenceElementFactory.h"

namespace cfem {

FEField::FEField(ConstructorToken, const std::string& name, std::shared_ptr<LagrangeSpace> space)
    : m_name(name), 
      m_space(std::move(space)), 
      m_data(PETSC_DECIDE, m_space->getGlobalNumDofs()) 
{
    // Initialisation du champ à zéro
    m_data.setAll(0.0);
}

void FEField::reallocate() {
    // 1. On récupère la nouvelle taille (après un upgrade du maillage par exemple)
    CfemInt newSize = m_space->getGlobalNumDofs();
    
    // 2. On reconstruit le vecteur PETSc. 
    // Si ton wrapper PetscVector gère bien l'opérateur d'affectation/move, 
    // l'ancien vecteur sera détruit proprement (RAII).
    m_data = PetscVector(PETSC_DECIDE, newSize);
    m_data.setAll(0.0);
    
    // 3. On invalide le cache local
    m_localValues.clear();
    
    std::cout << "[FEField] Field '" << m_name << "' reallocated to size " << newSize << std::endl;
}

void FEField::updateLocalBuffer() const {
    CfemInt totalDofs = m_space->getGlobalNumDofs();
    
    // On s'assure que le buffer CPU a la bonne taille
    if (m_localValues.size() != static_cast<size_t>(totalDofs)) {
        m_localValues.assign(totalDofs, 0.0);
    }

    // --- Logique PETSc Scatter ---
    VecScatter ctx;
    Vec localVec;
    
    // Création d'un vecteur séquentiel temporaire
    VecCreateSeq(PETSC_COMM_SELF, totalDofs, &localVec);
    VecScatterCreateToAll(m_data.getRaw(), &ctx, nullptr);

    VecScatterBegin(ctx, m_data.getRaw(), localVec, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, m_data.getRaw(), localVec, INSERT_VALUES, SCATTER_FORWARD);

    // Transfert vers std::vector
    PetscScalar* array;
    VecGetArray(localVec, &array);
    for (CfemInt i = 0; i < totalDofs; ++i) {
        m_localValues[i] = static_cast<CfemReal>(array[i]);
    }
    VecRestoreArray(localVec, &array);

    VecScatterDestroy(&ctx);
    VecDestroy(&localVec);
}

CfemReal FEField::getAtNode(CfemInt meshNodeId, int component) const {
    // On utilise l'accès par pointeur (->) au lieu de (.) car m_space est un shared_ptr
    CfemInt nodeIdx = m_space->getLocalNodeIdx(meshNodeId);
    
    if (nodeIdx == -1) {
        return 0.0;
    }
    // Logique de numérotation par BLOC
    CfemInt numNodes = m_space->getNumActiveNodes();
    CfemInt globalIdx = nodeIdx + (component * numNodes);

    // Sécurité : si le buffer n'est pas à jour
    if (globalIdx >= static_cast<CfemInt>(m_localValues.size())) {
        return 0.0; 
    }

    return m_localValues[globalIdx];
}

CfemReal FEField::getAt(CfemInt globalDofIdx, int component) const {
    // Si l'index est invalide ou négatif (cas des DOFs contraints/fixés selon ton implémentation)
    if (globalDofIdx < 0) return 0.0;
    
    // Le stockage local m_localValues est souvent entrelacé (interleaved) :
    // [u0_x, u0_y, u1_x, u1_y, ...]
    int numComponents = this->getDimension(); 
    
    return m_localValues[globalDofIdx * numComponents + component];
}

void FEField::evaluateAt(CfemInt cellId, const Vector3D& xi, std::vector<double>& outValues) const {
    int nComp = this->getDimension(); 
    outValues.assign(nComp, 0.0);

    // 1. DETERMINER L'ORDRE REEL DU CHAMP
    // On suppose que si c'est P1, on ne veut que 3 fonctions de forme.
    // On va "tricher" en récupérant l'élément de référence P1 explicitement.
    const Cell& cell = m_space->getMesh()->getCell(cellId);
    
    // Détection : est-ce un triangle ? un quad ?
    // On force un type "linéaire" pour l'évaluation si l'espace est P1
    CellType baseType = cell.type;
    if (m_space->getOrder() == Order::P1) {
        if (baseType == CellType::LINE3) baseType = CellType::LINE2;
        if (baseType == CellType::TRI6) baseType = CellType::TRI3;
        if (baseType == CellType::QUAD9) baseType = CellType::QUAD4;
        if (baseType == CellType::TET10) baseType = CellType::TET4;
        if (baseType == CellType::PRISM18) baseType = CellType::PRISM6;
        if (baseType == CellType::HEX27) baseType = CellType::HEX8;
    }

    const auto& refEl = m_space->getRefElement(baseType);
    ShapeInfo shape;
    
    refEl.evaluate(xi, ShapeEvaluationFlags::Value, shape);
    const auto& shapeValues = shape.values;
    
    // 2. INTERPOLATION STRICTE
    for (size_t i = 0; i < shapeValues.size(); ++i) {
        CfemInt mNodeId = cell.vertex_ids[i]; 
        CfemInt dofIdx = m_space->getLocalNodeIdx(mNodeId);
        
        if (dofIdx != -1) {
            CfemInt numNodes = m_space->getNumActiveNodes();
            for (int c = 0; c < nComp; ++c) {
                outValues[c] += (double)shapeValues[i] * m_localValues[dofIdx + (c * numNodes)];
            }
        }
    }
}

} // namespace cfem