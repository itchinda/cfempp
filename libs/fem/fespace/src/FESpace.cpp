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

#include "FESpace.h"
#include "Mesh.h"
#include "ReferenceElements.h"

namespace cfem
{
    /**
     * @brief Registers a finite element field (e.g., Temperature, Velocity) associated with this space.
     * * This method guarantees that:
     * 1. The field's name adheres to the proper naming policy rules.
     * 2. There are no global naming collisions across the mesh.
     * 3. The space keeps track of the field to notify it during mesh/topology updates.
     * * @param field Shared pointer to the FEField to be registered.
     */
    void FESpace::registerField(std::shared_ptr<FESpaceObserver> field)
    {
        // Enforce naming conventions
        const std::string& name = field->getFESpaceObserverName();
        NamingPolicy::validate(name);

        // Ensure the field name is unique within the context of the mesh
        if (m_mesh->nameAlreadyUsed(name))
        {
            CFEM_ERROR << "Global naming collision: Field with name '" << name << "' already exists.";
        }
        
        // Store a weak reference to the field to avoid circular dependency memory leaks
        m_observers.push_back(field);
        m_mesh->registerFieldName(name);
    }

    void FESpace::validate() const {
        if (m_mesh->getNumCells() == 0) {
            CFEM_ERROR << "FESpace Error: Mesh is empty.";
        }
        if (m_numComponents < 1) {
            CFEM_ERROR << "FESpace Error: Field dimension must be >= 1.";
        }
    }

    /**
     * @brief Entry point for the update cascade triggered by a mesh/topology modification.
     * * When the mesh changes (e.g., due to refinement or deformation), the continuous 
     * or discontinuous space must recalculate its global numbering. Consequently, all 
     * fields that rely on this space must resize their algebraic data structures (like PETSc vectors).
     * * The cascade works as follows:
     * 1. Validate the current state of the space.
     * 2. Re-enumerate the Degrees of Freedom (DOFs) via the specialized subclass implementation.
     * 3. Notify all living associated FEFields to reallocate their internal PETSc vectors.
     */
    void FESpace::updateAfterMeshChange() {
        this->validate();
        
        // Let the specific implementation (e.g., LagrangeSpace) handle DOF numbering
        this->enumerateDofs();

        // Notify all registered fields (Temperature, Velocity, etc.)
        for (auto it = m_observers.begin(); it != m_observers.end(); ) {
            if (auto obs = it->lock()) {
                obs->updateOnSpaceChange(); 
                ++it;
            } else {
                // The observer (FEField) has been destroyed elsewhere in the program.
                // Erase the expired weak pointer to prevent memory overhead.
                it = m_observers.erase(it);
            }
        }
    }

    /**
     * @brief Initializes all supported reference elements based on the mesh dimension.
     * Covers both Linear (P1/Q1) and Quadratic (P2/Q2) elements across 1D, 2D, and 3D.
     */
    void FESpace::initReferenceElements() {
        m_refElements.clear();
        CfemInt dim = m_mesh->getTopologicalDimension();

        // Utilisation de la même logique de filtrage que LagrangeSpace
        // Mais adaptée pour supporter les éléments DG (souvent P0 ou P1 discontinu)
        if (dim == 1) {
            m_refElements[CellType::LINE2] = std::make_unique<P1ReferenceSegment>();
            m_refElements[CellType::LINE3] = std::make_unique<P2ReferenceSegment>();
        } else if (dim == 2) {
            m_refElements[CellType::TRI3]  = std::make_unique<P1ReferenceTriangle>();
            m_refElements[CellType::QUAD4] = std::make_unique<Q1ReferenceQuad>();

            m_refElements[CellType::TRI6]  = std::make_unique<P2ReferenceTriangle>();
            m_refElements[CellType::QUAD9] = std::make_unique<Q2ReferenceQuad>();
        } else if (dim == 3) {
            m_refElements[CellType::TET4]    = std::make_unique<P1ReferenceTetrahedron>();
            m_refElements[CellType::PRISM6]  = std::make_unique<P1ReferencePrism>();
            m_refElements[CellType::HEX8]    = std::make_unique<Q1ReferenceHexahedron>();

            m_refElements[CellType::TET10]    = std::make_unique<P2ReferenceTetrahedron>();
            m_refElements[CellType::PRISM18]  = std::make_unique<P2ReferencePrism>();
            m_refElements[CellType::HEX27]    = std::make_unique<Q2ReferenceHexahedron>();
        }
    }


    const ReferenceElement& FESpace::getRefElement(CellType type) const {
        auto it = m_refElements.find(type);
        if (it == m_refElements.end()) {
            CFEM_ERROR << "FESpace: Element type not initialized for type " << static_cast<int>(type);
        }
        return *it->second;
    }

} // namespace cfem