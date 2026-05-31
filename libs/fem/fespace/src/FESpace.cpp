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

namespace cfem
{
    /**
     * @brief Registers a finite element field (e.g., Temperature, Velocity) associated with this space.
     * 
     * This method guarantees that:
     * 1. The field's name adheres to the proper naming policy rules.
     * 2. There are no global naming collisions across the mesh.
     * 3. The space keeps track of the field to notify it during mesh/topology updates.
     * 
     * @param field Shared pointer to the FEField to be registered.
     */
    void FESpace::registerField(std::shared_ptr<FESpaceObserver> field)
    {
        // Enforce naming conventions
        const std::string& name = field->getFESpaceObserverName();
        NamingPolicy::validate(name);

        // Ensure the field name is unique within the context of the mesh
        if (m_mesh->nameAlreadyUsed(name))
        {
            CFEM_THROW("Global naming collision: Field with name '" << name << "' already exists.");
        }
        
        // Store a weak reference to the field to avoid circular dependency memory leaks
        m_observers.push_back(field);
        m_mesh->registerFieldName(name);
    }

    void FESpace::validate() const {
        if (m_mesh->getNumCells() == 0) {
            CFEM_THROW("FESpace Error: Mesh is empty.");
        }
        if (m_numComponents < 1) {
            CFEM_THROW("FESpace Error: Field dimension must be >= 1.");
        }
    }

    /**
     * @brief Entry point for the update cascade triggered by a mesh/topology modification.
     * 
     * When the mesh changes (e.g., due to refinement or deformation), the continuous 
     * or discontinuous space must recalculate its global numbering. Consequently, all 
     * fields that rely on this space must resize their algebraic data structures (like PETSc vectors).
     * 
     * The cascade works as follows:
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

    const ReferenceElement& FESpace::getReferenceElement(CfemInt cellId) const{
        return ReferenceElementFactory::getReferenceElement(m_mesh->getCell(cellId).type, getOrderInt());
    }

} // namespace cfem