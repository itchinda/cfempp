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

/**
 * @file ReferenceElementFactory.h
 * @brief Factory class for managing singleton instances of reference finite elements.
 */

#ifndef CFEM_REFERENCE_ELEMENT_FACTORY_H
#define CFEM_REFERENCE_ELEMENT_FACTORY_H

#include "cfem_topology.h" // For CellType definition
#include "ReferenceElement.h"

namespace cfem {

// Utilitaire interne pour garantir une instance unique par type
template<typename T>
const T& singleton() {
    static T instance;
    return instance;
}

/**
 * @class ReferenceElementFactory
 * @brief Provides ultra-fast O(1) access to master ReferenceElement instances.
 * * This factory uses a pre-initialized lookup table to provide shape functions
 * and mapping data. It decouples the geometric cell type from the 
 * functional approximation order.
 */
class ReferenceElementFactory {
public:
    /**
     * @brief Retrieves a reference to the unique master element.
     * * @param type The geometric type of the cell (e.g., TRI3, HEX8).
     * @param order The polynomial approximation order (0, 1, 2 or Order::P1Mini).
     * @return const ReferenceElement& Persistent reference to the requested element.
     * @throw std::invalid_argument If the combination (type, order) is not implemented.
     */
    static const ReferenceElement& getReferenceElement(CellType type, CfemInt order);

    // Static provider only
    ReferenceElementFactory() = delete;

private:
    /**
     * @brief Internal helper to map the order (0, 1, 2, 101) to a dense array index [0..3].
     */
    static CfemInt orderToIndex(CfemInt order);
};

} // namespace cfem

#endif // CFEM_REFERENCE_ELEMENT_FACTORY_H