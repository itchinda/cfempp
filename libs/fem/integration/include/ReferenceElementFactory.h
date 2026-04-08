/**
 * @file ReferenceElementFactory.h
 * @brief Factory class for managing singleton instances of reference finite elements.
 */

#ifndef CFEM_REFERENCE_ELEMENT_FACTORY_H
#define CFEM_REFERENCE_ELEMENT_FACTORY_H

#include "ReferenceElement.h"
#include "cfem_topology.h" // Pour CellType
#include <map>
#include <memory>

namespace cfem {

/**
 * @class ReferenceElementFactory
 * @brief Provides access to master ReferenceElement instances.
 * This factory implements a Lazy Initialization pattern (Thread-safe Singleton) 
 * to provide shape functions, gradients, and nodal coordinates for various 
 * element types and approximation orders.
 * @note This class decouples the mesh geometry (CellType) from the 
 * functional approximation space (Order). For example, a TRI3 cell can be 
 * associated with a P0, P1, or P2 reference element.
 */
class ReferenceElementFactory {
public:
    /**
     * @brief Retrieves a reference to the unique master element.
     *        The first call for a specific (type, order) pair will instantiate the 
     * element. Subsequent calls return the same static instance.
     * @param type The geometric type of the cell (e.g., TRI3, HEX8).
     * @param order The desired polynomial approximation order:
     * - 0: Constant (DG P0)
     * - 1: Linear (P1/Q1)
     * - 2: Quadratic (P2/Q2)
     * @return const ReferenceElement& A persistent reference to the element.
     * @throw std::invalid_argument If the combination of type and order 
     * is not supported or implemented.
     * * @code
     * // Example: Get a linear triangle element
     * const auto& refElem = ReferenceElementFactory::getReferenceElement(CellType::TRI3, 1);
     * @endcode
     */
    static const ReferenceElement& getReferenceElement(CellType type, int order);

    // No instances should be created - this is a static provider
    ReferenceElementFactory() = delete;
};

} // namespace cfem

#endif