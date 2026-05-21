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

#ifndef CFEM_FESPACE_H
#define CFEM_FESPACE_H

#include "cfemutils.h"
#include "cfem_topology.h"
#include "ReferenceElementFactory.h"

#include <vector>
#include <memory>

namespace cfem {

    /**
     * @brief Unique identifiers for physical field types.
     */
    enum class FieldType : CfemInt { 
        Scalar,             // 1
        Vector2D,           // 2
        Vector3D,           // 3
        SymTensor2D,        // 3  [xx, yy, xy]
        SkewSymTensor3D,    // 3  [yz, xz, xy]
        Tensor2D,           // 4  [xx, xy, yx, yy]
        AxiSymTensor3D,     // 4  [xx, yy, zz, xy]
        SymTensor3D,        // 6  [xx, yy, zz, yz, xz, xy]
        Tensor3D            // 9  Full 3x3 matrix
    };

    /**
     * @brief Returns the number of independent components (DOF) for a given FieldType.
     * @param type The physical field type.
     * @return The size of the field.
     */
    [[nodiscard]] constexpr CfemInt getComponentCountFromType(FieldType type)
    {
        switch (type)
        {
            case FieldType::Scalar:          return 1;
            case FieldType::Vector2D:        return 2;
            case FieldType::Vector3D:        return 3;
            case FieldType::SymTensor2D:     return 3;
            case FieldType::SkewSymTensor3D: return 3;
            case FieldType::Tensor2D:        return 4;
            case FieldType::AxiSymTensor3D: return 4;
            case FieldType::SymTensor3D:     return 6;
            case FieldType::Tensor3D:        return 9;
            default: 
                CFEM_FATAL << "Unsupported field type.";
                return 0;
        }
    }

class Mesh;

/**
 * @brief Base class for managing finite element basis functions 
 * and Degrees of Freedom (DOFs) numbering.
 */
class FESpace : public MeshObserver, public std::enable_shared_from_this<FESpace> {
public:
    enum class Continuity { CG, DG };

    virtual ~FESpace() = default;

    // Basic Accessors
    constexpr CfemInt   getNumComponents()           const noexcept { return m_numComponents; }
    constexpr FieldType getFieldType()               const noexcept { return m_fieldType; }
    constexpr Order     getOrder()                   const noexcept { return m_order; }
    constexpr CfemInt   getOrderInt()                const noexcept { return static_cast<CfemInt>(m_order); }
    constexpr const std::shared_ptr<Mesh>& getMesh() const noexcept { return m_mesh; }

    constexpr Continuity getContinuity() const noexcept { return m_continuity; }

    // Geometry & Reference Elements
    const ReferenceElement& getReferenceElement(CfemInt cellId) const;

    // DOF (Degrees of Freedom) Management
    
    /** 
     * @brief Total number of DOFs in the global system (size of PETSc vectors). 
     */
    virtual CfemInt getGlobalNumDofs() const = 0;

    /** 
     * @brief Number of DOFs associated with a specific cell. 
     */
    virtual CfemInt getDofsPerCell(CfemInt cellId) const = 0;

    /** 
     * @brief Retrieves global DOF indices for a given cell. 
     */
    virtual void getElementDofs(CfemInt cellId, std::vector<CfemInt>& globalIndices) const = 0;

    /** 
     * @brief Assigns global indices to all DOFs in the space. 
     */
    virtual void enumerateDofs() = 0;

    // Abstraction of DOF Mapping (Agnostic to CG/DG and memory layout)

    /** 
     * @brief Number of unique locations where DOFs are evaluated (e.g., active nodes in CG). 
     */
    virtual CfemInt getNumUniquePositions() const = 0;

    /** 
     * @brief Maps a cell's local node to a unique position ID (useful for locks/visited arrays). 
     */
    virtual CfemInt getPositionId(CfemInt cellId, CfemInt localNodeIdx) const = 0;

    /** 
     * @brief Provides the exact global PETSc index for a given cell, local node, and component. 
     */
    virtual CfemInt getGlobalDofIdx(CfemInt cellId, CfemInt localNodeIdx, CfemInt component) const = 0;

    /** 
     * @brief Retrieves a global DOF index from a physical Mesh Vertex ID (Optional/Legacy for CG). 
     */
    virtual CfemInt getGlobalDofFromMeshNode(CfemInt meshNodeId, CfemInt component) const = 0;

    // Updates and Callbacks

    /**
     * @brief Called by the Mesh during a topology update.
     * Triggers the update cascade down to the registered fields.
     */
    void updateAfterMeshChange() override;

    /** 
     * @brief Evaluates basis functions at a local coordinate xi. 
     */
    virtual void evaluateBasis(CfemInt cellId, const Vector3D& xi, std::vector<CfemReal>& values) const = 0;
    
    /** 
     * @brief Evaluates gradients of basis functions (physical or local as needed). 
     */
    virtual void evaluateBasisGradients(CfemInt cellId, const Vector3D& xi, std::vector<Vector3D>& grads) const = 0;

    // Support for Assembler and PETSc

    /** 
     * @brief Computes the number of non-zeros (nnz) per row for PETSc matrix pre-allocation. 
     */
    virtual void computeSparsityPattern(std::vector<CfemInt>& nnz) const = 0;

    /** 
     * @brief Identifies DOFs located on a specific boundary (using Gmsh tags). 
     */
    virtual void getEntityDofs(CfemInt internalEntityId, std::vector<CfemInt>& boundaryDofs) const = 0;

    /** 
     * @brief Registers a field to be notified when the mesh or space changes.
     */
    void registerField(std::shared_ptr<FESpaceObserver> field);
    // void registerField(std::shared_ptr<FEField> field);

protected:
    FESpace(std::shared_ptr<Mesh> mesh, Order order, Continuity cont, FieldType fieldType)
        : m_mesh(std::move(mesh)), m_order(order), m_continuity(cont), m_fieldType(fieldType), m_numComponents(getComponentCountFromType(fieldType))
    { }

    std::shared_ptr<Mesh> m_mesh;
    Order m_order;
    Continuity m_continuity;
    FieldType m_fieldType;
    CfemInt m_numComponents; 
    CfemInt m_spatialDimension;


    /// Fields associated with this space
    std::vector<std::weak_ptr<FESpaceObserver>> m_observers;

    std::map<CellType, std::unique_ptr<ReferenceElement>> m_refElements;

    virtual void validate() const;
};

} // namespace cfem

#endif