#ifndef CFEM_FESPACE_H
#define CFEM_FESPACE_H

#include "cfemutils.h"
#include "cfem_topology.h"
#include "ReferenceElement.h"
#include <vector>
#include <memory>

namespace cfem {

enum class FieldType : CfemInt { Scalar = 1, Vector2D = 2, Vector3D = 3 };

class Mesh;

/**
 * @brief Base class for managing finite element basis functions 
 * and Degrees of Freedom (DOFs) numbering.
 */
class FESpace : public MeshObserver, public std::enable_shared_from_this<FESpace> {
public:
    enum class Continuity { CG, DG };

    virtual ~FESpace() = default;

    // --- Basic Accessors ---
    int getDimension() const { return m_dimension; }
    Order getOrder() const { return m_order; }
    CfemInt getOrderInt() const { return static_cast<CfemInt>(m_order); }
    const std::shared_ptr<Mesh>& getMesh() const { return m_mesh; }

    Continuity getContinuity() const { return m_continuity; }

    // --- Geometry & Reference Elements ---
    const ReferenceElement& getRefElement(CellType type) const;

    // --- DOF (Degrees of Freedom) Management ---
    
    /** @brief Total number of DOFs in the global system (size of PETSc vectors). */
    virtual CfemInt getGlobalNumDofs() const = 0;

    /** @brief Number of DOFs associated with a specific cell. */
    virtual int getDofsPerCell(CfemInt cellId) const = 0;

    /** @brief Retrieves global DOF indices for a given cell. */
    virtual void getElementDofs(CfemInt cellId, std::vector<CfemInt>& globalIndices) const = 0;

    /** @brief Assigns global indices to all DOFs in the space. */
    virtual void enumerateDofs() = 0;

    // --- NEW: Abstraction of DOF Mapping (Agnostic to CG/DG and memory layout) ---

    /** @brief Number of unique locations where DOFs are evaluated (e.g., active nodes in CG). */
    virtual CfemInt getNumUniquePositions() const = 0;

    /** @brief Maps a cell's local node to a unique position ID (useful for locks/visited arrays). */
    virtual CfemInt getPositionId(CfemInt cellId, int localNodeIdx) const = 0;

    /** @brief Provides the exact global PETSc index for a given cell, local node, and component. */
    virtual CfemInt getGlobalDofIdx(CfemInt cellId, int localNodeIdx, int component) const = 0;

    /** @brief Retrieves a global DOF index from a physical Mesh Vertex ID (Optional/Legacy for CG). */
    virtual CfemInt getGlobalDofFromMeshNode(CfemInt meshNodeId, int component) const = 0;

    // --- Updates and Callbacks ---

    /**
     * @brief Called by the Mesh during a topology update.
     * Triggers the update cascade down to the registered fields.
     */
    void updateAfterMeshChange() override;

    /** @brief Evaluates basis functions at a local coordinate xi. */
    virtual void evaluateBasis(CfemInt cellId, const Vector3D& xi, std::vector<CfemReal>& values) const = 0;
    
    /** @brief Evaluates gradients of basis functions (physical or local as needed). */
    virtual void evaluateGradients(CfemInt cellId, const Vector3D& xi, std::vector<Vector3D>& grads) const = 0;

    // --- Support for Assembler and PETSc ---

    /** @brief Computes the number of non-zeros (nnz) per row for PETSc matrix pre-allocation. */
    virtual void computeSparsityPattern(std::vector<CfemInt>& nnz) const = 0;

    /** @brief Identifies DOFs located on a specific boundary (using Gmsh tags). */
    virtual void getEntityDofs(int internalEntityId, std::vector<CfemInt>& boundaryDofs) const = 0;

    /** * @brief Registers a field to be notified when the mesh or space changes.
     */
    void registerField(std::shared_ptr<FESpaceObserver> field);
    // void registerField(std::shared_ptr<FEField> field);

protected:
    FESpace(std::shared_ptr<Mesh> mesh, Order order, Continuity cont, int dim)
        : m_mesh(std::move(mesh)), m_order(order), m_continuity(cont), m_dimension(dim) 
    {
        initReferenceElements();
    }

    std::shared_ptr<Mesh> m_mesh;
    Order m_order;
    Continuity m_continuity;
    int m_dimension; 

    /// Fields associated with this space
    std::vector<std::weak_ptr<FESpaceObserver>> m_observers;

    std::map<CellType, std::unique_ptr<ReferenceElement>> m_refElements;
    void initReferenceElements();

    virtual void validate() const;
};

} // namespace cfem

#endif