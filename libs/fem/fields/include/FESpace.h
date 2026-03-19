#ifndef CFEM_FESPACE_H
#define CFEM_FESPACE_H

#include "Mesh.h"
#include <vector>
#include <memory>

namespace cfem {


enum class FieldType : CfemInt { Scalar = 1, Vector2D = 2, Vector3D = 3 };
enum class Order { P0 = 0, P1 = 1, P2 = 2 };

class FEField;

/**
 * @brief Base class for managing finite element basis functions 
 * and Degrees of Freedom (DOFs) numbering.
 */
class FESpace : public std::enable_shared_from_this<FESpace> {
public:
    enum class Continuity { CG, DG };

    virtual ~FESpace() = default;

    // --- Basic Accessors ---
    int getDimension() const { return m_dimension; }
    Order getOrder() const { return m_order; }
    const std::shared_ptr<Mesh>& getMesh() const { return m_mesh; }

    // --- DOF (Degrees of Freedom) Management ---
    
    /** @brief Total number of DOFs in the global system (size of PETSc vectors). */
    virtual CfemInt getGlobalNumDofs() const = 0;

    /** @brief Number of DOFs associated with a specific cell. */
    virtual int getDofsPerCell(CfemInt cellId) const = 0;

    /** @brief Retrieves global DOF indices for a given cell. */
    virtual void getElementDofs(CfemInt cellId, std::vector<CfemInt>& globalIndices) const = 0;

    /** @brief Assigns global indices to all DOFs in the space. */
    virtual void enumerateDofs() = 0;

    /**
     * @brief Called by the Mesh during a topology update.
     * Triggers the update cascade down to the registered fields.
     */
    virtual void updateAfterMeshChange();

    /** @brief Evaluates basis functions at a local coordinate xi. */
    virtual void evaluateBasis(CfemInt cellId, const Vector3D& xi, std::vector<CfemReal>& values) const = 0;
    
    /** @brief Evaluates gradients of basis functions (physical or local as needed). */
    virtual void evaluateGradients(CfemInt cellId, const Vector3D& xi, std::vector<Vector3D>& grads) const = 0;

    // --- Support for Assembler and PETSc ---

    /** @brief Computes the number of non-zeros (nnz) per row for PETSc matrix pre-allocation. */
    virtual void computeSparsityPattern(std::vector<CfemInt>& nnz) const = 0;

    /** @brief Identifies DOFs located on a specific boundary (using Gmsh tags). */
    virtual void getBoundaryDofs(CfemInt tag, std::vector<CfemInt>& boundaryDofs) const = 0;

    /** 
     * @brief Registers a field to be notified when the mesh or space changes.
     */
    void registerField(std::shared_ptr<FEField> field);

protected:
    /**
     * @param mesh Shared pointer to the geometry.
     * @param order Polynomial order (P0, P1, P2).
     * @param cont Continuity (CG or DG).
     * @param dim Dimension of the physics (1 for scalar, 2/3 for vectors).
     */
    FESpace(std::shared_ptr<Mesh> mesh, Order order, Continuity cont, int dim)
        : m_mesh(std::move(mesh)), m_order(order), m_continuity(cont), m_dimension(dim) {}

    std::shared_ptr<Mesh> m_mesh;
    Order m_order;
    Continuity m_continuity;
    int m_dimension; // e.g., 1 for scalar, 3 for vector mechanics

    /// Fields associated with this space
    std::vector<std::weak_ptr<FEField>> m_fields;

    virtual void validate() const {
        if (m_mesh->getNumCells() == 0) {
            CFEM_ERROR << "FESpace Error: Mesh is empty.";
        }
        if (m_dimension < 1) {
            CFEM_ERROR << "FESpace Error: Field dimension must be >= 1.";
        }
    }
};

} // namespace cfem

#endif