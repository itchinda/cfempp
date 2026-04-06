#ifndef CFEM_LAGRANGESPACE_H
#define CFEM_LAGRANGESPACE_H

#include "FESpace.h"
#include "ReferenceElement.h"
#include "DynamicVector.h"
#include <vector>
#include <memory>
#include <map>

namespace cfem {

/**
 * @class LagrangeSpace
 * @brief Lagrange Finite Element Space with BLOCK-wise Degree of Freedom (DOF) numbering.
 * This class implements continuous Lagrange basis functions (P1, P2). 
 * The DOF numbering is organized by components (Block-sorting):
 * [All X-components, followed by all Y-components, followed by all Z-components].
 *
 */
class LagrangeSpace : public FESpace {
private:
    struct ConstructorToken {};

    /** @brief Assigns global indices to DOFs using block-wise ordering. */
    void enumerateDofs() override;
    
    /** @brief Maps mesh cell type to the number of nodes expected by the basis order. */
    int getNodesPerCell(CellType type) const;

    CfemInt m_totalDofs = 0;
    CfemInt m_numActiveNodes = 0; ///< Number of unique interpolation nodes in the mesh.
    
    /** @brief Mapping: Mesh Vertex ID -> Local node index (0 to m_numActiveNodes-1). */
    std::vector<CfemInt> m_nodeToLocalIdx; 
    
    /** @brief Collection of reference elements for basis function evaluation. */
    std::map<CellType, std::unique_ptr<ReferenceElement>> m_refElements;
    
    /** @brief Initializes reference elements based on the space order. */
    void initReferenceElements();

    std::vector<CfemInt> m_localToMeshNodeIdx;

    /** @brief Checks mesh compatibility (e.g., P2 space requires a Quadratic Mesh). */
    void validate() const override{
        FESpace::validate();
        MeshOrder meshOrder = m_mesh->getOrder();
        if (m_order == Order::P2 && (meshOrder==MeshOrder::None || meshOrder==MeshOrder::Linear) ) {

            CFEM_ERROR << "LagrangeSpace Order Mismatch: P2 Space requires a"
                       " Quadratic Mesh (order 2), but the provided mesh " <<\
                       (meshOrder == MeshOrder::Linear? "is linear." : "has not been assigned any order." );
        }
    }

public:
    /**
     * @brief Protected constructor used by the Factory.
     */
    LagrangeSpace(ConstructorToken, std::shared_ptr<Mesh> mesh, Order order, int dimension = 1);

    /**
     * @brief Factory method: The recommended way to instantiate a LagrangeSpace.
     * @param mesh Shared pointer to the geometry.
     * @param order Polynomial order (P1 or P2). P0 is not supported here.
     * @param field_type Scalar or Vectorial nature of the physics.
     * @return A shared pointer to the initialized and enumerated LagrangeSpace.
     * @throw CFEM_ERROR if order is P0 or if mesh cannot support requested order.
     */
    static std::shared_ptr<LagrangeSpace> create(std::shared_ptr<Mesh> mesh, 
                                                 Order order, 
                                                 FieldType field_type = FieldType::Scalar) {
        int dim = 1;
        if (field_type == FieldType::Vector3D) dim = 3;
        else if (field_type == FieldType::Vector2D) dim = 2;

        if (order == Order::P0) {
            CFEM_ERROR << "LagrangeSpace Error: P0 is a discontinuous space and is not supported by LagrangeSpace. "
                       "Please use DGSpace or upgrade the order to P1.";
        }

        if (order == Order::P2 && mesh->getOrder() == MeshOrder::Linear) {
            CFEM_INFO << "LagrangeSpace: Auto-upgrading mesh to quadratic to support P2 space.";
            mesh->upgradeToQuadratic(); 
        }

        // Instantiation & Initialization
        auto space = std::make_shared<LagrangeSpace>(ConstructorToken{}, mesh, order, dim);
        space->initReferenceElements();
        
        // Connect mesh updates to this space
        mesh->addObserver(space);
        
        space->validate();
        space->enumerateDofs();

        return space;
    }

    /**
     * @brief Create a Taylor-Hood finite element pair on a given mesh.
     * 
     * Velocity space: higher-order Lagrange (P2, vector field)
     * Pressure space: lower-order Lagrange (P1, scalar field)
     *
     * @param mesh The mesh on which to define the spaces.
     * @return std::pair<std::shared_ptr<LagrangeSpace>, std::shared_ptr<LagrangeSpace>>
     *         First: velocity space, Second: pressure space
     */
    static std::pair<std::shared_ptr<LagrangeSpace>, std::shared_ptr<LagrangeSpace>>
        TaylorHood(const std::shared_ptr<Mesh>& mesh) {
            auto velocitySpace = LagrangeSpace::create(mesh, Order::P2, FieldType::Vector3D);
            auto pressureSpace = LagrangeSpace::create(mesh, Order::P1, FieldType::Scalar);
            return {velocitySpace, pressureSpace}; 
        }


    ~LagrangeSpace() override = default;

    /** @brief Total size for PETSc vectors/matrices (m_numActiveNodes * dimension). */
    CfemInt getGlobalNumDofs() const override { return m_totalDofs; }
    
    int getDofsPerCell(CfemInt cellId) const override;
    
    void getElementDofs(CfemInt cellId, std::vector<CfemInt>& globalIndices) const override;
    
    /** @brief Converts a mesh vertex ID to the local nodal index (ignoring components). */
    CfemInt getLocalNodeIdx(CfemInt meshNodeId) const { return m_nodeToLocalIdx[meshNodeId]; }
    CfemInt getMeshNodeIdFromLocal(CfemInt localIdx) const { 
        return m_localToMeshNodeIdx[localIdx]; 
    }

    const ReferenceElement& getRefElement(CellType type) const;

    /** @brief Computes the number of non-zeros (nnz) per row for efficient PETSc matrix allocation. */
    void computeSparsityPattern(std::vector<CfemInt>& nnz) const override;
    
    /** @brief Retrieves global DOFs associated with a boundary entity. */
    void getBoundaryDofs(int internalEntityId, std::vector<CfemInt>& boundaryDofs) const override;
    
    /** @brief Returns the number of nodes before repeating for the next component. */
    CfemInt getNumActiveNodes() const { return m_numActiveNodes; }

    // Basis Function Evaluation
    
    /** @brief Evaluates Lagrange basis functions at local coordinate xi. */
    void evaluateBasis(CfemInt cellId, const Vector3D& xi, std::vector<CfemReal>& values) const override;
    
    /** @brief Evaluates gradients of basis functions at local coordinate xi. */
    void evaluateGradients(CfemInt cellId, const Vector3D& xi, std::vector<Vector3D>& grads) const override;
};

} // namespace cfem

#endif