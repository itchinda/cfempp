#ifndef CFEM_LAGRANGESPACE_H_2
#define CFEM_LAGRANGESPACE_H_2

#include "FESpace.h"
#include "DynamicVector.h"
#include <map>

namespace cfem {

/**
 * @class LagrangeSpace
 * @brief Continuous (CG) Finite Element Space with BLOCK-wise DOF numbering.
 */
class LagrangeSpace : public FESpace {
private:
    struct ConstructorToken {};

    void enumerateDofs() override;
    int getNodesPerCell(CellType type) const;
    void initReferenceElements();

    CfemInt m_totalDofs = 0;
    CfemInt m_numActiveNodes = 0; 
    
    std::vector<CfemInt> m_nodeToLocalIdx; 
    std::vector<CfemInt> m_localToMeshNodeIdx;
    
    void validate() const override {
        FESpace::validate();
        MeshOrder meshOrder = m_mesh->getOrder();
        if (m_order == Order::P2 && (meshOrder == MeshOrder::None || meshOrder == MeshOrder::Linear)) {
            CFEM_ERROR << "LagrangeSpace Order Mismatch: P2 Space requires a Quadratic Mesh.";
        }
    }

public:
    LagrangeSpace(ConstructorToken, std::shared_ptr<Mesh> mesh, Order order, int dimension = 1);
    ~LagrangeSpace() override = default;

    static std::shared_ptr<LagrangeSpace> create(std::shared_ptr<Mesh> mesh, 
                                                 Order order, 
                                                 FieldType field_type = FieldType::Scalar);

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

    // --- Core FESpace Overrides ---
    CfemInt getGlobalNumDofs() const override { return m_totalDofs; }
    int getDofsPerCell(CfemInt cellId) const override;
    void getElementDofs(CfemInt cellId, std::vector<CfemInt>& globalIndices) const override;
    void computeSparsityPattern(std::vector<CfemInt>& nnz) const override;
    void getEntityDofs(int internalEntityId, std::vector<CfemInt>& boundaryDofs) const override;

    // --- NEW: Abstraction Overrides ---
    CfemInt getNumUniquePositions() const override { return m_numActiveNodes; }
    CfemInt getPositionId(CfemInt cellId, int localNodeIdx) const override;
    CfemInt getGlobalDofIdx(CfemInt cellId, int localNodeIdx, int component) const override;
    CfemInt getGlobalDofFromMeshNode(CfemInt meshNodeId, int component) const override;

    // --- Basis Function Evaluation ---
    void evaluateBasis(CfemInt cellId, const Vector3D& xi, std::vector<CfemReal>& values) const override;
    void evaluateGradients(CfemInt cellId, const Vector3D& xi, std::vector<Vector3D>& grads) const override;

    // --- Legacy Accessors (Optional, strictly internal to Lagrange) ---
    CfemInt getLocalNodeIdx(CfemInt meshNodeId) const { return m_nodeToLocalIdx[meshNodeId]; }
    CfemInt getMeshNodeIdFromLocal(CfemInt localIdx) const { return m_localToMeshNodeIdx[localIdx]; }
};

} // namespace cfem

#endif