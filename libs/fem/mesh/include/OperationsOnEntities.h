#ifndef CFEM_OPERATIONSONENTITIES_
#define CFEM_OPERATIONSONENTITIES_

#include <vector>

#include "cfemutils.h"
#include "MeshEntity.h"

class Mesh;


class OperationsOnEntities {
public:

    explicit OperationsOnEntities(const Mesh* mesh,
                                  const std::vector<const MeshEntity*>& entities);

    explicit OperationsOnEntities(const Mesh* mesh,
                                  const MeshEntity* a,
                                  const MeshEntity* b);

    explicit OperationsOnEntities(const Mesh* mesh,
                                  const MeshEntity* entity);


    void computeUnion(MeshEntity& target);
    void computeIntersection(MeshEntity& target);
    void computeDifference(MeshEntity& target);
    void computeComplement(MeshEntity& target);

    CfemInt resultingDimension() const;


private:
    std::vector<const MeshEntity*> m_sources;
    const Mesh* m_mesh {nullptr};
};

#endif
