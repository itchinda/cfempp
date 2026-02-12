#ifndef CFEM_FEFIELD_H
#define CFEM_FEFIELD_H

#include "LagrangeSpace.h"
#include "PetscVector.h" 
#include <string>
#include <vector>
#include <memory> // Obligatoire pour shared_ptr
#include <set>

namespace cfem {



/**
 * @class FEField
 * @brief Gère les données d'un champ physique.
 * Hérite de enable_shared_from_this pour pouvoir s'enregistrer auprès de son espace.
 */
class FEField : public std::enable_shared_from_this<FEField> {

private:
    inline static std::set<std::string> s_globalFieldRegistry;
    struct ConstructorToken {};
    std::string m_name;
    std::shared_ptr<LagrangeSpace> m_space;
    PetscVector m_data; 
    mutable std::vector<CfemReal> m_localValues; 

public:

    FEField(ConstructorToken, const std::string& name, std::shared_ptr<LagrangeSpace> space);
    ~FEField() = default;

    /**
     * @brief Factory method: Creates an FEField and automatically registers it with its space.
     * @details This is the only way to instantiate an FEField. By using this factory:
     * 1. Memory is safely managed via std::shared_ptr.
     * 2. The field is automatically registered in the FESpace's update list.
     * 3. A strong link is established between the Lagrange mathematics and the physical data.
     * @param name  A unique identifier for the field (e.g., "Pressure", "Velocity").
     * @param space A shared pointer to the LagrangeSpace defining the field's topology.
     * @return      A shared pointer to the initialized and registered FEField.
     */
    static std::shared_ptr<FEField> create(const std::string& name, std::shared_ptr<LagrangeSpace> space) {
        auto field = std::make_shared<FEField>(ConstructorToken{}, name, space);

        NamingPolicy::validate(name);

        if (s_globalFieldRegistry.find(name) != s_globalFieldRegistry.end()) {
            CFEM_ERROR(" Global naming collision: Field with name '" + name + "' already exists.");
        }


        space->registerField(field);
        s_globalFieldRegistry.insert(name);

        return field;
    }

    // --- Accès aux données ---
    const std::string& getName() const { return m_name; }
    
    // On utilise ta nouvelle convention : retour par référence constante de shared_ptr
    const std::shared_ptr<LagrangeSpace>& getSpace() const { return m_space; }

    PetscVector& getData() { return m_data; }
    const PetscVector& getData() const { return m_data; }
    CfemInt getDimension() const { return m_space->getDimension(); }

    /**
     * @brief Réalloue le vecteur PETSc si la taille de l'espace change.
     * Cette méthode sera appelée par l'espace lors d'un upgrade du maillage.
     */
    void reallocate();

    void updateLocalBuffer() const;
    CfemReal getAtNode(CfemInt meshNodeId, int component = 0) const;
    CfemReal getAt(CfemInt globalDofIdx, int component) const;

    void evaluateAt(CfemInt cellId, const Vector3D& xi, std::vector<CfemReal>& outValues) const;

    /**
     * @brief Projects an analytical function onto the nodal degrees of freedom.
     * @param func A lambda or function: Vector3D (position) -> Vector3D or double
     */
    template <typename Func>
    void project(Func func) {
        auto space = getSpace(); // LagrangeSpace
        auto mesh = space->getMesh();  // The mesh Shared Pointer
        auto& vec = getData(); // The PetscVector
        
        CfemInt N = space->getNumActiveNodes();
        int dim = space->getDimension();
        CfemInt totalVertices = mesh->getNumVertices();

        for (CfemInt i = 0; i < totalVertices; ++i) {
            CfemInt nodeIdx = space->getLocalNodeIdx(i);
            if (nodeIdx == -1) continue;

            const auto& vertex = mesh->getVertex(i);
            
            // Execute the user's lambda
            auto result = func(vertex.coords);

            // Map the result to the PETSc vector based on dimension
            if constexpr (std::is_same_v<decltype(result), Vector3D>) {
                for (int d = 0; d < dim; ++d) {
                    // Block-wise mapping: [nodeIdx + d * N]
                    vec.setValue(nodeIdx + d * N, result[d], INSERT_VALUES);
                }
            } else {
                // Scalar case
                vec.setValue(nodeIdx, static_cast<double>(result), INSERT_VALUES);
            }
        }
        vec.assemble();
    }

};

} // namespace cfem

#endif