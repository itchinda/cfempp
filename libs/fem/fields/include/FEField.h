#ifndef CFEM_FEFIELD_H
#define CFEM_FEFIELD_H

#include "LagrangeSpace.h"
#include "PetscVector.h"
#include "Quadrature.h"
#include <string>
#include <vector>
#include <memory> // Obligatoire pour shared_ptr
#include <set>

namespace cfem
{

    /**
     * @class FEField
     * @brief Gère les données d'un champ physique.
     * Hérite de enable_shared_from_this pour pouvoir s'enregistrer auprès de son espace.
     */
    class FEField : public std::enable_shared_from_this<FEField>
    {

    private:
        inline static std::set<std::string> s_globalFieldRegistry;
        struct ConstructorToken
        {
        };
        std::string m_name;
        std::shared_ptr<LagrangeSpace> m_space;
        PetscVector m_data;
        mutable std::vector<CfemReal> m_localValues;

        // Helper: Square of scalar or Dot product of vector
        static inline CfemReal magnitude_sqr(CfemReal val) { return val * val; }
        static inline CfemReal magnitude_sqr(const Vector3D &val) { return val.dot(val); }

        // Helper: Uniform access to Field (via local eval) or Lambda (via physical pos)
        template <typename T>
        static Vector3D fetchValue(const T &obj, CfemInt cellId, const Vector3D &xi,
                                   const Vector3D &physPos, std::vector<CfemReal> &buffer)
        {
            // CASE A: The object is an FEField (P1, P2, etc.)
            if constexpr (std::is_convertible_v<T, const FEField &>)
            {
                obj.evaluateAt(cellId, xi, buffer);
                if (buffer.size() == 1)
                    return Vector3D{buffer[0], 0.0, 0.0};
                if (buffer.size() == 2)
                    return Vector3D{buffer[0], buffer[1], 0.0};
                return Vector3D{buffer[0], buffer[1], buffer[2]};
            }
            // CASE B: The object is a Lambda or Functor
            else
            {
                auto result = obj(physPos);
                using ReturnType = decltype(result);

                if constexpr (std::is_arithmetic_v<ReturnType>)
                {
                    // It's a scalar (double/float) -> Promote to Vector3D
                    return Vector3D{static_cast<CfemReal>(result), 0.0, 0.0};
                }
                else
                {
                    // It's already a Vector3D (or compatible) -> Use as is
                    return result;
                }
            }
        }

        struct CacheEntry
            {
                const QuadratureRule *quad = nullptr;
                DynamicVector<ShapeInfo> v_shapes; // [qp_index][node_index]
            };

    public:
        FEField(ConstructorToken, const std::string &name, std::shared_ptr<LagrangeSpace> space);
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
        static std::shared_ptr<FEField> create(const std::string &name, std::shared_ptr<LagrangeSpace> space)
        {
            auto field = std::make_shared<FEField>(ConstructorToken{}, name, space);

            NamingPolicy::validate(name);

            if (s_globalFieldRegistry.find(name) != s_globalFieldRegistry.end())
            {
                CFEM_ERROR(" Global naming collision: Field with name '" + name + "' already exists.");
            }

            space->registerField(field);
            s_globalFieldRegistry.insert(name);

            return field;
        }

        // --- Accès aux données ---
        const std::string &getName() const { return m_name; }

        // On utilise ta nouvelle convention : retour par référence constante de shared_ptr
        const std::shared_ptr<LagrangeSpace> &getSpace() const { return m_space; }

        PetscVector &getData() { return m_data; }
        const PetscVector &getData() const { return m_data; }
        CfemInt getDimension() const { return m_space->getDimension(); }

        /**
         * @brief Réalloue le vecteur PETSc si la taille de l'espace change.
         * Cette méthode sera appelée par l'espace lors d'un upgrade du maillage.
         */
        void reallocate();

        void updateLocalBuffer() const;
        CfemReal getAtNode(CfemInt meshNodeId, int component = 0) const;
        CfemReal getAt(CfemInt globalDofIdx, int component) const;

        void evaluateAt(CfemInt cellId, const Vector3D &xi, std::vector<CfemReal> &outValues) const;

        /**
         * @brief Projects an analytical function onto the nodal degrees of freedom.
         * @param func A lambda or function: Vector3D (position) -> Vector3D or double
         */
        template <typename Func>
        void project(Func func)
        {
            auto space = getSpace();      // LagrangeSpace
            auto mesh = space->getMesh(); // The mesh Shared Pointer
            auto &vec = getData();        // The PetscVector

            CfemInt N = space->getNumActiveNodes();
            int dim = space->getDimension();
            CfemInt totalVertices = mesh->getNumVertices();

            for (CfemInt i = 0; i < totalVertices; ++i)
            {
                CfemInt nodeIdx = space->getLocalNodeIdx(i);
                if (nodeIdx == -1)
                    continue;

                const auto &vertex = mesh->getVertex(i);

                // Execute the user's lambda
                auto result = func(vertex.coords);

                // Map the result to the PETSc vector based on dimension
                if constexpr (std::is_same_v<decltype(result), Vector3D>)
                {
                    for (int d = 0; d < dim; ++d)
                    {
                        // Block-wise mapping: [nodeIdx + d * N]
                        vec.setValue(nodeIdx + d * N, result[d], INSERT_VALUES);
                    }
                }
                else
                {
                    // Scalar case
                    vec.setValue(nodeIdx, static_cast<double>(result), INSERT_VALUES);
                }
            }
            vec.assemble();
        }

        /**
         * @brief Computes the L2 Norm of the field: ||u|| = sqrt( integral( u * u ) )
         */
        CfemReal normL2(IntegrationSchemeParams scheme_params = {}) const
        {
            return L2NormDistance(*this, [](const Vector3D &)
                                  { return 0.0; }, scheme_params);
        }

        // /**
        //  * @brief Computes ||u - v||_{L2}.
        //  * Handles combinations of (FEField vs FEField) and (FEField vs Lambda).
        //  * @tparam T1 Can be FEField or a Lambda/Functor.
        //  * @tparam T2 Can be FEField or a Lambda/Functor.
        //  */
        // template <typename T1, typename T2>
        // static CfemReal L2NormDistance(const T1& u, const T2& v, IntegrationSchemeParams scheme_params = {}) {

        //     CfemReal totalErrorSqr = 0.0;
        //     PROFILE_SCOPE("FEField::L2NormDistance");
        //     {
        //         // 1. Resolve Context: Find which argument provides the Mesh
        //         const FEField* contextField = nullptr;
        //         if constexpr (std::is_convertible_v<T1, const FEField&>)
        //         {
        //             contextField = &u; u.updateLocalBuffer();
        //         }
        //         else if constexpr (std::is_convertible_v<T2, const FEField&>) {
        //             contextField = &v;
        //             v.updateLocalBuffer();
        //         }
        //         else CFEM_ERROR("L2Diff requires at least one FEField argument.");

        //         auto mesh = contextField->getSpace()->getMesh();

        //         IntegrationScheme scheme;
        //         ShapeInfo shape;
        //         MappingInfo map;
        //         DynamicVector<Vector3D> physCoords; // Buffer for element vertices
        //         std::vector<CfemReal> fieldBuffer;  // Buffer for field evaluation

        //         for (CfemInt i = 0; i < mesh->getNumCells(); ++i) {
        //             const auto& cell = mesh->getCell(i);

        //             // Use high-order quadrature (e.g., 5-point) to capture P2*P2 interactions
        //             const auto& quad = scheme.getRule(cell.type, scheme_params.getNumIntgPt(cell.type) );

        //             mesh->getCellVerticesCoords(i, physCoords);
        //             const auto& refEl = contextField->getSpace()->getRefElement(cell.type);

        //             CfemReal cellIntegral = 0.0;

        //             for (const auto& qp : quad.getPoints()) {
        //                 // Compute Geometry (Reference -> Physical)
        //                 refEl.computeMapping(qp.xi, physCoords, ShapeEvaluationFlags::Value, shape, map);

        //                 // Fetch Values using the generic helper'

        //                 auto valU = fetchValue(u, i, qp.xi, map.physicalPosition, fieldBuffer);
        //                 auto valV = fetchValue(v, i, qp.xi, map.physicalPosition, fieldBuffer);

        //                 auto diff = valU - valV;

        //                 // Accumulate: works for Scalar (diff*diff) or Vector (diff.dot(diff))
        //                 cellIntegral += magnitude_sqr(diff) * std::abs(map.detJ) * qp.weight;
        //             }
        //             totalErrorSqr += cellIntegral;
        //         }
        //     }
        //     return std::sqrt(totalErrorSqr);
        // }

        /**
         * @brief Optimized L2 Norm distance calculation ||u - v||.
         * Uses a flat cache for basis functions and OpenMP for parallelization.
         */
        template <typename T1, typename T2>
        static CfemReal L2NormDistance(const T1 &u, const T2 &v, IntegrationSchemeParams scheme_params = {})
        {
            CfemReal totalErrorSqr = 0.0;
            PROFILE_SCOPE("FEField::L2NormDistance");

            // 1. Resolve Context: Find which argument provides the Mesh
            const FEField *contextField = nullptr;
            if constexpr (std::is_convertible_v<T1, const FEField &>)
            {
                contextField = &u;
                u.updateLocalBuffer();
            }
            else if constexpr (std::is_convertible_v<T2, const FEField &>)
            {
                contextField = &v;
                v.updateLocalBuffer();
            }
            else
            {
                CFEM_ERROR("L2Diff requires at least one FEField argument.");
            }

            auto mesh = contextField->getSpace()->getMesh();
            auto space = contextField->getSpace();
            IntegrationScheme scheme;

            // 2. Build the Flat Cache (O(1) access)
            // Structure to hold pre-computed reference element data
        

            // Allocate the full enum range (13 slots based on your enum)
            DynamicVector<CacheEntry> fastCache(static_cast<int>(CellType::COUNT));

            // Pre-fill only the types present in the mesh
            ShapeEvaluationFlags internalFields = ShapeEvaluationFlags::Value | ShapeEvaluationFlags::Gradient;
            
            for (CellType type : mesh->getCellTypesPresent())
            {
                int idx = static_cast<int>(type);
                const auto &quad = scheme.getRule(type, scheme_params.getNumIntgPt(type));
                const auto &refEl = space->getRefElement(type);

                const DynamicVector<QuadraturePoint>& lQuadPts = quad.getPoints();

                fastCache[idx].quad = &quad;
                fastCache[idx].v_shapes.resize(lQuadPts.size());

                // Pre-compute N_j(xi_q) for the reference element
                for (size_t q = 0; q < lQuadPts.size(); ++q)
                {
                    refEl.evaluate(lQuadPts[q].xi, internalFields, fastCache[idx].v_shapes[q]);
                    fastCache[idx].v_shapes[q].shapeMustBeUpdated = false;
                }


            }

// 3. Parallel Computation Loop
// reduction(+:totalErrorSqr) ensures thread-safe summation
#pragma omp parallel reduction(+ : totalErrorSqr)
            {
                // Thread-local buffers to prevent costly reallocations in the inner loop
                MappingInfo map;
                ShapeInfo shape;
                DynamicVector<Vector3D> physCoords;
                std::vector<CfemReal> fieldBuffer;

#pragma omp for
                for (CfemInt i = 0; i < mesh->getNumCells(); ++i)
                {
                    const auto &cell = mesh->getCell(i);
                    auto &cache = fastCache[static_cast<int>(cell.type)];
                    const auto &refEl = space->getRefElement(cell.type);

                    // Get physical coordinates of the cell vertices
                    mesh->getCellVerticesCoords(i, physCoords);
                    CfemReal cellIntegral = 0.0;

                    // Integration Loop
                    for (size_t q = 0; q < cache.quad->getPoints().size(); ++q)
                    {
                        const auto &qp = cache.quad->getPoints()[q];

                        // Geometry Mapping: Reference -> Physical
                        // Note: We use ShapeEvaluationFlags::None because we use our cache for N(xi)
                        refEl.computeMapping(qp.xi, physCoords, ShapeEvaluationFlags::None, cache.v_shapes[q], map);
                        // Fetch Values: FEField uses local evaluations, Lambdas use physical positions
                        auto valU = fetchValue(u, i, qp.xi, map.physicalPosition, fieldBuffer);
                        auto valV = fetchValue(v, i, qp.xi, map.physicalPosition, fieldBuffer);


                        auto diff = valU - valV;

                        // Accumulate: magnitude_sqr * |detJ| * weight
                        cellIntegral += magnitude_sqr(diff) * std::abs(map.detJ) * qp.weight;
                    }
                    totalErrorSqr += cellIntegral;
                }
            }

            return std::sqrt(totalErrorSqr);
        }
    };

} // namespace cfem

#endif