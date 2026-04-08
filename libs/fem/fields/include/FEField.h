#ifndef CFEM_FEFIELD_H_2
#define CFEM_FEFIELD_H_2

#include "CfemExpression.h"
#include "GradExpression.h"
#include "EvalContext.h"
#include "FESpace.h"         // Dependency inverted from LagrangeSpace to FESpace
#include "PetscVector.h"
#include "Quadrature.h"
#include <string>
#include <vector>
#include <memory>
#include <set>
#include <span>

namespace cfem
{
    class FEField : public sym::CfemExpression
    {
    private:
        struct ConstructorToken {};
        std::string m_name;
        CfemInt m_dimension;
        std::shared_ptr<FESpace> m_space; // <-- MAGIE: L'interface prend le relais
        PetscVector m_data;
        mutable std::vector<CfemReal> m_localValues;

        mutable int m_boundIndex = -1;

        struct CacheEntry {
            const QuadratureRule *quad = nullptr;
            DynamicVector<ShapeInfo> v_shapes; 
        };

    public:
        // Constructor now accepts FESpace
        FEField(ConstructorToken, const std::string &name, std::shared_ptr<FESpace> space);
        virtual ~FEField() = default;

        static std::shared_ptr<FEField> create(const std::string &name, std::shared_ptr<FESpace> space) {
            auto field = std::make_shared<FEField>(ConstructorToken{}, name, space);
            space->registerField(field);
            return field;
        }

        const std::string &getName() const { return m_name; }
        const std::shared_ptr<FESpace> &getSpace() const { return m_space; }
        PetscVector &getData() { return m_data; }
        const PetscVector &getData() const { return m_data; }
        CfemInt getDimension() const override { return m_dimension; }
        bool dependsOnFEField() const override { return true; }
        bool dependsOnShape() const override { return true; }

        ShapeEvaluationFlags requiredShapeFlags() const override {
            return ShapeEvaluationFlags::Value;
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
            spaces.insert(m_space);
        }

        /**
         * @brief Lie le champ à son index de cache O(1).
         */
        void bind(const std::map<const FESpace*, int>& spaceToId) const override {
            auto it = spaceToId.find(this->getSpace().get());
            if (it != spaceToId.end()) {
                m_boundIndex = it->second;
            } else {
                m_boundIndex = -1;
            }
        }
        
        void prepare() const override { updateLocalBuffer(); }
        void evaluate(CfemInt cellId, const Vector3D& xi, std::span<CfemReal> out, EvalContext& ctx ) const override;
        
        // --- Data Management ---
        void reallocate();
        void updateLocalBuffer() const;
        CfemReal getAtNode(CfemInt meshNodeId, int component = 0) const;
        CfemReal getAt(CfemInt globalDofIdx, int component) const;
        void evaluateWithoutShape(CfemInt cellId, const Vector3D& xi, std::span<CfemReal> out, EvalContext &ctx) const;
        void evaluateWithShape(CfemInt cellId, const ShapeInfo& shape, std::span<CfemReal> out) const;
        void evaluateGradient(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const;
        void evaluateGradientWithoutCachedMapping(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const;
        void evaluateGradientWithCachedMapping(CfemInt cellId, const MappingInfo& mapping, std::span<CfemReal> out) const;
        
        std::shared_ptr<const CfemExpression> timeDerivative() const override {
            CFEM_ERROR << "Time derivative of a FEField is not yet implemented.";
            return nullptr;
        }

        std::shared_ptr<const sym::CfemExpression> grad() const override {
            return std::make_shared<sym::GradExpression>(shared_from_this());
        }

        void project(const sym::CfemExpression& recipe);
        void project(std::shared_ptr<const sym::CfemExpression> recipe) {
            if (recipe) project(*recipe);
            else CFEM_ERROR << "project requires a valid CfemExpression object.";
        }

        template <typename Func>
        requires (!std::is_base_of_v<sym::CfemExpression, std::decay_t<Func>> && 
                  !std::is_convertible_v<std::decay_t<Func>, std::shared_ptr<const sym::CfemExpression>>)
        void project(const Func& func);

        sym::ExpressionType getType() const override { return sym::ExpressionType::FEField; }
        std::string toString() const override { return "FEField::" + m_name; }
    };
} 

#include "FEField.tpp"
#endif