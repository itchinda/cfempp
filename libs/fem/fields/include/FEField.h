#ifndef CFEM_FEFIELD_H
#define CFEM_FEFIELD_H

#include "CfemExpression.h"
#include "EvalContext.h"
#include "LagrangeSpace.h"
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
        std::shared_ptr<LagrangeSpace> m_space;
        PetscVector m_data;
        mutable std::vector<CfemReal> m_localValues;

        template <typename T>
        static Vector3D fetchValue(const T &obj, CfemInt cellId, const Vector3D &xi,
                                   const Vector3D &physPos, std::span<CfemReal> buffer);

        struct CacheEntry {
            const QuadratureRule *quad = nullptr;
            DynamicVector<ShapeInfo> v_shapes; 
        };

    public:
        FEField(ConstructorToken, const std::string &name, std::shared_ptr<LagrangeSpace> space);
        virtual ~FEField() = default;

        static std::shared_ptr<FEField> create(const std::string &name, std::shared_ptr<LagrangeSpace> space) {
            auto field = std::make_shared<FEField>(ConstructorToken{}, name, space);
            space->registerField(field);
            return field;
        }

        const std::string &getName() const { return m_name; }
        const std::shared_ptr<LagrangeSpace> &getSpace() const { return m_space; }
        PetscVector &getData() { return m_data; }
        const PetscVector &getData() const { return m_data; }
        CfemInt getDimension() const override { return m_space->getDimension(); }

        void prepare() const override { updateLocalBuffer(); }
        void evaluate(CfemInt cellId, const Vector3D& xi, std::span<CfemReal> out, EvalContext& ctx ) const override;
        
        // --- Data Management ---
        void reallocate();
        void updateLocalBuffer() const;
        CfemReal getAtNode(CfemInt meshNodeId, int component = 0) const;
        CfemReal getAt(CfemInt globalDofIdx, int component) const;
        void evaluateAt(CfemInt cellId, const Vector3D& xi, std::span<CfemReal> out) const;
        void evaluateWithShape(CfemInt cellId, const ShapeInfo& shape, std::span<CfemReal> out) const;
        void evaluateGradientAt(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const;
        void evaluateGradientWithCachedMapping(CfemInt cellId, const MappingInfo& mapping, std::span<CfemReal> out) const;
        std::shared_ptr<const CfemExpression> timeDerivative() const override{
            CFEM_ERROR << "You requested the time derivative of a FEField, but logic is not yet implemented.";
            return nullptr;
        }

        // --- Nodal projections ---
        void project(const sym::CfemExpression& recipe);
        template <typename Func>
        requires (!std::is_same_v<std::decay_t<Func>, sym::CfemExpression> &&
                !std::is_same_v<std::decay_t<Func>, std::shared_ptr<sym::CfemExpression>>)
        void project(const Func& func);

        sym::ExpressionType getType() const override { return sym::ExpressionType::FEField; }


    };
} 

#include "FEField.tpp"
#endif