#pragma once 

#include "cfemutils.h"
#include "Vector3D.h"
#include "TensorO2.h"
#include "CfemExpression.h"
#include "EvalContext.h"
#include <vector>
#include <span>
#include <functional>
#include <memory>

namespace cfem::sym
{

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class LambdaExpression : public CfemExpression
    {
        std::function<void(const Vector3D&, std::span<CfemReal>)> m_func;       //!
        CfemInt                                                   m_dim;        //!
        std::shared_ptr<const CfemExpression>                     m_gradExpr;   //!

    public:
        LambdaExpression(
                    std::function<void(const Vector3D&, std::span<CfemReal>)> f,
                    CfemInt                                                   dim,
                    std::shared_ptr<const CfemExpression>                     gradExpr = nullptr);

        CfemInt getDimension() const override { return m_dim; }
        bool isPointwise() const override { return true; }

        void evaluate(CfemInt, const Vector3D&, std::span<CfemReal> out, EvalContext& ctx) const override {
            CFEM_ASSERT(out.size() >= m_dim, "LambdaExpression: output span too small.");
            m_func(ctx.physPos, out);
        }

        std::shared_ptr<const CfemExpression> grad() const override { return m_gradExpr; }
        std::string toString() const override { return "LambdaExpression_dim_"+std::to_string(m_dim) ; }
        ExpressionType getType() const override { return ExpressionType::Lambda; }
        
        static std::shared_ptr<const LambdaExpression> fromScalar(
                                                    const std::function<CfemReal(const Vector3D&)>& f,
                                                    std::shared_ptr<const CfemExpression>           gradient = nullptr);

        static std::shared_ptr<const LambdaExpression> fromVector(
                                                    const std::function<Vector3D(const Vector3D&)>& f,
                                                    std::shared_ptr<const CfemExpression>           gradient = nullptr);

        static std::shared_ptr<const LambdaExpression> fromTensor(
                                                    const std::function<TensorO2(const Vector3D&)>& f,
                                                    std::shared_ptr<const CfemExpression>           gradient = nullptr);
                        
        
    };


}