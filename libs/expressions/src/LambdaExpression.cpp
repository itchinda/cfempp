#include "LambdaExpression.h"
#include "ZeroExpression.h"

namespace cfem::sym
{
    LambdaExpression::LambdaExpression( std::function<void(const Vector3D&, CfemReal, std::span<CfemReal>)> f,
                                        CfemInt                                                             dim,
                                        std::shared_ptr<const CfemExpression>                               gradExpr,
                                        std::shared_ptr<const CfemExpression>                               timeDeriv )
        : m_func(std::move(f)), m_dim(dim), m_gradExpr(std::move(gradExpr)), m_timeDeriv(std::move(timeDeriv)) 
    {
        if (m_gradExpr) {
            CFEM_ASSERT(m_gradExpr->getDimension() == m_dim * 3, 
                ("LambdaExpression: Dimension mismatch. Function dim is " + std::to_string(m_dim) + 
                ", so gradient/jacobian must be " + std::to_string(m_dim * 3) + ".").c_str());
        }

        if (m_timeDeriv) {
            CFEM_ASSERT(m_timeDeriv->getDimension() == m_dim, 
                ("LambdaExpression: Dimension mismatch. Function dim is " + std::to_string(m_dim) + 
                ", so time derivative must have dimension " + std::to_string(m_dim) + ".").c_str());
        }
    }

    LambdaExpression::LambdaExpression(
                    std::function<void(const Vector3D&, std::span<CfemReal>)> f,
                    CfemInt                                                   dim,
                    std::shared_ptr<const CfemExpression>                     gradExpr):
        LambdaExpression( [f](const Vector3D& X, CfemReal, std::span<CfemReal> out){ f(X, out); }, 
                          dim, 
                          gradExpr, 
                          std::make_shared<const ZeroExpression>(dim))      
    {
    }


    std::shared_ptr<const LambdaExpression> LambdaExpression::fromScalar(
        const std::function<CfemReal(const Vector3D&)>& f,
        std::shared_ptr<const CfemExpression> gradient) 
    {
        if (gradient) {
            CFEM_ASSERT(gradient->getDimension() == 3, 
                ("fromScalar: Gradient of a scalar must have dimension 3 (vector). Got " + 
                std::to_string(gradient->getDimension()) + ".").c_str());
        }

        return std::make_shared<LambdaExpression>(
            [f](const Vector3D& X, std::span<CfemReal> out) { out[0] = f(X); },
            1,
            gradient
        );
    }

    std::shared_ptr<const LambdaExpression> LambdaExpression::fromVector(
        const std::function<Vector3D(const Vector3D&)>& f,
        std::shared_ptr<const CfemExpression> gradient)
    {
        if (gradient) {
            CFEM_ASSERT(gradient->getDimension() == 9,( 
                "fromVector: Dimension mismatch. Gradient of a 3D vector must have dimension 9 (3x3 tensor). Got " + 
                std::to_string(gradient->getDimension()) + ".").c_str());
        }
        
        return std::make_shared<LambdaExpression>(
            [f](const Vector3D& X, std::span<CfemReal> out) {
                Vector3D v = f(X);
                out[0] = v.x; out[1] = v.y; out[2] = v.z;
            },
            3,
            gradient
        );
    }

    std::shared_ptr<const LambdaExpression> LambdaExpression::fromTensor(
        const std::function<TensorO2(const Vector3D&)>& f,
        std::shared_ptr<const CfemExpression> gradient)
    {
        if (gradient) {
            CFEM_ASSERT(gradient->getDimension() == 27,( 
                "fromTensor: Gradient of a 2nd order tensor must have dimension 27 (3x3x3). Got " + 
                std::to_string(gradient->getDimension()) + ".").c_str());
        }

        return std::make_shared<LambdaExpression>(
            [f](const Vector3D& X, std::span<CfemReal> out)
            {
                auto t = f(X);
                CFEM_ASSERT(out.size() >= 9, "LambdaExpression: tensor output span too small.");
                for (int i=0; i<3; ++i) {
                    out[3*i+0] = t(i,0); out[3*i+1] = t(i,1); out[3*i+2] = t(i,2);
                }
            },
            9,
            gradient
        );
    }

    std::shared_ptr<const LambdaExpression> LambdaExpression::fromScalar(
        const std::function<CfemReal(const Vector3D&, CfemReal)>& f,
        std::shared_ptr<const CfemExpression> gradient,
        std::shared_ptr<const CfemExpression> timeDeriv) 
    {
        if (gradient) {
            CFEM_ASSERT(gradient->getDimension() == 3, 
                ("fromScalar: Gradient of a scalar must have dimension 3 (vector). Got " + 
                std::to_string(gradient->getDimension()) + ".").c_str());
        }
        
        if (timeDeriv) {
            CFEM_ASSERT(timeDeriv->getDimension() == 1,( 
                "fromScalar: Dimension mismatch. Time derivative of a scalar must have dimension 1. Got " + 
                std::to_string(timeDeriv->getDimension()) + ".").c_str()); // FIXED TYPO
        }

        return std::make_shared<LambdaExpression>(
            [f](const Vector3D& X, CfemReal time, std::span<CfemReal> out) { out[0] = f(X, time); },
            1,
            gradient,
            timeDeriv
        );
    }

    std::shared_ptr<const LambdaExpression> LambdaExpression::fromVector(
        const std::function<Vector3D(const Vector3D&, CfemReal)>& f,
        std::shared_ptr<const CfemExpression> gradient,
        std::shared_ptr<const CfemExpression> timeDeriv)
    {
        if (gradient) {
            CFEM_ASSERT(gradient->getDimension() == 9,( 
                "fromVector: Dimension mismatch. Gradient of a 3D vector must have dimension 9 (3x3 tensor). Got " + 
                std::to_string(gradient->getDimension()) + ".").c_str());
        }

        if (timeDeriv) {
            CFEM_ASSERT(timeDeriv->getDimension() == 3,( 
                "fromVector: Dimension mismatch. Time derivative of a 3D vector must have dimension 3. Got " + 
                std::to_string(timeDeriv->getDimension()) + ".").c_str()); // FIXED TYPO
        }
        
        return std::make_shared<LambdaExpression>(
            [f](const Vector3D& X, CfemReal time, std::span<CfemReal> out) {
                Vector3D v = f(X, time);
                out[0] = v.x; out[1] = v.y; out[2] = v.z;
            },
            3,
            gradient,
            timeDeriv
        );
    }

    /**
     * @brief Wrap a 2nd-order tensor lambda as a CfemExpression.
     * @param f Lambda returning a 9-component array (3x3 tensor) from a 3D point.
     * @param gradient (Optional) Expression representing the 3rd-order gradient tensor.
     * Must have a dimension of 27.
     * @param timeDeriv (Optional) Expression representing the time derivative of the tensor.
     * @return A shared pointer to a CfemExpression of dimension 9.
     */
    std::shared_ptr<const LambdaExpression> LambdaExpression::fromTensor(
        const std::function<TensorO2(const Vector3D&, CfemReal)>& f,
        std::shared_ptr<const CfemExpression> gradient,
        std::shared_ptr<const CfemExpression> timeDeriv)
    {
        if (gradient) {
            CFEM_ASSERT(gradient->getDimension() == 27,( 
                "fromTensor: Gradient of a 2nd order tensor must have dimension 27 (3x3x3). Got " + 
                std::to_string(gradient->getDimension()) + ".").c_str());
        }

        if (timeDeriv) {
            CFEM_ASSERT(timeDeriv->getDimension() == 9,( 
                "fromTensor: Dimension mismatch. Time derivative of a second order tensor must have dimension 9. Got " + 
                std::to_string(timeDeriv->getDimension()) + ".").c_str()); // FIXED TYPO
        }

        return std::make_shared<LambdaExpression>(
            [f](const Vector3D& X, CfemReal time, std::span<CfemReal> out)
            {
                auto t = f(X, time);
                CFEM_ASSERT(out.size() >= 9, "LambdaExpression: tensor output span too small.");
                for (int i=0; i<3; ++i) {
                    out[3*i+0] = t(i,0); out[3*i+1] = t(i,1); out[3*i+2] = t(i,2);
                }
            },
            9,
            gradient,
            timeDeriv
        );
    }

}