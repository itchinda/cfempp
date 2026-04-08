// #pragma once

// #include "CfemExpression.h"
// #include "FESpace.h" // Ensure this is included or forward declared
// #include "EvalContext.h"
// #include <memory>
// #include <string>

// namespace cfem::sym {

// enum class ArgumentType { Trial, Test };

// class Argument : public CfemExpression {
// protected:
//     std::shared_ptr<FESpace> m_space;
//     ArgumentType m_argType;
//     std::string m_name;

// public:
//     Argument(std::shared_ptr<LagrangeSpace> space, ArgumentType type, const std::string& name)
//         : m_space(std::move(space)), m_argType(type), m_name(name) {}

//     virtual ~Argument() = default;

//     const std::string &getName() const { return m_name; }
//     const std::shared_ptr<FESpace>& getSpace() const { return m_space; }
//     ArgumentType getArgumentType() const { return m_argType; }

//     CfemInt getDimension() const override { return m_space->getDimension(); }
    
//     virtual ExpressionType getType() const override {
//         return (m_argType == ArgumentType::Trial) ? ExpressionType::Trial : ExpressionType::Test;
//     }

//     void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
//         spaces.insert(m_space);
//     }

//     std::string toString() const override { return m_name; }

//     std::shared_ptr<const CfemExpression> timeDerivative() const override {
//         CFEM_ERROR << "Time derivative not implemented for Argument (" << m_name << ").";
//         return nullptr;
//     }

//     // --- Standard Evaluation (Returns 0 for pure symbolic arguments) ---
//     void evaluate(CfemInt cellId, const Vector3D& xi, std::span<CfemReal> out, EvalContext& ctx) const override {
//         std::fill(out.begin(), out.end(), 0.0);
//     }

//     // --- The Explicit Basis Interface (Used by the Assembler) ---
//     virtual void evaluateBasis(CfemInt dof, std::span<CfemReal> out, EvalContext& ctx) const {
//         const ShapeInfo* shape = (m_argType == ArgumentType::Trial) ? ctx.trialShape : ctx.testShape;
//         if (!shape || dof < 0) {
//             CFEM_ERROR << "Argument (" << m_name << "): missing shape context.";
//         }

//         const CfemInt dim = m_space->getDimension();
//         const CfemInt numNodes = m_space->getNumActiveNodes();
//         const CfemInt comp = dof / numNodes;
//         const CfemInt node = dof % numNodes;

//         std::fill(out.begin(), out.end(), 0.0);
//         out[comp] = shape->values[node];
//     }

//     virtual void evaluateBasisGradient(CfemInt dof, std::span<CfemReal> out, EvalContext& ctx) const {
//         const MappingInfo* mapping = (m_argType == ArgumentType::Trial) ? ctx.trialMapping : ctx.testMapping;
//         if (!mapping || dof < 0) {
//             CFEM_ERROR << "Argument (" << m_name << "): missing mapping context.";
//         }

//         const CfemInt numNodes = m_space->getNumActiveNodes();
//         const CfemInt comp = dof / numNodes;
//         const CfemInt node = dof % numNodes;

//         const auto& g = mapping->physicalShapeGradients[node];

//         std::fill(out.begin(), out.end(), 0.0);
//         out[comp * 3 + 0] = g.x;
//         out[comp * 3 + 1] = g.y;
//         out[comp * 3 + 2] = g.z;
//     }
// };

// class TrialFunction : public Argument {
// public:
//     TrialFunction(std::shared_ptr<LagrangeSpace> space, const std::string& name = "u")
//         : Argument(std::move(space), ArgumentType::Trial, name) {}
//     bool dependsOnFEField() const override { return false; }
//     bool dependsOnShape() const override { return true; }
//     ShapeEvaluationFlags requiredShapeFlags() const override {
//         return ShapeEvaluationFlags::Value;
//     }

//     virtual ExpressionType getType() const override { return ExpressionType::Trial; }
// };

// class TestFunction : public Argument {
// public:
//     TestFunction(std::shared_ptr<LagrangeSpace> space, std::string name = "v")
//         : Argument(std::move(space), ArgumentType::Test, std::move(name)) {}

//     bool dependsOnFEField() const override { return false; }
//     bool dependsOnShape() const override { return true; }
//     ShapeEvaluationFlags requiredShapeFlags() const override {
//         return ShapeEvaluationFlags::Value;
//     }

//     virtual ExpressionType getType() const override { return ExpressionType::Test; }
// };

// } // namespace cfem::sym


#pragma once

#include "CfemExpression.h"
#include "FESpace.h" 
#include "EvalContext.h"
#include <memory>
#include <string>

namespace cfem::sym {

enum class ArgumentType { Trial, Test };

/**
 * @class Argument
 * @brief Represents a Trial (u) or Test (v) function in a variational form.
 * * Updated to be agnostic of the FESpace type (CG or DG).
 */
class Argument : public CfemExpression {
protected:
    std::shared_ptr<FESpace> m_space;
    ArgumentType m_argType;
    std::string m_name;

public:
    // Update: Accept base FESpace instead of LagrangeSpace
    Argument(std::shared_ptr<FESpace> space, ArgumentType type, const std::string& name)
        : m_space(std::move(space)), m_argType(type), m_name(name) {}

    virtual ~Argument() = default;

    const std::string &getName() const { return m_name; }
    const std::shared_ptr<FESpace>& getSpace() const { return m_space; }
    ArgumentType getArgumentType() const { return m_argType; }

    CfemInt getDimension() const override { return m_space->getDimension(); }
    
    virtual ExpressionType getType() const override {
        return (m_argType == ArgumentType::Trial) ? ExpressionType::Trial : ExpressionType::Test;
    }

    void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
        spaces.insert(m_space);
    }

    std::string toString() const override { return m_name; }

    std::shared_ptr<const CfemExpression> timeDerivative() const override {
        CFEM_ERROR << "Time derivative not implemented for Argument (" << m_name << ").";
        return nullptr;
    }

    // Standard evaluation returns 0 because an Argument is a basis function placeholder, 
    // not a field with values.
    void evaluate(CfemInt cellId, const Vector3D& xi, std::span<CfemReal> out, EvalContext& ctx) const override {
        std::fill(out.begin(), out.end(), 0.0);
    }

    /**
     * @brief Evaluates the basis function value for a specific local DOF.
     * @param localDofIdx The index of the DOF within the current cell [0 ... DofsPerCell-1]
     */
    virtual void evaluateBasis(CfemInt localDofIdx, std::span<CfemReal> out, EvalContext& ctx) const {
        const ShapeInfo* shape = (m_argType == ArgumentType::Trial) ? ctx.trialShape : ctx.testShape;
        
        if (!shape || localDofIdx < 0) {
            CFEM_ERROR << "Argument (" << m_name << "): invalid context or DOF index.";
        }

        const int dim = m_space->getDimension();
        // Crucial logic change:
        // We assume the assembler provides localDofIdx in [0, nodes * dim - 1]
        // mapped as: Node0_Comp0, Node0_Comp1... NodeN_CompM
        const int numNodes = static_cast<int>(shape->values.size());
        const int comp = static_cast<int>(localDofIdx % dim);
        const int node = static_cast<int>(localDofIdx / dim);

        std::fill(out.begin(), out.end(), 0.0);
        out[comp] = shape->values[node];
    }

    /**
     * @brief Evaluates the basis function gradient for a specific local DOF.
     */
    virtual void evaluateBasisGradient(CfemInt localDofIdx, std::span<CfemReal> out, EvalContext& ctx) const {
        const MappingInfo* mapping = (m_argType == ArgumentType::Trial) ? ctx.trialMapping : ctx.testMapping;
        
        if (!mapping || localDofIdx < 0) {
            CFEM_ERROR << "Argument (" << m_name << "): invalid mapping context.";
        }

        const int dim = m_space->getDimension();
        const int node = static_cast<int>(localDofIdx / dim);
        const int comp = static_cast<int>(localDofIdx % dim);

        const auto& g = mapping->physicalShapeGradients[node];

        std::fill(out.begin(), out.end(), 0.0);
        // Map to vector/tensor output span: comp*3 + spatial_dir
        out[comp * 3 + 0] = g.x;
        out[comp * 3 + 1] = g.y;
        out[comp * 3 + 2] = g.z;
    }
};

class TrialFunction : public Argument {
public:
    TrialFunction(std::shared_ptr<FESpace> space, const std::string& name = "u")
        : Argument(std::move(space), ArgumentType::Trial, name) {}
        
    bool dependsOnFEField() const override { return false; }
    bool dependsOnShape() const override { return true; }
    
    ShapeEvaluationFlags requiredShapeFlags() const override {
        return ShapeEvaluationFlags::Value;
    }
};

class TestFunction : public Argument {
public:
    TestFunction(std::shared_ptr<FESpace> space, const std::string& name = "v")
        : Argument(std::move(space), ArgumentType::Test, name) {}

    bool dependsOnFEField() const override { return false; }
    bool dependsOnShape() const override { return true; }
    
    ShapeEvaluationFlags requiredShapeFlags() const override {
        return ShapeEvaluationFlags::Value;
    }
};

}