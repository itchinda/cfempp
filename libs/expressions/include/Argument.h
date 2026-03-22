#pragma once

#include "CfemExpression.h"
#include "LagrangeSpace.h" // Ensure this is included or forward declared
#include "EvalContext.h"
#include <memory>
#include <string>

namespace cfem::sym {

enum class ArgumentType { Trial, Test };

class Argument : public CfemExpression {
protected:
    std::shared_ptr<LagrangeSpace> m_space;
    ArgumentType m_argType;
    std::string m_name;

public:
    Argument(std::shared_ptr<LagrangeSpace> space, ArgumentType type, const std::string& name)
        : m_space(std::move(space)), m_argType(type), m_name(name) {}

    virtual ~Argument() = default;

    const std::string &getName() const { return m_name; }
    const std::shared_ptr<LagrangeSpace>& getSpace() const { return m_space; }
    ArgumentType getArgumentType() const { return m_argType; }

    CfemInt getDimension() const override { return m_space->getDimension(); }
    
    virtual ExpressionType getType() const override {
        return (m_argType == ArgumentType::Trial) ? ExpressionType::Trial : ExpressionType::Test;
    }

    std::string toString() const override { return m_name; }

    std::shared_ptr<const CfemExpression> timeDerivative() const override {
        CFEM_ERROR << "Time derivative not implemented for Argument (" << m_name << ").";
        return nullptr;
    }

    // --- Standard Evaluation (Returns 0 for pure symbolic arguments) ---
    void evaluate(CfemInt cellId, const Vector3D& xi, std::span<CfemReal> out, EvalContext& ctx) const override {
        std::fill(out.begin(), out.end(), 0.0);
    }

    // --- The Explicit Basis Interface (Used by the Assembler) ---
    virtual void evaluateBasis(CfemInt dof, std::span<CfemReal> out, EvalContext& ctx) const {
        const ShapeInfo* shape = (m_argType == ArgumentType::Trial) ? ctx.trialShape : ctx.testShape;
        if (!shape || dof < 0) {
            CFEM_ERROR << "Argument (" << m_name << "): missing shape context.";
        }

        const CfemInt dim = m_space->getDimension();
        const CfemInt numNodes = m_space->getNumActiveNodes();
        const CfemInt comp = dof / numNodes;
        const CfemInt node = dof % numNodes;

        std::fill(out.begin(), out.end(), 0.0);
        out[comp] = shape->values[node];
    }

    virtual void evaluateBasisGradient(CfemInt dof, std::span<CfemReal> out, EvalContext& ctx) const {
        const MappingInfo* mapping = (m_argType == ArgumentType::Trial) ? ctx.trialMapping : ctx.testMapping;
        if (!mapping || dof < 0) {
            CFEM_ERROR << "Argument (" << m_name << "): missing mapping context.";
        }

        const CfemInt numNodes = m_space->getNumActiveNodes();
        const CfemInt comp = dof / numNodes;
        const CfemInt node = dof % numNodes;

        const auto& g = mapping->physicalShapeGradients[node];

        std::fill(out.begin(), out.end(), 0.0);
        out[comp * 3 + 0] = g.x;
        out[comp * 3 + 1] = g.y;
        out[comp * 3 + 2] = g.z;
    }
};

// Explicit subclasses so FEField can inherit cleanly
class TrialFunction : public Argument {
public:
    TrialFunction(std::shared_ptr<LagrangeSpace> space, const std::string& name = "u")
        : Argument(std::move(space), ArgumentType::Trial, name) {}
    virtual ExpressionType getType() const override { return ExpressionType::Trial; }
};

class TestFunction : public Argument {
public:
    TestFunction(std::shared_ptr<LagrangeSpace> space, std::string name = "v")
        : Argument(std::move(space), ArgumentType::Test, std::move(name)) {}

    virtual ExpressionType getType() const override { return ExpressionType::Test; }
};

} // namespace cfem::sym