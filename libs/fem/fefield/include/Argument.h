//************************************************************************
// --- CFEM++ Library
// ---
// --- Copyright 2024-2026 Ismaël Tchinda Ngueyong et al.
// --- ALL RIGHTS RESERVED
// ---
// --- This software is protected by international copyright laws.
// --- Use, distribution, or modification of this software in any form,
// --- source or binary, for personal, academic, or non-commercial 
// --- purposes is permitted free of charge, provided that this 
// --- copyright notice and this permission notice appear in all 
// --- copies and supporting documentation.
// ---
// --- The authors and contributors provide this code "as is" without 
// --- any express or implied warranty. In no event shall they be 
// --- held liable for any damages arising from the use of this software.
//************************************************************************

#pragma once

#include "cfemutils.h"
#include "CfemExpression.h"
#include "FESpace.h" 
#include "EvalContext.h"
#include <memory>
#include <string>
#include <map>

namespace cfem::sym {

/**
 * @enum ArgumentType
 * @brief Defines the role of an argument within a variational formulation.
 */
enum class ArgumentType { Trial, Test };

/**
 * @class Argument
 * @brief Represents a basis function, either Trial (u) or Test (v).
 * This class uses a resolved index (m_boundIndex) during the binding phase
 * to allow O(1) access to space data cached within the EvalContext.
 */
class Argument : public CfemExpression {
protected:

    std::shared_ptr<FESpace> m_space; ///< The Finite Element Space associated with this argument.
    ArgumentType m_argType;           ///< The type of argument (Trial or Test).
    std::string m_name;               ///< The symbolic name (e.g., "u", "v").
    
    /** 
     * @brief Direct index to EvalContext::m_spaceDataCache[m_boundIndex].
     * @note Marked as mutable because bind() is often called on const expressions.
     */
    mutable CfemInt m_boundIndex = -1;

    /**
     * @brief Prevents direct gradient calculation on a raw Argument.
     * @return std::shared_ptr<const CfemExpression> Always returns nullptr.
     * @note Triggers an internal error as basis functions are independent in this scope.
     */
    virtual std::shared_ptr<const CfemExpression> createGradient() const override {
        CFEM_ERROR << "You should not try to call grad on Argument expression";
        return nullptr; // Les fonctions de base sont indépendantes du temps
    }

    /**
     * @brief Prevents direct time derivative calculation on a raw Argument.
     * @return std::shared_ptr<const CfemExpression> Always returns nullptr.
     * @note Triggers an internal error as basis functions are time-independent.
     */
    virtual std::shared_ptr<const CfemExpression> createTimeDerivative() const override {
        CFEM_ERROR << "You should not try to call time derivative on Argument expression";
        return nullptr;
    }

public:

    /**
     * @brief Constructor for the Argument class.
     * @param space The FE Space this argument belongs to.
     * @param type Whether it is a Test or Trial function.
     * @param name The symbolic name for identification.
     */
    Argument(std::shared_ptr<FESpace> space, ArgumentType type, const std::string& name)
        : m_space(std::move(space)), m_argType(type), m_name(name) {}

    /**
     * @brief Virtual destructor.
     */
    virtual ~Argument() = default;

    /**
     * @brief Gets the symbolic name of the argument.
     * @return const std::string& The name string.
     */
    const std::string &getName() const { return m_name; }

    /**
     * @brief Gets the associated Finite Element Space.
     * @return const std::shared_ptr<FESpace>& Pointer to the FESpace.
     */
    const std::shared_ptr<FESpace>& getSpace() const { return m_space; }

    /**
     * @brief Gets the argument role (Trial or Test).
     * @return ArgumentType The type enum.
     */
    ArgumentType getArgumentType() const { return m_argType; }
    
    /**
     * @brief Gets the vector dimension of the associated space.
     * @return CfemInt The dimension (1 for scalars, >1 for vectors).
     */
    virtual CfemInt getNumComponents() const override { return m_space->getNumComponents(); }
    
    /**
     * @brief Gets the symbolic expression type.
     * @return ExpressionType Trial or Test depending on the argument role.
     */
    virtual ExpressionType getType() const override {
        return (m_argType == ArgumentType::Trial) ? ExpressionType::Trial : ExpressionType::Test;
    }

    /**
     * @brief Checks if the associated space is continuous (e.g., CG elements).
     * @return True if the space is continuous, false otherwise.
     */
    virtual bool isContinuous() const override {
        return m_space->getContinuity() == FESpace::Continuity::CG;
    }
    
    // --- Preparation Phase (Binding) ---

    /**
     * @brief Collects the FE spaces required by this expression.
     * @param spaces A set to store unique pointers to required FE spaces.
     */
    virtual void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
        spaces.insert(m_space);
    }

    /**
     * @brief Binds the argument to a specific index in the context cache.
     * @note Typically called once by the Integrator/Assembler before the quadrature loop.
     * @param spaceToId A map linking FE space pointers to their respective cache IDs.
     */
    virtual void bind(const std::map<const FESpace*, int>& spaceToId) const override {
        auto it = spaceToId.find(m_space.get());
        if (it != spaceToId.end()) {
            m_boundIndex = it->second;
        } else {
            m_boundIndex = -1;
            CFEM_ERROR << "Argument (" << m_name << "): Space not registered.";
        }
    }

    /**
     * @brief Standard symbolic evaluation.
     * @note An Argument cannot be evaluated as a stand-alone field because it represents 
     * a basis function. Use evaluateBasis() during assembly instead.
     * @throws std::runtime_error Always, to prevent invalid symbolic usage.
     */
    virtual void evaluate(CfemInt cellId, const Vector3D& xi, std::span<CfemReal> out, EvalContext& ctx) const override 
    {
        CFEM_ERROR << "Invalid call: Argument '" << m_name 
                   << "' cannot be evaluated as a field. "
                   << "Arguments are placeholders for basis functions and can only be "
                   << "evaluated via evaluateBasis() during assembly.";
    }

    /**
     * @brief Evaluates the basis function value for a specific local Degree of Freedom (DOF).
     * @param localDofIdx Local index within the cell (0 to DofsPerCell-1).
     * @param out Output span (size must be >= space dimension).
     * @param ctx Evaluation context containing precomputed shape data.
     */
    virtual void evaluateBasis(CfemInt localDofIdx, std::span<CfemReal> out, EvalContext& ctx) const {
        // O(1) Access via bound index
        CFEM_ASSERT(m_boundIndex >= 0, "Argument not bound before evaluation!");
        const auto& data = ctx.getSpaceData(m_boundIndex);

        const ShapeInfo* shape = data.shape;

        CFEM_ASSERT(shape != nullptr, "ShapeInfo is null in EvalContext!");
        CFEM_ASSERT(out.size() >= m_space->getNumComponents(), "Output span too small!");

        const CfemInt dim = m_space->getNumComponents();
        // Layout: Interleaved (N0_x, N0_y, N0_z, N1_x, N1_y...)
        const CfemInt node = static_cast<int>(localDofIdx / dim);
        const CfemInt comp = static_cast<int>(localDofIdx % dim);

        std::fill(out.begin(), out.end(), 0.0);
        out[comp] = shape->values[node];
    }

    /**
     * @brief Evaluates the physical gradient of the basis function for a specific local DOF.
     * @param localDofIdx Local index within the cell.
     * @param out Output span containing the full gradient (3 components per space dimension).
     * @param ctx Evaluation context containing precomputed mapping data.
     */
    virtual void evaluateBasisGradient(CfemInt localDofIdx, std::span<CfemReal> out, EvalContext& ctx) const {
        const auto& data = ctx.getSpaceData(m_boundIndex);
        const MappingInfo* mapping = data.mapping;

        const CfemInt dim = m_space->getNumComponents();
        const CfemInt node = static_cast<int>(localDofIdx / dim);
        const CfemInt comp = static_cast<int>(localDofIdx % dim);

        const auto& g = mapping->physicalShapeGradients[node];

        std::fill(out.begin(), out.end(), 0.0);
        // Output span contains full gradient: [dxUx, dyUx, dzUx, dxUy, dyUy, dzUy, dxUz, dyUz, dzUz]
        out[comp * 3 + 0] = g.x;
        out[comp * 3 + 1] = g.y;
        out[comp * 3 + 2] = g.z;
    }

    /**
     * @brief Converts the argument to a string representation.
     * @return std::string The name of the argument.
     */
    virtual std::string toString() const override { return m_name; }
};


/**
 * @class TestFunction
 * @brief Specialization of Argument for test functions (typically denoted 'v').
 */
class TestFunction : public Argument {
public:
    /**
     * @brief Constructor for TestFunction.
     * @param space The associated FE Space.
     * @param name Symbolic name (defaults to "v").
     */
    TestFunction(std::shared_ptr<FESpace> space, const std::string& name = "v")
        : Argument(std::move(space), ArgumentType::Test, name) {}

    /**
     * @brief Checks if this function depends on an FE Field value.
     * @return False (test functions are abstract variations).
     */
    bool dependsOnFEField() const override { return false; }

    /**
     * @brief Checks if this function requires shape function evaluation.
     * @return True.
     */
    bool dependsOnShape() const override { return true; }
    
    /**
     * @brief Returns required shape evaluation flags.
     * @return ShapeEvaluationFlags::Value.
     */
    ShapeEvaluationFlags requiredShapeFlags() const override {
        return ShapeEvaluationFlags::Value;
    }
};

/**
 * @class TrialFunction
 * @brief Specialization of Argument for unknown/trial functions (typically denoted 'u').
 */
class TrialFunction : public Argument {
public:
    /**
     * @brief Constructor for TrialFunction.
     * @param space The associated FE Space.
     * @param name Symbolic name (defaults to "u").
     */
    TrialFunction(std::shared_ptr<FESpace> space, const std::string& name = "u")
        : Argument(std::move(space), ArgumentType::Trial, name) {}
        
    /**
     * @brief Checks if this function depends on an FE Field value.
     * @return False (in the assembly context, TrialFunction is a basis, not an evaluated field).
     */
    bool dependsOnFEField() const override { return false; }

    /**
     * @brief Checks if this function requires shape function evaluation.
     * @return True.
     */
    bool dependsOnShape() const override { return true; }
    
    /**
     * @brief Returns required shape evaluation flags.
     * @return ShapeEvaluationFlags::Value.
     */
    ShapeEvaluationFlags requiredShapeFlags() const override {
        return ShapeEvaluationFlags::Value;
    }
};

} // namespace cfem::sym