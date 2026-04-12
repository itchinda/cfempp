#pragma once

#include "cfemutils.h"
#include "CfemExpression.h"
#include "FESpace.h" 
#include "EvalContext.h"
#include <memory>
#include <string>
#include <map>

namespace cfem::sym {

enum class ArgumentType { Trial, Test };

/**
 * @class Argument
 * @brief Représente une fonction de base Trial (u) ou Test (v).
 * * Cette classe utilise un index (m_boundIndex) résolu lors de la phase de préparation (bind)
 * pour accéder instantanément aux données de l'espace dans l'EvalContext.
 */
class Argument : public CfemExpression {
protected:
    std::shared_ptr<FESpace> m_space;
    ArgumentType m_argType;
    std::string m_name;
    
    // Index direct vers EvalContext::m_spaceDataCache[m_boundIndex]
    // mutable car bind() est souvent appelé sur des expressions const.
    mutable int m_boundIndex = -1;

    virtual std::shared_ptr<const CfemExpression> createGradient() const override {
        CFEM_ERROR << "You should not try to call grad on Argument expression";
        return nullptr; // Les fonctions de base sont indépendantes du temps
    }

    virtual std::shared_ptr<const CfemExpression> createTimeDerivative() const override {
        CFEM_ERROR << "You should not try to call time derivative on Argument expression";
        return nullptr; // Les fonctions de base sont indépendantes du temps
    }

public:
    Argument(std::shared_ptr<FESpace> space, ArgumentType type, const std::string& name)
        : m_space(std::move(space)), m_argType(type), m_name(name) {}

    virtual ~Argument() = default;

    // --- Identification ---
    const std::string &getName() const { return m_name; }
    const std::shared_ptr<FESpace>& getSpace() const { return m_space; }
    ArgumentType getArgumentType() const { return m_argType; }
    
    virtual CfemInt getDimension() const override { return m_space->getDimension(); }
    
    virtual ExpressionType getType() const override {
        return (m_argType == ArgumentType::Trial) ? ExpressionType::Trial : ExpressionType::Test;
    }

    virtual bool isContinuous() const override {
        return m_space->getContinuity() == FESpace::Continuity::CG;
    }
    // --- Phase de Préparation (Binding) ---

    virtual void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
        spaces.insert(m_space);
    }

    /**
     * @brief Lie l'argument à un index spécifique dans le tableau de cache du contexte.
     * Appelé une seule fois par le DomainIntegrator avant la boucle de quadrature.
     */
    virtual void bind(const std::map<const FESpace*, int>& spaceToId) const override {
        auto it = spaceToId.find(m_space.get());
        if (it != spaceToId.end()) {
            m_boundIndex = it->second;
        } else {
            m_boundIndex = -1;
            CFEM_ERROR << "Argument (" << m_name << "): Space not registered in integrator.";
        }
    }

    // --- Évaluation Symbolique (Standard) ---
    // Un argument pur n'a pas de valeur de champ, il retourne 0.
    virtual void evaluate(CfemInt cellId, const Vector3D& xi, std::span<CfemReal> out, EvalContext& ctx) const override {
        std::fill(out.begin(), out.end(), 0.0);
    }

    // --- Évaluation pour Assemblage (Bases) ---

    /**
     * @brief Évalue la valeur de la fonction de base pour un DOF local donné.
     * @param localDofIdx Index local au sein de la cellule (0 à DofsPerCell-1).
     */
    virtual void evaluateBasis(CfemInt localDofIdx, std::span<CfemReal> out, EvalContext& ctx) const {
        // Accès O(1) via l'index lié
        CFEM_ASSERT(m_boundIndex >= 0, "Argument not bound before evaluation!");
        const auto& data = ctx.getSpaceData(m_boundIndex);

        const ShapeInfo* shape = data.shape;

        CFEM_ASSERT(shape != nullptr, "ShapeInfo is null in EvalContext!");
        CFEM_ASSERT(out.size() >= m_space->getDimension(), "Output span too small!");

        const int dim = m_space->getDimension();
        // Layout: Interleaved (N0_x, N0_y, N0_z, N1_x, N1_y...)
        const int node = static_cast<int>(localDofIdx / dim);
        const int comp = static_cast<int>(localDofIdx % dim);

        std::fill(out.begin(), out.end(), 0.0);
        out[comp] = shape->values[node];
    }

    /**
     * @brief Évalue le gradient physique de la fonction de base pour un DOF local donné.
     */
    virtual void evaluateBasisGradient(CfemInt localDofIdx, std::span<CfemReal> out, EvalContext& ctx) const {
        const auto& data = ctx.getSpaceData(m_boundIndex);
        const MappingInfo* mapping = data.mapping;

        const int dim = m_space->getDimension();
        const int node = static_cast<int>(localDofIdx / dim);
        const int comp = static_cast<int>(localDofIdx % dim);

        const auto& g = mapping->physicalShapeGradients[node];

        std::fill(out.begin(), out.end(), 0.0);
        // Le span de sortie contient le gradient complet (3 composantes par dimension d'espace)
        // ex pour vecteur 3D: [dxUx, dyUx, dzUx, dxUy, dyUy, dzUy, dxUz, dyUz, dzUz]
        out[comp * 3 + 0] = g.x;
        out[comp * 3 + 1] = g.y;
        out[comp * 3 + 2] = g.z;
    }

    virtual std::string toString() const override { return m_name; }
};

/**
 * @brief Spécialisation pour les fonctions tests (v)
 */
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

/**
 * @brief Spécialisation pour les fonctions inconnues / Trial (u)
 */
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

} // namespace cfem::sym