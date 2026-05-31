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

#ifndef CFEM_FEFIELD_H_2
#define CFEM_FEFIELD_H_2

#include "cfemutils.h"
#include "CfemExpression.h"
#include "Argument.h"
#include "GradExpression.h"
#include "EvalContext.h"
#include "FESpace.h"
#include "PetscVector.h"
#include <string>
#include <vector>
#include <memory>
#include <set>
#include <span>

#include "Vector2D.h"
#include "Vector3D.h"
#include "Tensor2D.h"
#include "Tensor3D.h"

namespace cfem
{
    class FEField : public sym::TrialFunction, public FESpaceObserver, public std::enable_shared_from_this<FEField>
    {
    private:
        struct ConstructorToken {};
        PetscVector m_data;
        mutable std::vector<CfemReal> m_localValues;

        mutable CfemInt m_boundIndex = -1;

        std::shared_ptr<const CfemExpression> createTimeDerivative() const override {
            CFEM_ERROR << "Time derivative of a FEField is not yet implemented.";
            return nullptr;
        }

        std::shared_ptr<const sym::CfemExpression> createGradient() const override {
            return std::make_shared<sym::GradExpression>(shared_from_this());
        }

        void evaluateWithShape(CfemInt cellId, const ShapeInfo& shape, std::span<CfemReal> out) const;
        void evaluateWithoutShape(CfemInt cellId, const Vector3D& xi, std::span<CfemReal> out, EvalContext &ctx) const;
        void evaluateBatchWithShape(CfemInt cellId, CfemInt nQp, const ShapeInfoBatch& shapeBatch, std::span<CfemReal> outValues, EvalContext& ctx) const;
        void evaluateBatchWithoutShape(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const;


        void evaluateGradientWithCachedMapping(CfemInt cellId, const MappingInfo& mapping, std::span<CfemReal> out) const;
        void evaluateGradientWithoutCachedMapping(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const;
        void evaluateGradientBatchWithCachedMapping(CfemInt cellId, CfemInt nQp, const MappingInfoBatchView& mappingBatch, std::span<CfemReal> outValues, EvalContext& ctx) const;
        void evaluateGradientBatchWithoutCachedMapping(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const;

    public:
        FEField(ConstructorToken, std::shared_ptr<FESpace> space, const std::string &name);
        virtual ~FEField() = default;

        static std::shared_ptr<FEField> create(const std::string &name, std::shared_ptr<FESpace> space) {
            auto field = std::make_shared<FEField>(ConstructorToken{}, space, name);
            space->registerField(field);
            return field;
        }

        PetscVector &getData() { return m_data; }
        const PetscVector &getData() const { return m_data; }
        bool dependsOnFEField() const override { return true; }
        bool dependsOnShape() const override { return true; }

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
        
        // Data Management
        void updateOnSpaceChange() override { reallocate(); }
        void reallocate();
        void updateLocalBuffer() const;
        CfemReal getAtNode(CfemInt meshNodeId, CfemInt component = 0) const;

        CfemReal        getScalarFieldAtNode         (CfemInt meshNodeId) const;
        Vector2D        get2DVectorFieldAtNode       (CfemInt meshNodeId) const;
        Tensor2D        get2DTensorFieldAtNode       (CfemInt meshNodeId) const;
        SymTensor2D     get2DSymTensorFieldAtNode    (CfemInt meshNodeId) const;
        Vector3D        get3DVectorFieldAtNode       (CfemInt meshNodeId) const;
        Tensor3D        get3DTensorFieldAtNode       (CfemInt meshNodeId) const;
        SymTensor3D     get3DSymTensorFieldAtNode    (CfemInt meshNodeId) const;
        SkewSymTensor3D get3DSkewSymTensorFieldAtNode(CfemInt meshNodeId) const;

        void extractCellDofs(CfemInt cellId, DynamicMatrix<CfemReal>& localDofs, EvalContext& ctx) const;
        void evaluateGradient(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override;

        // ====================================================================
        // --- Vectorized Batch (SoA) Evaluation Extensions
        // ====================================================================
        void evaluateBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const override;
        void evaluateGradientBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const override;

        void interpolate(const sym::CfemExpression& recipe);
        void interpolate(std::shared_ptr<const sym::CfemExpression> recipe) {
            if (recipe) interpolate(*recipe);
            else CFEM_ERROR << "interpolate requires a valid CfemExpression object.";
        }

        template <typename Func>
        requires (!std::is_base_of_v<sym::CfemExpression, std::decay_t<Func>> && 
                  !std::is_convertible_v<std::decay_t<Func>, std::shared_ptr<const sym::CfemExpression>>)
        void interpolate(const Func& func);

        sym::ExpressionType getType() const override { return sym::ExpressionType::FEField; }
        std::string toString() const override { return "FEField::" + m_name; }
    };
} 

#include "FEField.tpp"
#endif