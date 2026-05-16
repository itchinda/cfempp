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

#include "FEField.h"
#include "ReferenceElementFactory.h"
#include "Mesh.h"

#include <type_traits>

namespace cfem
{
    /**
     * @brief Projects a functional lambda onto the FE Space via nodal interpolation.
     * 
     * This method evaluates a user-provided function (e.g., a mathematical formula 
     * or analytic solution) at each unique DOF location and assigns the values 
     * to the field.
     * 
     * @tparam Func A callable type with signature (Vector3D) -> CfemReal or array-like.
     * @param func The function to evaluate at physical node positions.
     * @note This method is OpenMP-parallelized. It uses an atomic visiting map 
     * to ensure that shared nodes in Continuous Galerkin (CG) spaces are 
     * only evaluated once.
     */
    template <typename Func>
    requires (!std::is_base_of_v<sym::CfemExpression, std::decay_t<Func>> && 
              !std::is_convertible_v<std::decay_t<Func>, std::shared_ptr<const sym::CfemExpression>>)
    void FEField::interpolate(const Func& func)
    {
        CFEM_INFO << "Projecting functional lambda/constant into field '" << getName() << "'";

        const auto& mesh = m_space->getMesh();
        Vec& petscVec = m_data.getRaw();
        CfemInt nComp = this->getNumComponents();

        // Compile-time type analysis of the input
        constexpr bool is_constant = std::is_arithmetic_v<std::decay_t<Func>>;
        constexpr bool is_callable = std::is_invocable_v<Func, Vector3D>;
        static_assert(is_constant || is_callable, 
            "FEField::interpolate requires either a numeric constant or a callable taking a Vector3D.");

        // MPI Ownership Range
        CfemInt low, high;
        VecGetOwnershipRange(petscVec, &low, &high);
        m_localValues.assign(high - low, 0.0);

        CfemInt numPositions = m_space->getNumUniquePositions();
        std::vector<uint8_t> visited(numPositions, 0);

        #pragma omp parallel
        {
            EvalContext evctx;
            
            #pragma omp for schedule(guided)
            for (CfemInt cellId = 0; cellId < mesh->getNumCells(); ++cellId)
            {
                auto cell_type = mesh->getCell(cellId).type;
                const auto& refEl = ReferenceElementFactory::getReferenceElement(cell_type, m_space->getOrderInt());
                const auto& refElGeo = ReferenceElementFactory::getReferenceElement(cell_type, mesh->getOrderInt());

                // Load coordinates into EvalContext scratchpad
                mesh->getCellVerticesCoords(cellId, evctx.scratchCoords);

                for (CfemInt i = 0; i < (CfemInt)refEl.getNumNodes(); ++i)
                {
                    CfemInt posId = m_space->getPositionId(cellId, i);
                    if (posId == -1) continue;

                    // Atomic visiting check to handle shared nodes (CG)
                    if (visited[posId]) continue;
                    if (__atomic_test_and_set(&visited[posId], __ATOMIC_RELAXED)) continue;

                    if constexpr (is_constant) 
                    {
                        // --- FAST PATH: Constant Value ---
                        // No geometry mapping needed. Just apply the constant to the DOFs.
                        for (CfemInt d = 0; d < nComp; ++d)
                        {
                            CfemInt globalIdx = m_space->getGlobalDofIdx(cellId, i, d);
                            if (globalIdx >= low && globalIdx < high) 
                            {
                                m_localValues[globalIdx - low] = static_cast<CfemReal>(func);
                            }
                        }
                    } 
                    else 
                    {
                        // --- STANDARD PATH: Callable Lambda ---
                        // Precise Geometry: Find physical position of the node
                        Vector3D xi = refEl.getNodeCoords(i);
                        refElGeo.evaluate(xi, ShapeEvaluationFlags::Value, evctx.scratchShape);
                        refElGeo.computeMapping(
                            xi, 
                            evctx.scratchCoords, 
                            ShapeEvaluationFlags::PhysicalPosition, 
                            static_cast<const ShapeInfo&>(evctx.scratchShape), 
                            evctx.scratchMapping);

                        // Evaluate lambda
                        auto result = func(evctx.scratchMapping.physicalPosition);
                        
                        // Safely deduce return type (array vs scalar) inside the isolated branch
                        using LocalResultType = decltype(result);
                        constexpr bool returns_scalar = std::is_arithmetic_v<LocalResultType>;

                        // Global Assembly (Local Buffer)
                        for (CfemInt d = 0; d < nComp; ++d)
                        {
                            CfemInt globalIdx = m_space->getGlobalDofIdx(cellId, i, d);
                            if (globalIdx >= low && globalIdx < high) 
                            {
                                if constexpr (returns_scalar) {
                                    m_localValues[globalIdx - low] = static_cast<CfemReal>(result);
                                } else {
                                    m_localValues[globalIdx - low] = static_cast<CfemReal>(result[d]);
                                }
                            }
                        }
                    }
                }
            }
        }

        // Direct PETSc update
        PetscScalar *array;
        VecGetArray(petscVec, &array);
        std::copy(m_localValues.begin(), m_localValues.end(), array);
        VecRestoreArray(petscVec, &array);
        
        m_data.assemble();
    }
}