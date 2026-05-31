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
        CFEM_INFO << "Interpolating functional lambda/constant into field '" << getName() << "'";

        const auto& mesh = m_space->getMesh();
        Vec& petscVec = m_data.getRaw();
        const CfemInt nComp = this->getNumComponents();

        // Compile-time type analysis of the input parameter
        constexpr bool is_constant = std::is_arithmetic_v<std::decay_t<Func>>;
        constexpr bool is_callable = std::is_invocable_v<Func, Vector3D>;
        
        static_assert(is_constant || is_callable, 
            "FEField::interpolate requires either a numeric constant or a callable taking a Vector3D.");

        // Retrieve local partition layout boundaries for the parallel MPI distributed vector
        CfemInt low = 0;
        CfemInt high = 0;
        VecGetOwnershipRange(petscVec, &low, &high);
        m_localValues.assign(high - low, 0.0);

        const CfemInt numPositions = m_space->getNumUniquePositions();
        std::vector<uint8_t> visited(numPositions, 0);

        #pragma omp parallel
        {
            // Allocate thread-local context workspace safely isolating scratch matrix buffers
            EvalContext evctx;
            
            #pragma omp for schedule(guided)
            for (CfemInt cellId = 0; cellId < mesh->getNumCells(); ++cellId)
            {
                const auto cell_type = mesh->getCell(cellId).type;
                const auto& refEl = ReferenceElementFactory::getReferenceElement(cell_type, m_space->getOrderInt());
                const auto& refElGeo = ReferenceElementFactory::getReferenceElement(cell_type, mesh->getOrderInt());

                // Load contiguous cell node physical coordinates into context scratchpad
                mesh->getCellVerticesCoords(cellId, evctx.scratchCoords);

                const CfemInt nNodes = static_cast<CfemInt>(refEl.getNumNodes());
                for (CfemInt i = 0; i < nNodes; ++i)
                {
                    const CfemInt posId = m_space->getPositionId(cellId, i);
                    if (posId == -1) continue;

                    // Atomic boundary test-and-set to guarantee unique thread evaluation per geometric node
                    if (visited[posId]) continue;
                    if (__atomic_test_and_set(&visited[posId], __ATOMIC_RELAXED)) continue;

                    if constexpr (is_constant) 
                    {
                        for (CfemInt d = 0; d < nComp; ++d)
                        {
                            const CfemInt globalIdx = m_space->getGlobalDofIdx(cellId, i, d);
                            if (globalIdx >= low && globalIdx < high) 
                            {
                                m_localValues[globalIdx - low] = static_cast<CfemReal>(func);
                            }
                        }
                    } 
                    else 
                    {
                        // STANDARD PATH: Callable Analytic Function 
                        const Vector3D xi = refEl.getNodeCoords(i);

                        refElGeo.computeMapping(
                            xi, 
                            evctx.scratchCoords,
                            ShapeEvaluationFlags::Value,
                            MappingEvaluationFlags::PhysicalPosition, 
                            evctx.singleWS.shape,
                            evctx.singleWS.mapping
                        );

                        // Evaluate user-provided function lambda at computed physical spatial coordinates
                        auto result = func(evctx.singleWS.mapping.physicalPosition);
                        
                        using LocalResultType = decltype(result);
                        constexpr bool returns_scalar = std::is_arithmetic_v<LocalResultType>;

                        // Map the calculated physical properties to the local distributed ownership buffer
                        for (CfemInt d = 0; d < nComp; ++d)
                        {
                            const CfemInt globalIdx = m_space->getGlobalDofIdx(cellId, i, d);
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

        // Unify distributed MPI vectors into the global parallel algebraic architecture
        PetscScalar *array = nullptr;
        VecGetArray(petscVec, &array);
        std::copy(m_localValues.begin(), m_localValues.end(), array);
        VecRestoreArray(petscVec, &array);
        
        m_data.assemble();
    }
}