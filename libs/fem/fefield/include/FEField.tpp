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
#include "Mesh.h"

namespace cfem
{
    template <typename Func>
    requires (!std::is_base_of_v<sym::CfemExpression, std::decay_t<Func>> && 
              !std::is_convertible_v<std::decay_t<Func>, std::shared_ptr<const sym::CfemExpression>>)
    void FEField::project(const Func& func)
    {
        auto mesh = m_space->getMesh();
        int dim = getDimension();
        
        // Use a "visited" array only if the space is Continuous (CG)
        // In DG, we want to write every local DOF regardless of others.
        bool isContinuous = (m_space->getContinuity() == FESpace::Continuity::CG);
        std::vector<uint8_t> visited;
        if (isContinuous) {
            visited.assign(m_space->getNumUniquePositions(), 0);
        }

        // Cell-centric loop: Works for both CG and DG
        for (CfemInt cellId = 0; cellId < mesh->getNumCells(); ++cellId)
        {
            const auto& cell = mesh->getCell(cellId);
            const auto& refEl = m_space->getRefElement(cell.type);
            int nNodes = refEl.getNumNodes();

            for (int i = 0; i < nNodes; ++i)
            {
                // In CG, posId is the Vertex ID. In DG, posId is unique per cell/node pair.
                CfemInt posId = m_space->getPositionId(cellId, i);
                
                if (isContinuous) {
                    if (visited[posId]) continue;
                    visited[posId] = 1;
                }

                // Get physical coordinates of this local node
                // Note: We evaluate the node position in the physical world
                Vector3D physPtr = mesh->getVertex( mesh->getCellConnectivity(cellId)[i] ).coords;

                auto result = func(physPtr);

                for (int d = 0; d < dim; ++d)
                {
                    CfemInt globalIdx = m_space->getGlobalDofIdx(cellId, i, d);
                    if (globalIdx != -1) {
                        m_data.setValue(globalIdx, (dim > 1) ? result[d] : static_cast<CfemReal>(result), INSERT_VALUES);
                    }
                }
            }
        }

        m_data.assemble();
        CFEM_INFO << "Universal projection onto '" << m_name << "' completed.";
    }
}