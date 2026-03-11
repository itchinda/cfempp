#pragma once

#include "cfemutils.h"
#include "Vector3D.h"
#include "TensorO2.h"
#include "ReferenceElement.h"
#include "DynamicVector.h"

#include <deque>
#include <vector>
#include <span>

namespace cfem{
    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    struct EvalContext {
        Vector3D physPos;
        CfemReal currentTime;
        const ShapeInfo* shapeP1 = nullptr;
        const ShapeInfo* shapeP2 = nullptr;
        
        const MappingInfo* mappingP2 = nullptr; 
        const MappingInfo* mappingP1 = nullptr; 
        
        CfemInt currentCellId = -1;
        Vector3D currentXi = {-1e9, -1e9, -1e9};

        ShapeInfo scratchShape;
        MappingInfo scratchMapping;
        DynamicVector<Vector3D> scratchCoords;

private:
        static constexpr size_t PAGE_SIZE = 4096; // 32 Ko par page (4096 * 8 octets)
        std::deque<std::vector<CfemReal>> m_pages; 
        size_t m_currentPageIdx = 0;
        size_t m_currentOffset = 0;

    public:
        EvalContext() {
            m_pages.emplace_back(PAGE_SIZE);
        }

        std::span<CfemReal> getBuffer(CfemInt size) {
            // Si la demande dépasse la taille d'une page standard, on crée une page dédiée
            if (size > PAGE_SIZE) {
                m_pages.insert(m_pages.begin() + m_currentPageIdx, std::vector<CfemReal>(size));
                auto& specialPage = m_pages[m_currentPageIdx];
                m_currentPageIdx++; // On "saute" cette page pour la suite
                return std::span<CfemReal>(specialPage.data(), size);
            }

            // Si on dépasse la page actuelle, on passe à la suivante
            if (m_currentOffset + size > PAGE_SIZE) {
                m_currentPageIdx++;
                m_currentOffset = 0;
                if (m_currentPageIdx >= m_pages.size()) {
                    m_pages.emplace_back(PAGE_SIZE);
                }
            }

            auto span = std::span<CfemReal>(m_pages[m_currentPageIdx].data() + m_currentOffset, size);
            m_currentOffset += size;
            return span;
        }

        struct Bookmark { size_t page; size_t offset; };

        Bookmark getBookmark() const {
            return { m_currentPageIdx, m_currentOffset };
        }

        void releaseToBookmark(Bookmark bm) {
            m_currentPageIdx = bm.page;
            m_currentOffset = bm.offset;
        }

        void resetBuffers() {
            m_currentPageIdx = 0;
            m_currentOffset = 0;
        }
    
    };
}