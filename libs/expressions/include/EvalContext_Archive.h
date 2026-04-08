#pragma once

#include "cfemutils.h"
#include "Vector3D.h"
#include "ReferenceElement.h"
#include "DynamicVector.h"

#include <deque>
#include <vector>
#include <span>

namespace cfem {

    /**
     * @brief Entry for the linear cache. 
     * id: The memory address of the expression.
     * data: A span pointing to the results stored in the EvalContext pages.
     */
    struct CacheEntry {
        const void* id;
        std::span<CfemReal> data;
    };

    struct EvalContext {
        Vector3D physPos = {};
        CfemReal currentTime = 0;
        const ShapeInfo* shapeP1 = nullptr;
        const ShapeInfo* shapeP2 = nullptr;
        const MappingInfo* mappingP2 = nullptr; 
        const MappingInfo* mappingP1 = nullptr; 
        
        CfemInt currentCellId = -1;
        Vector3D currentXi = {-1e9, -1e9, -1e9};

        const ShapeInfo* trialShape = nullptr;
        const ShapeInfo* testShape  = nullptr;

        const MappingInfo* trialMapping = nullptr;
        const MappingInfo* testMapping  = nullptr;

        CfemInt trialDof = -1;
        CfemInt testDof  = -1;

        ShapeInfo scratchShape;
        MappingInfo scratchMapping;
        DynamicVector<Vector3D> scratchCoords;
        std::vector<CfemInt> scratchDofs;

    private:
        static constexpr size_t PAGE_SIZE = 4096; 
        std::deque<std::vector<CfemReal>> m_pages; 
        size_t m_currentPageIdx = 0;
        size_t m_currentOffset = 0;

        // HIGH PERFORMANCE CACHE: Linear search is faster than a Map for small N.
        std::vector<CacheEntry> m_expressionCache;

    public:
        EvalContext() {
            m_pages.emplace_back(PAGE_SIZE);
            m_expressionCache.reserve(16); // Pre-reserve to avoid reallocations
        }

        std::span<CfemReal> getBuffer(CfemInt size) {
            if (size > PAGE_SIZE) {
                m_pages.insert(m_pages.begin() + m_currentPageIdx, std::vector<CfemReal>(size));
                auto& specialPage = m_pages[m_currentPageIdx];
                m_currentPageIdx++; 
                return std::span<CfemReal>(specialPage.data(), size);
            }

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

        /**
         * @brief Check if expression result is already in cache for this point.
         */
        bool tryGetCache(const void* exprPtr, std::span<CfemReal> out) const {
            for (const auto& entry : m_expressionCache) {
                if (entry.id == exprPtr) {
                    std::copy(entry.data.begin(), entry.data.end(), out.begin());
                    return true;
                }
            }
            return false;
        }

        /**
         * @brief Store result in the fast page buffer and index it in the linear cache.
         */
        void storeInCache(const void* exprPtr, std::span<CfemReal> result) {
            auto cacheStorage = getBuffer(static_cast<CfemInt>(result.size()));
            std::copy(result.begin(), result.end(), cacheStorage.begin());
            m_expressionCache.push_back({exprPtr, cacheStorage});
        }

        void resetBuffers() {
            m_currentPageIdx = 0;
            m_currentOffset = 0;
            m_expressionCache.clear(); // Resets size to 0 but keeps capacity
        }

        struct Bookmark { size_t page; size_t offset; };
        Bookmark getBookmark() const { return { m_currentPageIdx, m_currentOffset }; }
        void releaseToBookmark(Bookmark bm) {
            m_currentPageIdx = bm.page;
            m_currentOffset = bm.offset;
        }
    };
}