// //************************************************************************
// // --- CFEM++ Library
// // ---
// // --- Copyright 2024-2026 Ismaël Tchinda Ngueyong et al.
// // --- ALL RIGHTS RESERVED
// // ---
// // --- This software is protected by international copyright laws.
// // --- Use, distribution, or modification of this software in any form,
// // --- source or binary, for personal, academic, or non-commercial 
// // --- purposes is permitted free of charge, provided that this 
// // --- copyright notice and this permission notice appear in all 
// // --- copies and supporting documentation.
// // ---
// // --- The authors and contributors provide this code "as is" without 
// // --- any express or implied warranty. In no event shall they be 
// // --- held liable for any damages arising from the use of this software.
// //************************************************************************

// #pragma once

// #include "cfemutils.h"
// #include "Vector3D.h"
// #include "DynamicVector.h"
// #include "DynamicMatrix.h"

// #include "ShapeInfo.h"
// #include "ShapeInfoBatch.h"
// #include "MappingInfo.h"
// #include "MappingInfoBatchView.h"

// #include <vector>
// #include <span>
// #include <array>
// #include <deque>
// #include <iostream>
// #include <algorithm>

// namespace cfem {

//     class FESpace; // Forward declaration
    
//     // Forward declarations des structures Batch (définies dans votre architecture SoA)
//     struct ShapeInfoBatch;
//     struct MappingInfoBatch;
    
//     struct MappingInfoBatch;

//     /**
//      * @struct SpaceDataEntry
//      * @brief Ultra-fast link between an FE Space and its evaluations at the current integration cell.
//      * Uses explicit pointers to safely separate legacy/single-point data from vectorized batch data.
//      */
//     struct SpaceDataEntry {
//         const ShapeInfo* shapeSingle   = nullptr; 
//         const MappingInfo* mappingSingle = nullptr; 
        
//         const ShapeInfoBatch* shapeBatch    = nullptr; 
//         const MappingInfoBatchView* mappingBatch  = nullptr;
//     };

//     /**
//      * @struct CacheEntry
//      * @brief Stores a pointer-based identifier and its associated data span.
//      */
//     struct CacheEntry {
//         const void* id;           
//         std::span<CfemReal> data; 
//     };

//     /**
//      * @struct EvalContext
//      * @brief High-performance evaluation context for Finite Element assembly.
//      * Provides O(1) access to geometric data and an Arena Memory system.
//      * Structured to support both Single-Point (Probing) and Batch (Assembly) evaluation modes.
//      */
//     struct EvalContext {

//         // ====================================================================
//         // GLOBAL SIMULATION STATE
//         // ====================================================================
//         CfemReal currentTime   = 0.0; /**< Current simulation time (for transient problems). */
//         CfemInt  currentCellId = -1;  /**< Index of the cell currently being integrated. */
//         CfemInt currentQpIndex = 0;

//         // ====================================================================
//         // SINGLE-POINT WORKSPACE (Legacy / Probing)
//         // ====================================================================
//         struct SingleWorkspace {
//             Vector3D    physPos = {};                /**< Physical coordinates. */
//             Vector3D    xi      = {-1e9, -1e9, -1e9};/**< Reference coordinates. */
//             ShapeInfo   shape;                       /**< Buffer for base functions. */
//             MappingInfo mapping;                     /**< Buffer for geometric mapping. */
//         } singleWS;

//         // ====================================================================
//         // BATCH WORKSPACE (Vectorized Element Assembly)
//         // ====================================================================
//         struct BatchWorkspace {
//             DynamicVector<Vector3D>* physPos;     /**< Physical coordinates for all QPs. */
//             DynamicVector<Vector3D>* xi;          /**< Reference coordinates for all QPs. */
//             ShapeInfoBatch* shape   = nullptr;   /**< Pointer to allocated batch shape memory. */
//             MappingInfoBatch* mapping = nullptr; /**< Pointer to allocated batch metric memory. */
//         } batchWS;

//         // ====================================================================
//         // MULTI-SPACE CACHE & SHORTCUTS (O(1) Access)
//         // ====================================================================
//         static constexpr CfemInt MAX_CACHED_SPACES = 16;

//         #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
//             std::vector<SpaceDataEntry> m_spaceDataCache;
//         #else
//             std::array<SpaceDataEntry, MAX_CACHED_SPACES> m_spaceDataCache;
//         #endif

//         const ShapeInfoBatch* trialShapeBatch     = nullptr;
//         const ShapeInfoBatch* testShapeBatch      = nullptr;
//         const MappingInfoBatch* trialMappingBatch   = nullptr;
//         const MappingInfoBatch* testMappingBatch    = nullptr;
//         const MappingInfoBatch* geometricMappingBatch = nullptr;

        
//         const MappingInfo* geometricMapping = nullptr;

//         CfemInt trialDof = -1; 
//         CfemInt testDof  = -1; 

//         // ====================================================================
//         // PRE-ALLOCATED SCRATCHPADS
//         // ====================================================================
//         DynamicVector<Vector3D> scratchCoords;   /**< Buffer for cell vertex coordinates. */
//         std::vector<CfemInt>    scratchDofs;     /**< Buffer for gathering local DOF indices. */
//         DynamicMatrix<CfemReal> scratchCellDofs; /**< Buffer for cell DOF values. */
//         MappingInfo             scratchMapping;  /**< Buffer for cell DOF values. */
//         ShapeInfo               scratchShape;    /**< Buffer for cell DOF values. */
//         // ====================================================================
//         // ARENA MEMORY SYSTEM
//         // ====================================================================
//         static constexpr size_t PAGE_SIZE = 4096;

//         std::deque<std::vector<CfemReal>> m_pages;  
//         size_t m_currentPageIdx = 0;
//         size_t m_currentOffset   = 0;

//         std::vector<CacheEntry> m_expressionCache; 

//         // ====================================================================
//         // API METHODS
//         // ====================================================================

//         EvalContext() 
//         {
//             m_pages.emplace_back(PAGE_SIZE);
//             m_expressionCache.reserve(16);
// #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
//             m_spaceDataCache.reserve(8);
// #endif
//         }

// #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
//         inline void prepareSpaceCache(CfemInt n) {
//             m_spaceDataCache.resize(n);
//         }
// #endif

//         inline const SpaceDataEntry& getSpaceData(CfemInt index) const 
//         {
//         #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
//             CFEM_ASSERT(index >= 0 && index < static_cast<int>(m_spaceDataCache.size()), "EvalContext::getSpaceData index out of bounds");
//         #else
//             CFEM_ASSERT(index >= 0 && index < MAX_CACHED_SPACES, "EvalContext::getSpaceData index out of bounds");
//         #endif
//             return m_spaceDataCache[index];
//         }

//         inline void setSpaceData(CfemInt index, const ShapeInfo* shape, const MappingInfo* mapping) 
//         {
//             CFEM_ASSERT(index >= 0, "Negative space index");
// #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
//             CFEM_ASSERT(index < static_cast<int>(m_spaceDataCache.size()), "Space cache not pre-sized");
// #else
//             CFEM_ASSERT(index < MAX_CACHED_SPACES, "Space index exceeds MAX_CACHED_SPACES");
// #endif
//             auto& entry = m_spaceDataCache[index];
//             entry.shapeSingle = shape;
//             entry.mappingSingle = mapping;
//         }

//         inline void setSpaceDataBatch(CfemInt index, const ShapeInfoBatch* shape, const MappingInfoBatch* mapping) 
//         {
//             CFEM_ASSERT(index >= 0, "Negative space index");
// #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
//             CFEM_ASSERT(index < static_cast<int>(m_spaceDataCache.size()), "Space cache not pre-sized");
// #else
//             CFEM_ASSERT(index < MAX_CACHED_SPACES, "Space index exceeds MAX_CACHED_SPACES");
// #endif
//             auto& entry = m_spaceDataCache[index];
//             entry.shapeBatch = shape;
//             entry.mappingBatch = mapping;
//         }

//         inline const MappingInfoBatch* getMappingBatch(CfemInt index) const 
//         {
//         #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
//             CFEM_ASSERT(index >= 0 && index < static_cast<int>(m_spaceDataCache.size()), "EvalContext::getMappingBatch index out of bounds");
//         #else
//             CFEM_ASSERT(index >= 0 && index < MAX_CACHED_SPACES, "EvalContext::getMappingBatch index out of bounds");
//         #endif
//             const auto* mapping = m_spaceDataCache[index].mappingBatch;
//             CFEM_ASSERT(mapping != nullptr, "MappingInfoBatch not set");
//             return mapping;
//         }

//         inline const MappingInfo* getMapping(CfemInt index) const 
//         {
//         #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
//             CFEM_ASSERT(index >= 0 && index < static_cast<int>(m_spaceDataCache.size()), "EvalContext::getMappingSingle index out of bounds");
//         #else
//             CFEM_ASSERT(index >= 0 && index < MAX_CACHED_SPACES, "EvalContext::getMappingSingle index out of bounds");
//         #endif
//             const auto* mapping = m_spaceDataCache[index].mappingSingle;
//             CFEM_ASSERT(mapping != nullptr, "MappingInfo not set");
//             return mapping;
//         }

//         /**
//          * @brief Retrieves precomputed BATCH shape info for a given space index.
//          * @param index Space identifier.
//          * @return Pointer to ShapeInfoBatch.
//          * @pre index is valid and batch shape has been set.
//          */
//         inline const ShapeInfoBatch* getShapeBatch(CfemInt index) const 
//         {
//         #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
//             CFEM_ASSERT(index >= 0 && index < static_cast<int>(m_spaceDataCache.size()), "EvalContext::getShapeBatch index out of bounds");
//         #else
//             CFEM_ASSERT(index >= 0 && index < MAX_CACHED_SPACES, "EvalContext::getShapeBatch index out of bounds");
//         #endif
//             const auto* shape = m_spaceDataCache[index].shapeBatch;
//             CFEM_ASSERT(shape != nullptr, "ShapeInfoBatch not set. Did you call setSpaceDataBatch?");
//             return shape;
//         }

//         /**
//          * @brief Retrieves precomputed SINGLE shape info for a given space index.
//          * @param index Space identifier.
//          * @return Pointer to ShapeInfo (Legacy/Probing).
//          * @pre index is valid and single shape has been set.
//          */
//         inline const ShapeInfo* getShape(CfemInt index) const 
//         {
//         #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
//             CFEM_ASSERT(index >= 0 && index < static_cast<int>(m_spaceDataCache.size()), "EvalContext::getShapeSingle index out of bounds");
//         #else
//             CFEM_ASSERT(index >= 0 && index < MAX_CACHED_SPACES, "EvalContext::getShapeSingle index out of bounds");
//         #endif
//             const auto* shape = m_spaceDataCache[index].shapeSingle;
//             CFEM_ASSERT(shape != nullptr, "ShapeInfo not set. Did you call setSpaceDataSingle?");
//             return shape;
//         }

//         void resetBuffers() {
//             m_currentPageIdx = 0;
//             m_currentOffset = 0;
//             m_expressionCache.clear();
//         }

//         /**
//          * @brief Allocates a buffer from the Arena memory.
//          * @param size Number of CfemReal elements needed.
//          * @return A span pointing to the allocated memory.
//          */
//         std::span<CfemReal> getBuffer(CfemInt size) {
//             CFEM_ASSERT(size >= 0, "Negative buffer size");
//             CFEM_ASSERT(m_currentPageIdx < m_pages.size(),
//                         "Invalid arena state: page index out of bounds");

//             size_t usize = static_cast<size_t>(size);

//             if (m_currentOffset + usize > m_pages[m_currentPageIdx].size()){
//                 // Check if the current page has enough remaining space
//                 m_currentPageIdx++;
//                 m_currentOffset = 0;
                
//                 // If we ran out of pages, add a new one
//                 if (m_currentPageIdx >= m_pages.size()) {
//                     m_pages.emplace_back(std::max(usize, PAGE_SIZE));
//                 } 
//                 // If a page exists from a previous run but is too small, resize it
//                 else if (m_pages[m_currentPageIdx].size() < usize) {
//                     m_pages[m_currentPageIdx].resize(usize);
//                 }
//             }
            
//             auto span = std::span<CfemReal>(m_pages[m_currentPageIdx].data() + m_currentOffset, usize);
//             m_currentOffset += usize;
//             return span;
//         }

//         bool tryGetCache(const void* exprPtr, std::span<CfemReal> out) const {
//             for (const auto& entry : m_expressionCache) {
//                 if (entry.id == exprPtr) {
//                     CFEM_ASSERT(out.size() >= entry.data.size(), "Output buffer too small in tryGetCache");
//                     std::copy(entry.data.begin(), entry.data.end(), out.begin());
//                     return true;
//                 }
//             }
//             return false;
//         }

//         void storeInCache(const void* exprPtr, std::span<CfemReal> result) {
//             auto storage = getBuffer(static_cast<CfemInt>(result.size()));
//             std::copy(result.begin(), result.end(), storage.begin());
//             m_expressionCache.push_back({exprPtr, storage});
//         }

//         struct Bookmark { size_t page; size_t offset; size_t cacheSize; };
        
//         Bookmark getBookmark() const {
//             return { m_currentPageIdx, m_currentOffset, m_expressionCache.size() }; 
//         }
        
//         void releaseToBookmark(Bookmark bm) {
//             m_currentPageIdx = bm.page;
//             m_currentOffset = bm.offset;
//             m_expressionCache.resize(bm.cacheSize);
//         }

//         /**
//          * @brief Calculates total memory allocated by the Arena in bytes.
//          */
//         size_t getArenaMemoryBytes() const {
//             size_t total = 0;
//             total += sizeof(m_pages);
//             for (const auto& page : m_pages) {
//                 total += sizeof(page); // objet vector
//                 total += page.capacity() * sizeof(CfemReal);
//             }
//             return total;
//         }

//         /**
//          * @brief Debug print of the Arena Memory state to the standard output.
//          */
//         void printArenaMemory() const {
//             size_t total = 0;
//             std::cout << "---- Arena Memory Debug ----\n";
//             std::cout << "Number of pages: " << m_pages.size() << "\n";

//             for (size_t i = 0; i < m_pages.size(); ++i) {
//                 const auto& page = m_pages[i];
//                 size_t mem = sizeof(page) + page.capacity() * sizeof(CfemReal);
//                 std::cout << "Page " << i
//                           << " | size: " << page.size()
//                           << " | capacity: " << page.capacity()
//                           << " | mem: " << mem / 1024.0 << " KB\n";
//                 total += mem;
//             }
//             std::cout << "TOTAL arena memory: " << total / (1024.0 * 1024.0) << " MB\n";
//         }
//     };
// }


#pragma once

#include "cfemutils.h"
#include "Vector3D.h"
#include "ShapeInfo.h"
#include "ShapeInfoBatch.h"
#include "MappingInfo.h"
#include "MappingInfoBatch.h"
#include "MappingInfoBatchView.h" 
#include "DynamicVector.h"
#include "DynamicMatrix.h"

#include <vector>
#include <span>
#include <array>
#include <deque>
#include <iostream>
#include <algorithm>

namespace cfem {

    class FESpace; // Forward declaration

    /**
     * @struct SpaceDataEntry
     * @brief Ultra-fast link between an FE Space and its evaluations at the current integration cell.
     */
    struct SpaceDataEntry {
        const ShapeInfo* shapeSingle   = nullptr; 
        const MappingInfo* mappingSingle = nullptr; 
        
        const ShapeInfoBatch* shapeBatch    = nullptr; 
        const MappingInfoBatchView* mappingBatch  = nullptr; // HPC View
    };

    /**
     * @struct CacheEntry
     * @brief Stores a pointer-based identifier and its associated data span.
     */
    struct CacheEntry {
        const void* id;           
        std::span<CfemReal> data; 
    };

    /**
     * @struct EvalContext
     * @brief High-performance evaluation context for Finite Element assembly.
     */
    struct EvalContext {

        // ====================================================================
        // GLOBAL SIMULATION STATE
        // ====================================================================
        CfemReal currentTime   = 0.0;
        CfemInt  currentCellId = -1;  
        CfemInt currentQpIndex = 0;

        // ====================================================================
        // SINGLE-POINT WORKSPACE (Legacy / Probing)
        // ====================================================================
        struct SingleWorkspace {
            Vector3D    physPos = {};                
            Vector3D    xi      = {-1e9, -1e9, -1e9};
            ShapeInfo   shape;                       
            MappingInfo mapping;                     
        } singleWS;

        // ====================================================================
        // BATCH WORKSPACE (Vectorized Element Assembly)
        // ====================================================================
        struct BatchWorkspace {
            DynamicVector<Vector3D>* physPos = nullptr; 
            DynamicVector<Vector3D>* xi = nullptr;      
            ShapeInfoBatch* shape = nullptr;       
            MappingInfoBatchView* mapping = nullptr; // HPC View
        } batchWS;

        // ====================================================================
        // MULTI-SPACE CACHE & SHORTCUTS (O(1) Access)
        // ====================================================================
        static constexpr CfemInt MAX_CACHED_SPACES = 16;

        #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
            std::vector<SpaceDataEntry> m_spaceDataCache;
        #else
            std::array<SpaceDataEntry, MAX_CACHED_SPACES> m_spaceDataCache;
        #endif

        const ShapeInfoBatch* trialShapeBatch     = nullptr;
        const ShapeInfoBatch* testShapeBatch      = nullptr;
        const MappingInfoBatchView* trialMappingBatch   = nullptr;
        const MappingInfoBatchView* testMappingBatch    = nullptr;
        const MappingInfoBatchView* geometricMappingBatch = nullptr;

        const MappingInfo* geometricMapping = nullptr;

        CfemInt trialDof = -1; 
        CfemInt testDof  = -1; 

        // ====================================================================
        // PRE-ALLOCATED SCRATCHPADS
        // ====================================================================
        DynamicVector<Vector3D> scratchCoords;   
        std::vector<CfemInt>    scratchDofs;     
        DynamicMatrix<CfemReal> scratchCellDofs; 
        MappingInfo             scratchMapping;  
        ShapeInfo               scratchShape;  
        ShapeInfoBatch          scratchShapeBatch;
        MappingInfoBatch        scratchMappingBatch;  
        
        // ====================================================================
        // ARENA MEMORY SYSTEM
        // ====================================================================
        static constexpr size_t PAGE_SIZE = 4096;

        std::deque<std::vector<CfemReal>> m_pages;  
        size_t m_currentPageIdx = 0;
        size_t m_currentOffset   = 0;

        std::vector<CacheEntry> m_expressionCache; 

        // ====================================================================
        // API METHODS
        // ====================================================================

        EvalContext() 
        {
            m_pages.emplace_back(PAGE_SIZE);
            m_expressionCache.reserve(16);
#ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
            m_spaceDataCache.reserve(8);
#endif
        }

#ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
        inline void prepareSpaceCache(CfemInt n) {
            m_spaceDataCache.resize(n);
        }
#endif

        inline const SpaceDataEntry& getSpaceData(CfemInt index) const 
        {
        #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
            CFEM_ASSERT(index >= 0 && index < static_cast<int>(m_spaceDataCache.size()), "EvalContext::getSpaceData index out of bounds");
        #else
            CFEM_ASSERT(index >= 0 && index < MAX_CACHED_SPACES, "EvalContext::getSpaceData index out of bounds");
        #endif
            return m_spaceDataCache[index];
        }

        inline void setSpaceData(CfemInt index, const ShapeInfo* shape, const MappingInfo* mapping) 
        {
            CFEM_ASSERT(index >= 0, "Negative space index");
#ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
            CFEM_ASSERT(index < static_cast<int>(m_spaceDataCache.size()), "Space cache not pre-sized");
#else
            CFEM_ASSERT(index < MAX_CACHED_SPACES, "Space index exceeds MAX_CACHED_SPACES");
#endif
            auto& entry = m_spaceDataCache[index];
            entry.shapeSingle = shape;
            entry.mappingSingle = mapping;
        }

        inline void setSpaceDataBatch(CfemInt index, const ShapeInfoBatch* shape, const MappingInfoBatchView* mapping) 
        {
            CFEM_ASSERT(index >= 0, "Negative space index");
#ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
            CFEM_ASSERT(index < static_cast<int>(m_spaceDataCache.size()), "Space cache not pre-sized");
#else
            CFEM_ASSERT(index < MAX_CACHED_SPACES, "Space index exceeds MAX_CACHED_SPACES");
#endif
            auto& entry = m_spaceDataCache[index];
            entry.shapeBatch = shape;
            entry.mappingBatch = mapping;
        }

        inline const MappingInfoBatchView* getMappingBatch(CfemInt index) const 
        {
        #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
            CFEM_ASSERT(index >= 0 && index < static_cast<int>(m_spaceDataCache.size()), "EvalContext::getMappingBatch index out of bounds");
        #else
            CFEM_ASSERT(index >= 0 && index < MAX_CACHED_SPACES, "EvalContext::getMappingBatch index out of bounds");
        #endif
            const auto* mapping = m_spaceDataCache[index].mappingBatch;
            CFEM_ASSERT(mapping != nullptr, "MappingInfoBatchView not set");
            return mapping;
        }

        inline const MappingInfo* getMapping(CfemInt index) const 
        {
        #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
            CFEM_ASSERT(index >= 0 && index < static_cast<int>(m_spaceDataCache.size()), "EvalContext::getMappingSingle index out of bounds");
        #else
            CFEM_ASSERT(index >= 0 && index < MAX_CACHED_SPACES, "EvalContext::getMappingSingle index out of bounds");
        #endif
            const auto* mapping = m_spaceDataCache[index].mappingSingle;
            CFEM_ASSERT(mapping != nullptr, "MappingInfo not set");
            return mapping;
        }

        inline const ShapeInfoBatch* getShapeBatch(CfemInt index) const 
        {
        #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
            CFEM_ASSERT(index >= 0 && index < static_cast<int>(m_spaceDataCache.size()), "EvalContext::getShapeBatch index out of bounds");
        #else
            CFEM_ASSERT(index >= 0 && index < MAX_CACHED_SPACES, "EvalContext::getShapeBatch index out of bounds");
        #endif
            const auto* shape = m_spaceDataCache[index].shapeBatch;
            CFEM_ASSERT(shape != nullptr, "ShapeInfoBatch not set. Did you call setSpaceDataBatch?");
            return shape;
        }

        inline const ShapeInfo* getShape(CfemInt index) const 
        {
        #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
            CFEM_ASSERT(index >= 0 && index < static_cast<int>(m_spaceDataCache.size()), "EvalContext::getShapeSingle index out of bounds");
        #else
            CFEM_ASSERT(index >= 0 && index < MAX_CACHED_SPACES, "EvalContext::getShapeSingle index out of bounds");
        #endif
            const auto* shape = m_spaceDataCache[index].shapeSingle;
            CFEM_ASSERT(shape != nullptr, "ShapeInfo not set. Did you call setSpaceDataSingle?");
            return shape;
        }

        void resetBuffers() {
            m_currentPageIdx = 0;
            m_currentOffset = 0;
            m_expressionCache.clear();
        }

        std::span<CfemReal> getBuffer(CfemInt size) {
            CFEM_ASSERT(size >= 0, "Negative buffer size");
            CFEM_ASSERT(m_currentPageIdx < m_pages.size(),
                        "Invalid arena state: page index out of bounds");

            size_t usize = static_cast<size_t>(size);

            if (m_currentOffset + usize > m_pages[m_currentPageIdx].size()){
                m_currentPageIdx++;
                m_currentOffset = 0;
                
                if (m_currentPageIdx >= m_pages.size()) {
                    m_pages.emplace_back(std::max(usize, PAGE_SIZE));
                } 
                else if (m_pages[m_currentPageIdx].size() < usize) {
                    m_pages[m_currentPageIdx].resize(usize);
                }
            }
            
            auto span = std::span<CfemReal>(m_pages[m_currentPageIdx].data() + m_currentOffset, usize);
            m_currentOffset += usize;
            return span;
        }

        bool tryGetCache(const void* exprPtr, std::span<CfemReal> out) const {
            for (const auto& entry : m_expressionCache) {
                if (entry.id == exprPtr) {
                    CFEM_ASSERT(out.size() >= entry.data.size(), "Output buffer too small in tryGetCache");
                    std::copy(entry.data.begin(), entry.data.end(), out.begin());
                    return true;
                }
            }
            return false;
        }

        void storeInCache(const void* exprPtr, std::span<CfemReal> result) {
            auto storage = getBuffer(static_cast<CfemInt>(result.size()));
            std::copy(result.begin(), result.end(), storage.begin());
            m_expressionCache.push_back({exprPtr, storage});
        }

        struct Bookmark { size_t page; size_t offset; size_t cacheSize; };
        
        Bookmark getBookmark() const {
            return { m_currentPageIdx, m_currentOffset, m_expressionCache.size() }; 
        }
        
        void releaseToBookmark(Bookmark bm) {
            m_currentPageIdx = bm.page;
            m_currentOffset = bm.offset;
            m_expressionCache.resize(bm.cacheSize);
        }

        size_t getArenaMemoryBytes() const {
            size_t total = 0;
            total += sizeof(m_pages);
            for (const auto& page : m_pages) {
                total += sizeof(page); 
                total += page.capacity() * sizeof(CfemReal);
            }
            return total;
        }

        void printArenaMemory() const {
            size_t total = 0;
            std::cout << "---- Arena Memory Debug ----\n";
            std::cout << "Number of pages: " << m_pages.size() << "\n";

            for (size_t i = 0; i < m_pages.size(); ++i) {
                const auto& page = m_pages[i];
                size_t mem = sizeof(page) + page.capacity() * sizeof(CfemReal);
                std::cout << "Page " << i
                          << " | size: " << page.size()
                          << " | capacity: " << page.capacity()
                          << " | mem: " << mem / 1024.0 << " KB\n";
                total += mem;
            }
            std::cout << "TOTAL arena memory: " << total / (1024.0 * 1024.0) << " MB\n";
        }
    };
}