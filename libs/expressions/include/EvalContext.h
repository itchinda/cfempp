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
#include "Vector3D.h"
#include "ReferenceElement.h"
#include "DynamicVector.h"

#include <vector>
#include <span>
#include <array>
#include <deque>

namespace cfem {

    class FESpace; // Forward declaration

    /**
     * @struct SpaceDataEntry
     * @brief Ultra-fast link between an FE Space and its evaluations at the current integration point.
     * * This structure uses raw pointers to avoid the overhead of smart pointers during 
     * heavy quadrature loops. It groups shape functions and geometric mapping data.
     */
    struct SpaceDataEntry {
        const ShapeInfo* shape = nullptr;     /**< Pointer to precomputed shape function values/gradients. */
        const MappingInfo* mapping = nullptr; /**< Pointer to geometric mapping info (Jacobians, etc.). */
    };

    /**
     * @struct CacheEntry
     * @brief Stores a pointer-based identifier and its associated data span.
     * Used for caching symbolic expression results to avoid redundant computations.
     */
    struct CacheEntry {
        const void* id;           /**< Unique identifier (usually the expression's 'this' pointer). */
        std::span<CfemReal> data; /**< Cached numerical values. */
    };

    /**
     * @struct EvalContext
     * @brief High-performance evaluation context for Finite Element assembly.
     * This structure centralizes all data required at a specific integration point (quadrature point).
     * It provides:
     * 1. O(1) access to geometric and shape data via direct indexing.
     * 2. An Arena Memory system to avoid dynamic allocations during assembly.
     * 3. A symbolic expression cache to optimize complex variational forms.
     */
    struct EvalContext {

        // --- Integration Point Context ---
        Vector3D physPos = {};                   /**< Physical coordinates of the point. */
        CfemReal currentTime = 0;                /**< Current simulation time (for transient problems). */
        CfemInt  currentCellId = -1;             /**< Index of the cell currently being integrated. */
        Vector3D currentXi = {-1e9, -1e9, -1e9}; /**< Local/Reference coordinates. */

        // --- Multi-Space Cache (O(1) Access) ---
        static constexpr int MAX_CACHED_SPACES  = 16; /**< Maximum number of FE spaces allowed in a single variational term (e.g., for multi-physics). */

        /**
         * @name Space Data Cache
         * @brief Storage for FE space evaluations.
         * @details Defaults to a fixed-size array for O(1) stack-like speed.
         * Define CFEM_USE_DYNAMIC_SPACE_CACHE to use std::vector for unlimited spaces.
         * @{
         */
        #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
            std::vector<SpaceDataEntry> m_spaceDataCache;
        #else
            std::array<SpaceDataEntry, MAX_CACHED_SPACES > m_spaceDataCache;
        #endif
        /** @} */

        // --- Variational Form Shortcuts (Trial/Test) ---
        /** 
         * @name Variational Form Shortcuts
         * @brief Shortcuts for standard u (Trial) and v (Test) functions to simplify kernel code.
         * @{ 
         */
        const ShapeInfo*   trialShape   = nullptr;
        const ShapeInfo*   testShape    = nullptr;
        const MappingInfo* trialMapping = nullptr;
        const MappingInfo* testMapping  = nullptr;
        /** @} */

        CfemInt trialDof = -1; /**< Current local Trial DOF index being processed. */
        CfemInt testDof  = -1; /**< Current local Test DOF index being processed. */

        /** 
         * @name Pre-allocated Scratchpads
         * @brief Buffers used for temporary local computations to avoid heap allocations 
         * inside the element-level loops.
         * @{ 
         */
        /** @brief Buffer for local basis function values and gradients. */
        ShapeInfo scratchShape;

        /** @brief Buffer for local geometric transformation data (Jacobians, determinants). */
        MappingInfo scratchMapping;

        /** @brief Pre-allocated storage for cell vertex coordinates. */
        DynamicVector<Vector3D> scratchCoords;

        /** @brief Reusable buffer for gathering local DOF indices from the global system. */
        std::vector<CfemInt> scratchDofs;
        /** @} */

        static constexpr size_t PAGE_SIZE = 4096;  /**< Size of a single memory page in doubles. */

        /** 
         * @name Arena Memory System
         * @brief Memory management to avoid heap allocations during quadrature loops.
         * @{
         */
        /** @brief Deque of pages; ensures pointers remain valid when new pages are allocated. */
        std::deque<std::vector<CfemReal>> m_pages;  
        
        /** @brief Index of the page currently being used for allocation. */
        size_t m_currentPageIdx = 0;
        
        /** @brief Current position (in CfemReal units) within the active page. */
        size_t m_currentOffset = 0;
        /** @} */

        std::vector<CacheEntry> m_expressionCache; /** Storage for expression caching results. */

        /**
         * @brief Initializes the EvalContext with a primary memory page.
         */
        EvalContext() 
        {
            m_pages.emplace_back(PAGE_SIZE);
            m_expressionCache.reserve(16);
#ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
            m_spaceDataCache.reserve(8);
#endif
        }

#ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
        /**
         * @brief Prepares the internal space cache for use.
         *
         * Must be called before any call to setSpaceData() when using
         * dynamic space cache mode.
         */
        inline void prepareSpaceCache(CfemInt n) {
            m_spaceDataCache.resize(n);
        }
#endif

        /**
         * @brief Direct access to space data by index.
         * @param index The ID assigned by the Integrator before the assembly loop.
         * @return A reference to the SpaceDataEntry (O(1)).
         */
        inline const SpaceDataEntry& getSpaceData(int index) const 
        {
        #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
            CFEM_ASSERT(index >= 0 && index < static_cast<int>(m_spaceDataCache.size()),
                "EvalContext::getSpaceData index out of bounds (dynamic cache)");
        #else
            CFEM_ASSERT(index >= 0 && index < MAX_CACHED_SPACES,
                "EvalContext::getSpaceData index out of bounds (fixed cache)");
        #endif

            return m_spaceDataCache[index];
        }

        /**
         * @brief Binds FE space data (shape/mapping) to a specific index.
         * @param index Slot index assigned to the FE Space.
         * @param shape Pointer to shape function evaluations at the current point.
         * @param mapping Pointer to geometric mapping info at the current point.
         * 
         * @note Typically called before evaluating terms at a new integration point.
         * @warning Pointers must remain valid for the duration of the evaluation.
         */
        inline void setSpaceData(int index, const ShapeInfo* shape, const MappingInfo* mapping) 
        {
            CFEM_ASSERT(index >= 0, "Negative space index");

#ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
            CFEM_ASSERT(index < static_cast<int>(m_spaceDataCache.size()),
                "Space cache not pre-sized");
#else
            CFEM_ASSERT(index < MAX_CACHED_SPACES,
                "Space index exceeds MAX_CACHED_SPACES");
#endif

            auto& entry = m_spaceDataCache[index];
            entry.shape = shape;
            entry.mapping = mapping;
        }

        /**
         * @brief Retrieves precomputed mapping info for a given space index.
         * @param index Space identifier.
         * @return Pointer to MappingInfo.
         * @pre index is valid and mapping has been set.
         */
        inline const MappingInfo* getMapping(int index) const 
        {
        #ifdef CFEM_USE_DYNAMIC_SPACE_CACHE
            CFEM_ASSERT(index >= 0 && index < static_cast<int>(m_spaceDataCache.size()),
                "EvalContext::getMapping index out of bounds");
        #else
            CFEM_ASSERT(index >= 0 && index < MAX_CACHED_SPACES,
                "EvalContext::getMapping index out of bounds");
        #endif

            const auto* mapping = m_spaceDataCache[index].mapping;
            CFEM_ASSERT(mapping != nullptr, "MappingInfo not set");

            return mapping;
        }

        // --- Arena Memory Management ---

        /**
         * @brief Resets all temporary buffers and caches.
         * 
         * Called before moving to a new integration point or cell.
         */
        void resetBuffers() {
            m_currentPageIdx = 0;
            m_currentOffset = 0;
            m_expressionCache.clear();
        }

        /**
         * @brief Allocates a buffer from the Arena memory.
         * @param size Number of CfemReal elements needed.
         * @return A span pointing to the allocated memory.
         */
        std::span<CfemReal> getBuffer(CfemInt size) {
            CFEM_ASSERT(size >= 0, "Negative buffer size");
            CFEM_ASSERT(m_currentPageIdx < m_pages.size(),
                        "Invalid arena state: page index out of bounds");

            size_t usize = static_cast<size_t>(size);

            if (m_currentOffset + usize > m_pages[m_currentPageIdx].size()){
                // Check if the current page has enough remaining space
                m_currentPageIdx++;
                m_currentOffset = 0;
                
                // If we ran out of pages, add a new one
                if (m_currentPageIdx >= m_pages.size()) {
                    m_pages.emplace_back(std::max(usize, PAGE_SIZE));
                } 
                // If a page exists from a previous run but is too small, resize it
                else if (m_pages[m_currentPageIdx].size() < usize) {
                    m_pages[m_currentPageIdx].resize(usize);
                }
            }
            
            auto span = std::span<CfemReal>(m_pages[m_currentPageIdx].data() + m_currentOffset, usize);
            m_currentOffset += usize;
            return span;
        }


        // --- Expression Caching (CfemExpression) ---

        /**
         * @brief Attempts to retrieve a result from the expression cache.
         * @param exprPtr Pointer to the expression (id).
         * @param out Span to receive the cached data.
         * @return True if found and copied, false otherwise.
         */
        bool tryGetCache(const void* exprPtr, std::span<CfemReal> out) const {
            for (const auto& entry : m_expressionCache) {
                if (entry.id == exprPtr) {
                    CFEM_ASSERT(out.size() >= entry.data.size(),
                        "Output buffer too small in tryGetCache");
                    std::copy(entry.data.begin(), entry.data.end(), out.begin());
                    return true;
                }
            }
            return false;
        }

        /**
         * @brief Stores an expression result in the cache for later reuse.
         * @param exprPtr Expression identifier.
         * @param result Data to be cached.
         */
        void storeInCache(const void* exprPtr, std::span<CfemReal> result) {
            auto storage = getBuffer(static_cast<CfemInt>(result.size()));
            std::copy(result.begin(), result.end(), storage.begin());
            m_expressionCache.push_back({exprPtr, storage});
        }

        /**
         * @struct Bookmark
         * @brief Represents a state in the Arena Memory to allow partial rollbacks.
         */
        struct Bookmark { size_t page; size_t offset; size_t cacheSize; };
        
        /**
         * @brief Captures the current state of the memory Arena.
         */
        Bookmark getBookmark() const {
            return { m_currentPageIdx, m_currentOffset, m_expressionCache.size() }; 
        }
        
        /**
         * @brief Releases memory and expression cache back to a previously saved bookmark.
         * 
         * This allows partial deallocation within a single integration point. It restores 
         * the Arena to the saved offset and automatically removes any expression cache 
         * entries created after the bookmark was taken to prevent data corruption.
         * 
         * @param bm The bookmark captured via getBookmark().
         */
        void releaseToBookmark(Bookmark bm) {
            m_currentPageIdx = bm.page;
            m_currentOffset = bm.offset;
            m_expressionCache.resize(bm.cacheSize);
        }

        /**
         * @brief Calculates total memory allocated by the Arena in bytes.
         */
        size_t getArenaMemoryBytes() const {
            size_t total = 0;

            total += sizeof(m_pages);

            for (const auto& page : m_pages) {
                total += sizeof(page); // objet vector
                total += page.capacity() * sizeof(CfemReal);
            }

            return total;
        }

        /**
         * @brief Debug print of the Arena Memory state to the standard output.
         */
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

            std::cout << "TOTAL arena memory: "
                    << total / (1024.0 * 1024.0)
                    << " MB\n";
        }
    };
}