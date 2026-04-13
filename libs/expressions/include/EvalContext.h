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

// #pragma once

// #include "cfemutils.h"
// #include "Vector3D.h"
// #include "ReferenceElement.h"
// #include "DynamicVector.h"

// #include <deque>
// #include <vector>
// #include <span>
// #include <array>
// #include <memory>

// namespace cfem {

//     class FESpace; // Forward declaration

//     /**
//      * @brief Entry for the linear expression cache. 
//      * id: The memory address of the evaluated expression.
//      * data: A span pointing to the results stored in the EvalContext pages.
//      */
//     struct CacheEntry {
//         const void* id;
//         std::span<CfemReal> data;
//     };

//     /**
//      * @brief High-performance entry linking a finite element space to its 
//      * evaluated shape functions and physical mappings at the current integration point.
//      */
//     struct SpaceDataEntry {
//         const FESpace* space = nullptr;
//         const ShapeInfo* shape = nullptr;
//         const MappingInfo* mapping = nullptr;
//     };

//     /**
//      * @brief The Evaluation Context passed around during matrix assembly and field projection.
//      * It acts as a thread-local memory arena and a cache for geometric/algebraic evaluations
//      * to prevent heap allocations during deep computational loops.
//      */
//     struct EvalContext {
//         // --- Integration Point Context ---
//         Vector3D physPos = {};                    ///< Physical coordinates of the current evaluation point
//         CfemReal currentTime = 0;                 ///< Current simulation time (for transient problems)
//         CfemInt currentCellId = -1;               ///< ID of the mesh cell currently being processed
//         Vector3D currentXi = {-1e9, -1e9, -1e9};  ///< Reference coordinates (xi, eta, zeta) of the point

//         // --- Multi-Space Geometry Cache ---
//         static constexpr int MAX_CACHED_SPACES = 8; ///< Maximum number of coupled spaces evaluated simultaneously
//         std::array<SpaceDataEntry, MAX_CACHED_SPACES> m_spaceDataCache;
//         int m_numCachedSpaces = 0;

//         // --- Assembler Shortcuts (Bilinear/Linear Forms) ---
//         const ShapeInfo* trialShape = nullptr;
//         const ShapeInfo* testShape  = nullptr;
//         const MappingInfo* trialMapping = nullptr;
//         const MappingInfo* testMapping  = nullptr;
//         CfemInt trialDof = -1;
//         CfemInt testDof  = -1;

//         // --- Scratchpads (Pre-allocated buffers to avoid mallocs) ---
//         ShapeInfo scratchShape;
//         MappingInfo scratchMapping;
//         DynamicVector<Vector3D> scratchCoords;
//         std::vector<CfemInt> scratchDofs;

//     private:
//         // --- Memory Arena (Pages) ---
//         static constexpr size_t PAGE_SIZE = 4096; 
//         std::deque<std::vector<CfemReal>> m_pages; 
//         size_t m_currentPageIdx = 0;
//         size_t m_currentOffset = 0;

//         // --- Expression Cache ---
//         // HIGH PERFORMANCE CACHE: Linear search is faster than a Map for small N.
//         std::vector<CacheEntry> m_expressionCache;

//     public:
//         EvalContext() {
//             m_pages.emplace_back(PAGE_SIZE);
//             m_expressionCache.reserve(16); // Pre-reserve to avoid reallocations
//         }

//         // =====================================================================
//         // FESpace Geometry Management (The New Abstraction)
//         // =====================================================================

//         /**
//          * @brief Registers the evaluated shape and mapping for a specific FESpace.
//          * Typically called by the Assembler before evaluating user expressions.
//          */
//         void registerSpaceData(const std::shared_ptr<const FESpace>& space, const ShapeInfo* shape, const MappingInfo* mapping = nullptr) {
//             if (m_numCachedSpaces < MAX_CACHED_SPACES) {
//                 m_spaceDataCache[m_numCachedSpaces++] = {space.get(), shape, mapping};
//             } else {
//                 CFEM_ERROR << "EvalContext: Exceeded maximum number of coupled spaces (" << MAX_CACHED_SPACES << ").";
//             }
//         }

//         /**
//          * @brief Checks if shape functions for the requested space have been computed and cached.
//          */
//         bool hasValidShape(const std::shared_ptr<const FESpace>& space) const {
//             for (int i = 0; i < m_numCachedSpaces; ++i) {
//                 if (m_spaceDataCache[i].space == space.get() && m_spaceDataCache[i].shape != nullptr) {
//                     return true;
//                 }
//             }
//             return false;
//         }

//         /**
//          * @brief Retrieves the cached shape functions for the given space.
//          * @warning Call `hasValidShape` first, or ensure the assembler populated the context.
//          */
//         const ShapeInfo& getShape(const std::shared_ptr<const FESpace>& space) const {
//             for (int i = 0; i < m_numCachedSpaces; ++i) {
//                 if (m_spaceDataCache[i].space == space.get()) {
//                     return *(m_spaceDataCache[i].shape);
//                 }
//             }
//             throw std::runtime_error("EvalContext: Requested ShapeInfo not found for the given FESpace.");
//         }

//         /**
//          * @brief Retrieves the cached mapping (Jacobians, physical gradients) for the given space.
//          */
//         const MappingInfo& getMapping(const std::shared_ptr<const FESpace>& space) const {
//             for (int i = 0; i < m_numCachedSpaces; ++i) {
//                 if (m_spaceDataCache[i].space == space.get() && m_spaceDataCache[i].mapping != nullptr) {
//                     return *(m_spaceDataCache[i].mapping);
//                 }
//             }
//             throw std::runtime_error("EvalContext: Requested MappingInfo not found for the given FESpace.");
//         }

//         inline const SpaceDataEntry& getSpaceData(int index) const {
//             return m_spaceDataCache[index];
//         }

//         // On garde registerSpaceData mais on l'utilise différemment
//         void setSpaceData(int index, const FESpace* space, const ShapeInfo* shape, const MappingInfo* mapping) {
//             m_spaceDataCache[index] = {space, shape, mapping};
//             if (index >= m_numCachedSpaces) m_numCachedSpaces = index + 1;
//         }

//         // =====================================================================
//         // Memory Arena & Expression Cache Methods (Untouched Core Logic)
//         // =====================================================================

//         /**
//          * @brief Requests a contiguous block of memory from the arena.
//          * If the current page is full, it moves to the next one.
//          */
//         std::span<CfemReal> getBuffer(CfemInt size) {
//             if (size > PAGE_SIZE) {
//                 m_pages.insert(m_pages.begin() + m_currentPageIdx, std::vector<CfemReal>(size));
//                 auto& specialPage = m_pages[m_currentPageIdx];
//                 m_currentPageIdx++; 
//                 return std::span<CfemReal>(specialPage.data(), size);
//             }

//             if (m_currentOffset + size > PAGE_SIZE) {
//                 m_currentPageIdx++;
//                 m_currentOffset = 0;
//                 if (m_currentPageIdx >= m_pages.size()) {
//                     m_pages.emplace_back(PAGE_SIZE);
//                 }
//             }

//             auto span = std::span<CfemReal>(m_pages[m_currentPageIdx].data() + m_currentOffset, size);
//             m_currentOffset += size;
//             return span;
//         }

//         /**
//          * @brief Checks if the mathematical expression has already been evaluated at this point.
//          * Useful to avoid re-evaluating shared sub-expressions (e.g., repeating `sin(x)`).
//          */
//         bool tryGetCache(const void* exprPtr, std::span<CfemReal> out) const {
//             for (const auto& entry : m_expressionCache) {
//                 if (entry.id == exprPtr) {
//                     std::copy(entry.data.begin(), entry.data.end(), out.begin());
//                     return true;
//                 }
//             }
//             return false;
//         }

//         /**
//          * @brief Stores the evaluation result in the fast page buffer and indexes it.
//          */
//         void storeInCache(const void* exprPtr, std::span<CfemReal> result) {
//             auto cacheStorage = getBuffer(static_cast<CfemInt>(result.size()));
//             std::copy(result.begin(), result.end(), cacheStorage.begin());
//             m_expressionCache.push_back({exprPtr, cacheStorage});
//         }

//         /**
//          * @brief Resets the memory arena and caches. 
//          * Must be called at the beginning of every new integration point (Quadrature point).
//          */
//         void resetBuffers() {
//             m_currentPageIdx = 0;
//             m_currentOffset = 0;
//             m_expressionCache.clear(); // Resets size to 0 but keeps capacity
//             m_numCachedSpaces = 0;     // Flushes the geometric space cache
//         }

//         /**
//          * @brief Bookmark system to free memory within the same integration point
//          * when temporary buffers are no longer needed.
//          */
//         struct Bookmark { size_t page; size_t offset; };
        
//         Bookmark getBookmark() const { return { m_currentPageIdx, m_currentOffset }; }
        
//         void releaseToBookmark(Bookmark bm) {
//             m_currentPageIdx = bm.page;
//             m_currentOffset = bm.offset;
//         }
//     };
// } // namespace cfem


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
     * @brief Lien ultra-rapide entre un espace et ses évaluations au point courant.
     * On utilise des pointeurs bruts pour éviter le surcoût des shared_ptr.
     */
    struct SpaceDataEntry {
        const ShapeInfo* shape = nullptr;
        const MappingInfo* mapping = nullptr;
    };

    struct CacheEntry {
        const void* id;
        std::span<CfemReal> data;
    };

    /**
     * @brief Version haute-performance de EvalContext.
     * Les accès aux données géométriques sont désormais en O(1) via indexation directe.
     */
    struct EvalContext {
        // --- Contexte du point d'intégration ---
        Vector3D physPos = {}; 
        CfemReal currentTime = 0;
        CfemInt currentCellId = -1;
        Vector3D currentXi = {-1e9, -1e9, -1e9};

        // --- Cache Multi-Espaces (Accès O(1)) ---
        static constexpr int MAX_CACHED_SPACES = 8;
        std::array<SpaceDataEntry, MAX_CACHED_SPACES> m_spaceDataCache;

        // --- Raccourcis pour les Formes Variationnelles (Trial/Test) ---
        // Optionnel : peut être bindé comme les autres espaces, mais utile pour la clarté
        const ShapeInfo* trialShape = nullptr;
        const ShapeInfo* testShape   = nullptr;
        const MappingInfo* trialMapping = nullptr;
        const MappingInfo* testMapping  = nullptr;

        CfemInt trialDof = -1;
        CfemInt testDof  = -1;

        // --- Scratchpads ---
        ShapeInfo scratchShape;
        MappingInfo scratchMapping;
        DynamicVector<Vector3D> scratchCoords;
        std::vector<CfemInt> scratchDofs;

        // --- Arena Memory ---
        static constexpr size_t PAGE_SIZE = 4096;
        std::deque<std::vector<CfemReal>> m_pages;
        size_t m_currentPageIdx = 0;
        size_t m_currentOffset = 0;
        std::vector<CacheEntry> m_expressionCache;

        EvalContext() {
            m_pages.emplace_back(PAGE_SIZE);
            m_expressionCache.reserve(16);
        }

        /**
         * @brief Accès direct à une donnée d'espace par index (Performance maximale).
         * @param index L'ID attribué par l'intégrateur AVANT la boucle.
         */
        inline const SpaceDataEntry& getSpaceData(int index) const {
            return m_spaceDataCache[index];
        }

        /**
         * @brief Enregistre les données pour un index spécifique.
         * Utilisé dans la boucle de quadrature.
         */
        inline void setSpaceData(int index, const ShapeInfo* shape, const MappingInfo* mapping) {
            m_spaceDataCache[index].shape = shape;
            m_spaceDataCache[index].mapping = mapping;
        }

        /**
         * @brief Récupère le mapping pré-calculé pour un index donné.
         * @return pointeur vers MappingInfo ou nullptr si non défini pour cet index.
         */
        inline const MappingInfo* getMapping(int index) const {
            if (index < 0 || index >= MAX_CACHED_SPACES) return nullptr;
            return m_spaceDataCache[index].mapping;
        }

        // --- Gestion de la Mémoire Arena ---

        void resetBuffers() {
            m_currentPageIdx = 0;
            m_currentOffset = 0;
            m_expressionCache.clear();
            // On ne reset pas m_spaceDataCache, il est écrasé à chaque point.
        }

        std::span<CfemReal> getBuffer(CfemInt size) {
            if (size > (CfemInt)PAGE_SIZE) {
                m_pages.insert(m_pages.begin() + m_currentPageIdx, std::vector<CfemReal>(size));
                return std::span<CfemReal>(m_pages[m_currentPageIdx++].data(), size);
            }
            if (m_currentOffset + size > PAGE_SIZE) {
                m_currentPageIdx++;
                m_currentOffset = 0;
                if (m_currentPageIdx >= m_pages.size()) m_pages.emplace_back(PAGE_SIZE);
            }
            auto span = std::span<CfemReal>(m_pages[m_currentPageIdx].data() + m_currentOffset, size);
            m_currentOffset += size;
            return span;
        }

        // --- Cache d'expressions (CfemExpression) ---

        bool tryGetCache(const void* exprPtr, std::span<CfemReal> out) const {
            for (const auto& entry : m_expressionCache) {
                if (entry.id == exprPtr) {
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

        /**
         * @brief Bookmark system to free memory within the same integration point
         * when temporary buffers are no longer needed.
         */
        struct Bookmark { size_t page; size_t offset; };
        
        Bookmark getBookmark() const {
            return { m_currentPageIdx, m_currentOffset }; 
        }
        
        void releaseToBookmark(Bookmark bm) {
            m_currentPageIdx = bm.page;
            m_currentOffset = bm.offset;
        }
    };
}