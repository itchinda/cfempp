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

        size_t getArenaMemoryBytes() const {
            size_t total = 0;

            // Taille du deque lui-même (approx)
            total += sizeof(m_pages);

            // Taille de chaque page
            for (const auto& page : m_pages) {
                total += sizeof(page); // objet vector
                total += page.capacity() * sizeof(CfemReal); // buffer réel
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

            std::cout << "TOTAL arena memory: "
                    << total / (1024.0 * 1024.0)
                    << " MB\n";
        }
    };
}