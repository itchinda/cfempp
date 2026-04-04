#ifndef CFEM_TOPOLOGY_H
#define CFEM_TOPOLOGY_H

#include "cfemutils.h"
#include "Vector3D.h"
#include <algorithm>
#include <cstddef>
#include <iostream>
#include <span>

namespace cfem{

enum class MeshOrder{ None, Linear = 1, Quadratic = 2};

struct Vertex
{
    cfem::Vector3D coords;
    //std::vector<CfemInt> connected_cells;
    int tag_list_id = -1;

    const CfemInt* connected_cells_ptr = nullptr;
    CfemInt num_connected_cells = 0;

    // Helper for easier looping
    auto getConnectedCells() const {
        return std::span<const CfemInt>(connected_cells_ptr, num_connected_cells);
    }
};

enum class CellType
{
    LINE2, TRI3,QUAD4,TET4,PRISM6,HEX8,
    LINE3, TRI6, QUAD9, TET10, PRISM18, HEX27,
    NODE,
    COUNT
};

// On utilise un tableau indexé par l'entier de l'enum CellType
static const std::vector<CellType> cellFacetMapping[static_cast<int>(CellType::COUNT)] = {
    /* [CellType::LINE2]        */ { CellType::NODE},
    /* [CellType::TRI3]         */ { CellType::LINE2 },
    /* [CellType::QUAD4]        */ { CellType::LINE2 },
    /* [CellType::TET4]         */ { CellType::TRI3 },
    /* [CellType::PRISM6]       */ { CellType::TRI3, CellType::QUAD4 },
    /* [CellType::HEX8]         */ { CellType::QUAD4 },
    /* [CellType::LINE3]        */ { CellType::NODE },
    /* [CellType::TRI6]         */ { CellType::LINE3 },
    /* [CellType::QUAD9]        */ { CellType::LINE3 },
    /* [CellType::TET10]        */ { CellType::TRI6 },
    /* [CellType::PRISM18]      */ { CellType::TRI6, CellType::QUAD9 },
    /* [CellType::HEX27]        */ { CellType::QUAD9 },
    /* [CellType::NODE]         */ { }
};

struct Cell {
    CfemInt id{-1};
    CellType type;
    int tag_list_id{-1};

    size_t offset{0};
    size_t node_count{0};

    // --- Consistent with Vertex: View-based pointers ---
    const CfemInt* neighbor_cells_ptr = nullptr;
    const CfemInt* facet_ids_ptr = nullptr;
    CfemInt num_faces = 0;

    /** 
     * @brief Helpers for easier looping
     * @code
     *      for (CfemInt neighborId : cell.getNeighborCells()) {
     *          if (neighborId != -1) {
     *              // Process internal neighbor
     *          } else {
     *              // This face is on the boundary
     *          }
     *      }
     */
    auto getNeighborCells() const {
        return std::span<const CfemInt>(neighbor_cells_ptr, num_faces);
    }

    auto getCellFacets() const {
        return std::span<const CfemInt>(facet_ids_ptr, num_faces);
    }

    inline static CellType getLinearEquivalent(CellType type) {
        switch (type) {
            case CellType::LINE3:  return CellType::LINE2;
            case CellType::TRI6:  return CellType::TRI3;
            case CellType::QUAD9: return CellType::QUAD4;
            case CellType::TET10: return CellType::TET4;
            case CellType::PRISM18: return CellType::PRISM6;
            case CellType::HEX27: return CellType::HEX8;
            default: return type;
        }
    }

    inline static CellType getQuadraticEquivalent(CellType type) {
        switch (type) {
            case CellType::LINE2:  return CellType::LINE3;
            case CellType::TRI3:  return CellType::TRI6;
            case CellType::QUAD4: return CellType::QUAD9;
            case CellType::TET4: return CellType::TET10;
            case CellType::PRISM6: return CellType::PRISM18;
            case CellType::HEX8: return CellType::HEX27;
            default: return type;
        }
    }

    inline static std::span<const CellType> getFacetTypes(CellType type) {
        int idx = static_cast<int>(type);
        if (idx < 0 || idx >= static_cast<int>(CellType::COUNT)) {
            return {};
        }
        return cellFacetMapping[idx];
    }
};

struct Facet {
    CfemInt id{-1};
    CellType type;
    CfemInt left_cell{-1};
    CfemInt right_cell{-1};
    
    // Index-based storage
    size_t offset{0};
    size_t node_count{0};

    bool isBoundary() const { return right_cell == -1; }
};


} // namespace
#endif