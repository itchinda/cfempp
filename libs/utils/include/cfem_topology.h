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
    std::vector<CfemInt> connected_cells;
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


struct Cell {
    CfemInt id{-1};
    CellType type;
    int tag_list_id{-1};

    size_t offset{0};
    size_t node_count{0};

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