#ifndef CFEM_TOPOLOGY_H
#define CFEM_TOPOLOGY_H

#include "cfemutils.h"
#include "Vector3D.h"

namespace cfem{

struct Vertex
{
    cfem::Vector3D coords;
    std::vector<CfemInt> connected_cells;
    int tag_list_id = -1;
};

enum class CellType
{
    LINE2, TRI3,QUAD4,TET4,PRISM6,HEX8,
    LINE3, TRI6, QUAD9, TET10, PRISM18, HEX27,
    NODE
};

struct Cell
{
    CellType type;
    std::vector<CfemInt> vertex_ids{};
    int tag_list_id = -1;

    inline CellType getLinearEquivalent(CellType type) {
    switch (type) {
        case CellType::LINE3:  return CellType::LINE2;
        case CellType::TRI6:  return CellType::TRI3;
        case CellType::QUAD9: return CellType::QUAD4;
        case CellType::TET10: return CellType::TET4;
        case CellType::PRISM6: return CellType::PRISM18;
        case CellType::HEX27: return CellType::HEX8;
        default: return type;
    }
}
};

struct Facet
{
    CfemInt id{-1};
    CellType type;
    std::vector<CfemInt> vertex_ids{};
    CfemInt number_of_vertices{0};
    CfemInt left_cell{-1};
    CfemInt right_cell{-1};
    int tag_list_id = -1;

    bool isBoundary() const { return right_cell < 0; }
};

} // namespace
#endif