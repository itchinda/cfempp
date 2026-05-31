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

#ifndef CFEM_TOPOLOGY_H
#define CFEM_TOPOLOGY_H

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <span>

#include "cfemutils.h"
#include "Vector3D.h"

namespace cfem{

enum class MeshOrder{ None, Linear = 1, Quadratic = 2};

struct Vertex
{
    cfem::Vector3D coords;
    CfemInt tag_list_id = -1;

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
        CfemInt tag_list_id{-1};

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

        inline static CfemInt getNodeCount(CellType type)
        {
            switch (type)
            {
                case CellType::LINE2:   return 2;
                case CellType::TRI3:    return 3;
                case CellType::QUAD4:   return 4;
                case CellType::TET4:    return 4;
                case CellType::PRISM6:  return 6;
                case CellType::HEX8:    return 8;
                case CellType::LINE3:   return 3;
                case CellType::TRI6:    return 6;
                case CellType::QUAD9:   return 9;
                case CellType::TET10:   return 10;
                case CellType::PRISM18: return 18;
                case CellType::HEX27:   return 27;
                default:
                    return 0;
            }
        }

        inline static CfemInt getEdgeCount(CellType type)
        {
            switch (type)
            {
                case CellType::LINE2: 
                case CellType::LINE3:    return 1;
                case CellType::TRI3:
                case CellType::TRI6:     return 3;
                case CellType::QUAD4:
                case CellType::QUAD9:    return 4;
                case CellType::TET4:
                case CellType::TET10:    return 6;
                case CellType::PRISM6:
                case CellType::PRISM18:  return 9;
                case CellType::HEX8:
                case CellType::HEX27:    return 12;
                default:
                    CFEM_ERROR << "Unsupported Cell Type";
                    return 0;
            }
        }

        inline static CfemInt getFaceCount(CellType type)
        {
            switch (type)
            {
                case CellType::LINE2:
                case CellType::LINE3:    return 2;
                case CellType::TRI3:
                case CellType::TRI6:     return 3;
                case CellType::QUAD4:
                case CellType::QUAD9:    return 4;
                case CellType::TET4:
                case CellType::TET10:    return 4;
                case CellType::PRISM6:
                case CellType::PRISM18:  return 5;
                case CellType::HEX8:
                case CellType::HEX27:    return 6;
                default:
                    CFEM_ERROR << "Unsupported Cell Type";
                    return 0;
            }
        }

        inline static std::span<const CellType> getFacetTypes(CellType type) {
            CfemInt idx = static_cast<CfemInt>(type);
            if (idx < 0 || idx >= static_cast<CfemInt>(CellType::COUNT)) {
                return {};
            }
            return cellFacetMapping[idx];
        }

        inline static CellType getFacetType(CellType type, CfemInt localFacetId) {
            auto facetTypes = getFacetTypes(type);
            
            if (facetTypes.empty()) return CellType::NODE; 

            if (type == CellType::PRISM6 || type == CellType::PRISM18) {
                // Les deux premières faces (0 et 1) sont des Triangles (index 0 dans cellFacetMapping)
                // Les trois faces suivantes (2, 3 et 4) sont des Quads (index 1 dans cellFacetMapping)
                return (localFacetId < 2) ? facetTypes[0] : facetTypes[1];
            }

            // Pour tous les autres types d'éléments (Tetra, Hexa, Tri, Quad), 
            // toutes les faces partagent le même type.
            return facetTypes[0];
        }

        inline static CfemInt getLocalFacesCount(CellType type)
        {
            switch (type)
            {
                case CellType::LINE2:
                case CellType::LINE3:
                    return 2;
                case CellType::TRI3:
                case CellType::TRI6:
                    return 3;
                case CellType::QUAD4:
                case CellType::QUAD9:
                    return 4;
                case CellType::TET4:
                case CellType::TET10:
                    return 4;
                case CellType::PRISM6:
                case CellType::PRISM18:
                    return 5;
                case CellType::HEX8:
                case CellType::HEX27:
                    return 6;

                default:
                    CFEM_ERROR << "Unsupported cell type."; 
                    return -1;
            }
        }

        inline static CfemInt getLocalEdgesCount(CellType type)
        {
            switch (type)
            {
                case CellType::LINE2:
                case CellType::LINE3:
                    return 1;
                case CellType::TRI3:
                case CellType::TRI6:
                    return 3;
                case CellType::QUAD4:
                case CellType::QUAD9:
                    return 4;
                case CellType::TET4:
                case CellType::TET10:
                    return 6;
                case CellType::PRISM6:
                case CellType::PRISM18:
                    return 9;
                case CellType::HEX8:
                case CellType::HEX27:
                    return 12;

                default:
                    CFEM_ERROR << "Unsupported cell type."; 
                    return -1;
            }
        }

        inline static CfemInt getLocalQuadFacesCount(CellType type)
        {
            switch (type)
            {
                case CellType::LINE2:
                case CellType::LINE3:
                case CellType::TRI3:
                case CellType::TRI6:
                case CellType::TET4:
                case CellType::TET10:
                    return 0;
                case CellType::QUAD4:
                case CellType::QUAD9:
                    return 1;
                case CellType::PRISM6:
                case CellType::PRISM18:
                    return 3;
                case CellType::HEX8:
                case CellType::HEX27:
                    return 6;
                default:
                    CFEM_ERROR << "Unsupported cell type."; 
                    return -1;
            }
        }

        inline static const std::vector<std::vector<CfemInt>>& getLocalEdgesIndices(CellType type)
        {
            static const std::vector<std::vector<CfemInt>> empty = {};
            static const std::vector<std::vector<CfemInt>> line2 = {{0, 1}};
            static const std::vector<std::vector<CfemInt>> line3 = {{0, 1, 2}};
            static const std::vector<std::vector<CfemInt>> tri3  = {{0, 1}, {1, 2}, {2, 0}};
            static const std::vector<std::vector<CfemInt>> tri6  = {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}};
            static const std::vector<std::vector<CfemInt>> quad4 = {{0, 1}, {1, 2}, {2, 3}, {3, 0}}; 
            static const std::vector<std::vector<CfemInt>> quad9 = {{0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7}};
            static const std::vector<std::vector<CfemInt>> tet4  = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};
            static const std::vector<std::vector<CfemInt>> tet10 = {{0, 1, 4}, {1, 2, 5}, {2, 0, 6}, {0, 3, 7}, {1, 3, 8}, {2, 3, 9}};
            static const std::vector<std::vector<CfemInt>> prism6 = {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}};
            static const std::vector<std::vector<CfemInt>> prism18 = {{0, 1, 6}, {1, 2, 7}, {2, 0, 8}, {3, 4, 9}, {4, 5, 10}, {5, 3, 11}, {0, 3, 12}, {1, 4, 13}, {2, 5, 14}};
            static const std::vector<std::vector<CfemInt>> hex8  = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7} };
            static const std::vector<std::vector<CfemInt>> hex27 = {{0, 1, 8}, {1, 2, 9}, {2, 3, 10}, {3, 0, 11}, {4, 5, 12}, {5, 6, 13}, {6, 7, 14}, {7, 4, 15}, {0, 4, 16}, {1, 5, 17}, {2, 6, 18}, {3, 7, 19} };

            switch (type)
            {
                case CellType::LINE2:   return line2;
                case CellType::LINE3:   return line3;
                case CellType::TRI3:    return tri3;
                case CellType::TRI6:    return tri6;
                case CellType::QUAD4:   return quad4;
                case CellType::QUAD9:   return quad9;
                case CellType::TET4:    return tet4;
                case CellType::TET10:   return tet10;
                case CellType::PRISM6:  return prism6;
                case CellType::PRISM18: return prism18;
                case CellType::HEX8:    return hex8;
                case CellType::HEX27:   return hex27;

                default:
                    CFEM_ERROR << "Unsupported cell type."; 
                    return empty;
            }
        }

        inline static const std::vector<std::vector<CfemInt>>& getLocalFacesIndices(CellType type)
        {
            static const std::vector<std::vector<CfemInt>> empty = {};
            
            // --- 1D: Les faces sont les noeuds d'extrémité (0D) ---
            static const std::vector<std::vector<CfemInt>> line2 = {{0}, {1}};
            static const std::vector<std::vector<CfemInt>> line3 = {{0}, {1}};
            
            // --- 2D: Les faces sont les arêtes (1D) ---
            static const std::vector<std::vector<CfemInt>> tri3  = {{0, 1}, {1, 2}, {2, 0}};
            static const std::vector<std::vector<CfemInt>> tri6  = {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}};
            static const std::vector<std::vector<CfemInt>> quad4 = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
            static const std::vector<std::vector<CfemInt>> quad9 = {{0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7}};

            // --- 3D: Les faces sont des surfaces (2D) ---
            
            static const std::vector<std::vector<CfemInt>> tet4  = {
                {0, 2, 1},
                {0, 1, 3}, 
                {1, 2, 3}, 
                {2, 0, 3}
            };
            
            static const std::vector<std::vector<CfemInt>> tet10 = {
                {0, 2, 1, 6, 5, 4}, 
                {0, 1, 3, 4, 8, 7}, 
                {1, 2, 3, 5, 9, 8}, 
                {2, 0, 3, 6, 7, 9}
            };
            // Prismes (2 Faces Triangulaires + 3 Faces Quadrilatérales)
            // L'ordre respecte ta méthode getFacetType (Triangles en index 0 et 1)
            static const std::vector<std::vector<CfemInt>> prism6 = {
                {0, 2, 1}, {3, 4, 5},                        // Faces 0, 1 : Triangles
                {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}     // Faces 2, 3, 4 : Quads
            };
            
            static const std::vector<std::vector<CfemInt>> prism18 = {
                {0, 2, 1, 8, 7, 6}, {3, 4, 5, 9, 10, 11},    // Faces 0, 1 : Triangles P2
                {0, 1, 4, 3, 6, 13, 9, 12, 15},              // Faces 2, 3, 4 : Quads P2
                {1, 2, 5, 4, 7, 14, 10, 13, 16}, 
                {2, 0, 3, 5, 8, 12, 11, 14, 17}
            };

            // Hexaèdres (Faces Quadrilatérales)
            static const std::vector<std::vector<CfemInt>> hex8 = {
                {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, 
                {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}
            };
            
            static const std::vector<std::vector<CfemInt>> hex27 = {
                {0, 3, 2, 1, 11, 10, 9, 8, 24}, 
                {4, 5, 6, 7, 12, 13, 14, 15, 25},
                {0, 1, 5, 4, 8, 17, 12, 16, 22},
                {1, 2, 6, 5, 9, 18, 13, 17, 21},
                {2, 3, 7, 6, 10, 19, 14, 18, 23},
                {3, 0, 4, 7, 11, 16, 15, 19, 20}
            };

            switch (type)
            {
                case CellType::LINE2:   return line2;
                case CellType::LINE3:   return line3;
                case CellType::TRI3:    return tri3;
                case CellType::TRI6:    return tri6;
                case CellType::QUAD4:   return quad4;
                case CellType::QUAD9:   return quad9;
                case CellType::TET4:    return tet4;
                case CellType::TET10:   return tet10;
                case CellType::PRISM6:  return prism6;
                case CellType::PRISM18: return prism18;
                case CellType::HEX8:    return hex8;
                case CellType::HEX27:   return hex27;

                default:
                    CFEM_ERROR << "Unsupported cell type."; 
                    return empty;
            }
        }

        inline static const std::vector<std::vector<CfemInt>>& getLocalQuadFacesIndices(CellType type)
        {
            static const std::vector<std::vector<CfemInt>> empty = {};
            static const std::vector<std::vector<CfemInt>> quad4 = {{0, 1, 2, 3}};
            static const std::vector<std::vector<CfemInt>> quad9 = {{0, 1, 2, 3, 4, 5, 6, 7, 8}};
            
            static const std::vector<std::vector<CfemInt>> prism6 = {
                {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}
            };
            
            static const std::vector<std::vector<CfemInt>> prism18 = {
                {0, 1, 4, 3, 6, 13, 9, 12, 15}, 
                {1, 2, 5, 4, 7, 14, 10, 13, 16}, 
                {2, 0, 3, 5, 8, 12, 11, 14, 17}
            };
            
            static const std::vector<std::vector<CfemInt>> hex8 = {
                {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, 
                {1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}
            };
            
            static const std::vector<std::vector<CfemInt>> hex27 = {
                {0, 3, 2, 1, 11, 10, 9, 8, 24}, 
                {4, 5, 6, 7, 12, 13, 14, 15, 25},
                {0, 1, 5, 4, 8, 17, 12, 16, 22},
                {1, 2, 6, 5, 9, 18, 13, 17, 21},
                {2, 3, 7, 6, 10, 19, 14, 18, 23},
                {3, 0, 4, 7, 11, 16, 15, 19, 20}
            };

            switch (type)
            {
                case CellType::LINE2: 
                case CellType::LINE3:
                case CellType::TRI3:
                case CellType::TRI6:
                case CellType::TET4:
                case CellType::TET10:
                    return empty;
                    
                case CellType::QUAD4:
                    return quad4;
                case CellType::QUAD9:
                    return quad9;
                    
                case CellType::PRISM6:
                    return prism6;
                case CellType::PRISM18:
                    return prism18;
                    
                case CellType::HEX8:
                    return hex8;
                case CellType::HEX27:
                    return hex27;

                default:
                    CFEM_ERROR << "Unsupported cell type."; 
                    return empty;
            }
        }

        /** ************************************************************************************************************************************
         * @brief Convert a cell type enum to a human-readable string.
         *
         * @param type The cell type (CellType enum).
         * @return A string representing the cell type, e.g., "TET4" or "QUAD9".
         *         Returns "UNKNOWN_CELL_TYPE" if the type is unrecognized.
         ***************************************************************************************************************************************/
        inline static std::string cellTypeToString(CellType type)
        {
            switch (type)
            {
            case CellType::LINE2:   return "LINE2";
            case CellType::TRI3:    return "TRI3";
            case CellType::QUAD4:   return "QUAD4";
            case CellType::TET4:    return "TET4";
            case CellType::PRISM6:  return "PRISM6";
            case CellType::HEX8:    return "HEX8";
            case CellType::LINE3:   return "LINE3";
            case CellType::TRI6:    return "TRI6";
            case CellType::QUAD9:   return "QUAD9";
            case CellType::TET10:   return "TET10";
            case CellType::PRISM18: return "PRISM18";
            case CellType::HEX27:   return "HEX27";
            // ... others?
            default:
                CFEM_ERROR << "UNKNOWN_CELL_TYPE";
                return "UNKNOWN_CELL_TYPE";
            }
        }

    };

    struct Facet {
        CfemInt id{-1};
        CellType type;
        CfemInt left_cell{-1};
        CfemInt right_cell{-1};

        CfemInt left_cell_facet_id{-1};  //! local id of the Facet in the left cell
        CfemInt right_cell_facet_id{-1}; //! local id of the Facet in the right ell
        
        // Index-based storage
        size_t offset{0};
        size_t node_count{0};

        bool isBoundary() const { return right_cell == -1; }
    };


    class MeshObserver {
    public:
        virtual ~MeshObserver() = default;
        virtual void updateAfterMeshChange() = 0;
    };

    class FESpaceObserver {
        std::string m_field_name;
    public:
        FESpaceObserver(const std::string& field_name): m_field_name(field_name) {}
        virtual ~FESpaceObserver() = default;
        const std::string& getFESpaceObserverName() const { return m_field_name; }
        virtual void updateOnSpaceChange() = 0;
    };


} // namespace
#endif