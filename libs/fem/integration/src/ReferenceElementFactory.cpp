#include "ReferenceElementFactory.h"
#include "ReferenceSegment.h"
#include "ReferenceTriangle.h"
#include "ReferenceQuad.h"
#include "ReferenceTetrahedron.h"
#include "ReferencePrism.h"
#include "ReferenceHexahedron.h"


namespace cfem 
{

    const ReferenceElement& ReferenceElementFactory::getReferenceElement(CellType type, int order) {
        switch (type) {
            case CellType::LINE2:
            case CellType::LINE3:
            {
                if (order == 0) {
                    static P0ReferenceSegment instance;
                    return instance;
                } else if (order == 1) {
                    static P1ReferenceSegment instance;
                    return instance;
                } else if (order == 2) {
                    static P2ReferenceSegment instance;
                    return instance;
                }
                else{
                    CFEM_ERROR << "Element order not supported in Line (segment)";
                }
                break;
            }
            case CellType::TRI3:
            case CellType::TRI6:
            {
                if (order == 0) {
                    static P0ReferenceTriangle instance;
                    return instance;
                } else if (order == 1) {
                    static P1ReferenceTriangle instance;
                    return instance;
                } else if (order == 2) {
                    static P2ReferenceTriangle instance;
                    return instance;
                }
                else{
                    CFEM_ERROR << "Order " << order << " not supported for Triangle";
                }
                break;
            }
            case CellType::QUAD4:
            case CellType::QUAD9:
            {
                if (order == 0) {
                    static Q0ReferenceQuad instance;
                    return instance;
                } else if (order == 1) {
                    static Q1ReferenceQuad instance;
                    return instance;
                } else if (order == 2) {
                    static Q2ReferenceQuad instance;
                    return instance;
                }
                else{
                    CFEM_ERROR << "Order " << order << " not supported for Quad";
                }
            }
            case CellType::TET4: 
            case CellType::TET10:
            {
                if (order == 0) {
                    static P0ReferenceTetrahedron instance;
                    return instance;
                } else if (order == 1) {
                    static P1ReferenceTetrahedron instance;
                    return instance;
                } else if (order == 2) {
                    static P2ReferenceTetrahedron instance;
                    return instance;
                }
                else{
                    CFEM_ERROR << "Order " << order << " not supported for Tetrahedron";
                }
                break;
            }
            case CellType::PRISM6: 
            case CellType::PRISM18:
            {
                if (order == 0) {
                    static P0ReferencePrism instance;
                    return instance;
                } else if (order == 1) {
                    static P1ReferencePrism instance;
                    return instance;
                } else if (order == 2) {
                    static P2ReferencePrism instance;
                    return instance;
                }
                else{
                    CFEM_ERROR << "Order " << order << " not supported for Edge (prism)";
                }
                break;
            }
            case CellType::HEX8:
            case CellType::HEX27:
            {
                if (order == 0) {
                    static Q0ReferenceHexahedron instance;
                    return instance;
                } else if (order == 1) {
                    static Q1ReferenceHexahedron instance;
                    return instance;
                } else if (order == 2) {
                    static Q2ReferenceHexahedron instance;
                    return instance;
                }
                else{
                    CFEM_ERROR << "Order " << order << " not supported for Hexahedron";
                }
                break;
            }
            default: 
                CFEM_ERROR << "CellType not supported in ReferenceElementFactory";
        }
        throw std::invalid_argument("Unsupported CellType");// Just to make the compiler hapy. Should never reach here
    }
} // namespace cfem