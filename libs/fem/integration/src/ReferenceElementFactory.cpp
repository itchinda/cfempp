#include "ReferenceElementFactory.h"
#include "ReferenceSegment.h"
#include "ReferenceTriangle.h"
#include "ReferenceQuad.h"
#include "ReferenceTetrahedron.h"
#include "ReferencePrism.h"
#include "ReferenceHexahedron.h"


namespace cfem {

// std::map<CellType, std::unique_ptr<ReferenceElement>> ReferenceElementFactory::m_instances;

// const ReferenceElement& ReferenceElementFactory::getReferenceElement(CellType type) {
//     auto it = m_instances.find(type);
    
//     if (it == m_instances.end()) {
//         switch (type) {
//             case CellType::TRI3:  m_instances[type] = std::make_unique<P1ReferenceTriangle>(); break;
//             case CellType::TRI6:  m_instances[type] = std::make_unique<P2ReferenceTriangle>(); break;
//             case CellType::QUAD4: m_instances[type] = std::make_unique<Q1ReferenceQuad>(); break;
//             case CellType::QUAD9: m_instances[type] = std::make_unique<Q2ReferenceQuad>(); break;
//             case CellType::TET4: m_instances[type] = std::make_unique<P1ReferenceTetrahedron>(); break;
//             case CellType::TET10: m_instances[type] = std::make_unique<P2ReferenceTetrahedron>(); break;
//             case CellType::PRISM6: m_instances[type] = std::make_unique<P1ReferencePrism>(); break;
//             case CellType::PRISM18: m_instances[type] = std::make_unique<P2ReferencePrism>(); break;
//             case CellType::HEX8: m_instances[type] = std::make_unique<Q1ReferenceHexahedron>(); break;
//             case CellType::HEX27: m_instances[type] = std::make_unique<Q2ReferenceHexahedron>(); break;
//             default:
//                 CFEM_ERROR << "CellType not supported in ReferenceElementFactory";
//         }
//         return *m_instances[type];
//     }
    
//     return *(it->second);
// }


const ReferenceElement& ReferenceElementFactory::getReferenceElement(CellType type) {
    switch (type) {
        case CellType::LINE2: {
            static P1ReferenceSegment instance;
            return instance;
        }
        case CellType::LINE3: {
            static P2ReferenceSegment instance;
            return instance;
        }
        case CellType::TRI3: {
            static P1ReferenceTriangle instance;
            return instance;
        }
        case CellType::TRI6: {
            static P2ReferenceTriangle instance;
            return instance;
        }
        case CellType::QUAD4: {
            static Q1ReferenceQuad instance;
            return instance;
        }
        case CellType::QUAD9: {
            static Q2ReferenceQuad instance;
            return instance;
        }
        case CellType::TET4: {
            static P1ReferenceTetrahedron instance;
            return instance;
        }
        case CellType::TET10: {
            static P2ReferenceTetrahedron instance;
            return instance;
        }
        case CellType::PRISM6: {
            static P1ReferencePrism instance;
            return instance;
        }
        case CellType::PRISM18: {
            static P2ReferencePrism instance;
            return instance;
        }
        case CellType::HEX8: {
            static Q1ReferenceHexahedron instance;
            return instance;
        }
        case CellType::HEX27: {
            static Q2ReferenceHexahedron instance;
            return instance;
        }
        default:
            CFEM_ERROR << "CellType not supported in ReferenceElementFactory";
            // Juste pour le compilateur, il n'atteindra jamais ce point
            throw std::runtime_error("Unsupported CellType");
    }
}

} // namespace cfem