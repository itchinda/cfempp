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

// #include "ReferenceElementFactory.h"
// #include "ReferenceSegment.h"
// #include "ReferenceTriangle.h"
// #include "ReferenceQuad.h"
// #include "ReferenceTetrahedron.h"
// #include "ReferencePrism.h"
// #include "ReferenceHexahedron.h"


// namespace cfem 
// {

//     const ReferenceElement& ReferenceElementFactory::getReferenceElement(CellType type, int order) {
//         switch (type) {
//             case CellType::LINE2:
//             case CellType::LINE3:
//             {
//                 if (order == 0) {
//                     static P0ReferenceSegment instance;
//                     return instance;
//                 } else if (order == 1) {
//                     static P1ReferenceSegment instance;
//                     return instance;
//                 } else if (order == 2) {
//                     static P2ReferenceSegment instance;
//                     return instance;
//                 }
//                 else{
//                     CFEM_ERROR << "Element order not supported in Line (segment)";
//                 }
//                 break;
//             }
//             case CellType::TRI3:
//             case CellType::TRI6:
//             {
//                 if (order == 0) {
//                     static P0ReferenceTriangle instance;
//                     return instance;
//                 } else if (order == 1) {
//                     static P1ReferenceTriangle instance;
//                     return instance;
//                 } else if (order == 2) {
//                     static P2ReferenceTriangle instance;
//                     return instance;
//                 }
//                 else{
//                     CFEM_ERROR << "Order " << order << " not supported for Triangle";
//                 }
//                 break;
//             }
//             case CellType::QUAD4:
//             case CellType::QUAD9:
//             {
//                 if (order == 0) {
//                     static Q0ReferenceQuad instance;
//                     return instance;
//                 } else if (order == 1) {
//                     static Q1ReferenceQuad instance;
//                     return instance;
//                 } else if (order == 2) {
//                     static Q2ReferenceQuad instance;
//                     return instance;
//                 }
//                 else{
//                     CFEM_ERROR << "Order " << order << " not supported for Quad";
//                 }
//             }
//             case CellType::TET4: 
//             case CellType::TET10:
//             {
//                 if (order == 0) {
//                     static P0ReferenceTetrahedron instance;
//                     return instance;
//                 } else if (order == 1) {
//                     static P1ReferenceTetrahedron instance;
//                     return instance;
//                 } else if (order == 2) {
//                     static P2ReferenceTetrahedron instance;
//                     return instance;
//                 }
//                 else{
//                     CFEM_ERROR << "Order " << order << " not supported for Tetrahedron";
//                 }
//                 break;
//             }
//             case CellType::PRISM6: 
//             case CellType::PRISM18:
//             {
//                 if (order == 0) {
//                     static P0ReferencePrism instance;
//                     return instance;
//                 } else if (order == 1) {
//                     static P1ReferencePrism instance;
//                     return instance;
//                 } else if (order == 2) {
//                     static P2ReferencePrism instance;
//                     return instance;
//                 }
//                 else{
//                     CFEM_ERROR << "Order " << order << " not supported for Edge (prism)";
//                 }
//                 break;
//             }
//             case CellType::HEX8:
//             case CellType::HEX27:
//             {
//                 if (order == 0) {
//                     static Q0ReferenceHexahedron instance;
//                     return instance;
//                 } else if (order == 1) {
//                     static Q1ReferenceHexahedron instance;
//                     return instance;
//                 } else if (order == 2) {
//                     static Q2ReferenceHexahedron instance;
//                     return instance;
//                 }
//                 else{
//                     CFEM_ERROR << "Order " << order << " not supported for Hexahedron";
//                 }
//                 break;
//             }
//             default: 
//                 CFEM_ERROR << "CellType not supported in ReferenceElementFactory";
//         }
//         throw std::invalid_argument("Unsupported CellType");// Just to make the compiler hapy. Should never reach here
//     }
// } // namespace cfem


#include "ReferenceElementFactory.h"
#include "ReferenceElements.h"
#include "FESpace.h"

namespace cfem {

namespace constants{
    /** 
     * @brief Constants for the lookup table dimensions.
     * We use a fixed size based on the CellType enum max value.
     */
    constexpr int MAX_CELL_TYPES = 30; 
    constexpr int MAX_SUPPORTED_ORDERS = 4; // P0, P1, P2, P1Mini
}

int ReferenceElementFactory::orderToIndex(int order) {
    if (order >= 0 && order <= 2) return order; // 0, 1, 2
    if (order == static_cast<int>(Order::Mini)) return 3; 
    return -1;
}

const ReferenceElement& ReferenceElementFactory::getReferenceElement(CellType type, int order) {
    // 1. Static table of pointers initialized to nullptr
    static const ReferenceElement* s_lookupTable[constants::MAX_CELL_TYPES][constants::MAX_SUPPORTED_ORDERS] = { nullptr };
    static bool s_isInitialized = false;

    // 2. Thread-safe initialization of the table (happens only once)
    if (!s_isInitialized) {
        // --- SEGMENTS ---
        int seg = static_cast<int>(CellType::LINE2);
        s_lookupTable[seg][0] = &singleton<P0ReferenceSegment>();
        s_lookupTable[seg][1] = &singleton<P1ReferenceSegment>();
        s_lookupTable[seg][2] = &singleton<P2ReferenceSegment>();
        s_lookupTable[seg][3] = &singleton<P2ReferenceSegment>(); // P1Mini on Line = P2

        // --- TRIANGLES ---
        int tri = static_cast<int>(CellType::TRI3);
        s_lookupTable[tri][0] = &singleton<P0ReferenceTriangle>();
        s_lookupTable[tri][1] = &singleton<P1ReferenceTriangle>();
        s_lookupTable[tri][2] = &singleton<P2ReferenceTriangle>();
        s_lookupTable[tri][3] = &singleton<MiniReferenceTriangle>();

        // --- QUADS ---
        int quad = static_cast<int>(CellType::QUAD4);
        s_lookupTable[quad][0] = &singleton<Q0ReferenceQuad>();
        s_lookupTable[quad][1] = &singleton<Q1ReferenceQuad>();
        s_lookupTable[quad][2] = &singleton<Q2ReferenceQuad>();
        s_lookupTable[quad][3] = &singleton<MiniReferenceQuad>();

        // --- TETRAHEDRA ---
        int tet = static_cast<int>(CellType::TET4);
        s_lookupTable[tet][0] = &singleton<P0ReferenceTetrahedron>();
        s_lookupTable[tet][1] = &singleton<P1ReferenceTetrahedron>();
        s_lookupTable[tet][2] = &singleton<P2ReferenceTetrahedron>();
        s_lookupTable[tet][3] = &singleton<MiniReferenceTetrahedron>();

        // --- HEXAHEDRA ---
        int hex = static_cast<int>(CellType::HEX8);
        s_lookupTable[hex][0] = &singleton<Q0ReferenceHexahedron>();
        s_lookupTable[hex][1] = &singleton<Q1ReferenceHexahedron>();
        s_lookupTable[hex][2] = &singleton<Q2ReferenceHexahedron>();
        s_lookupTable[hex][3] = &singleton<MiniReferenceHexahedron>();

        // --- PRISMS ---
        int prism = static_cast<int>(CellType::PRISM6);
        s_lookupTable[prism][0] = &singleton<P0ReferencePrism>();
        s_lookupTable[prism][1] = &singleton<P1ReferencePrism>();
        s_lookupTable[prism][2] = &singleton<P2ReferencePrism>();
        s_lookupTable[prism][3] = &singleton<MiniReferencePrism>();

        // Support for quadratic topology types (TRI6, TET10, etc.)
        // They point to the same reference elements as their linear counterparts
        // --- Mapping des types quadratiques (TRI6, TET10, etc.) ---
        auto copyRow = [&](CellType source, CellType target) {
            for (int i = 0; i < constants::MAX_SUPPORTED_ORDERS; ++i) {
                s_lookupTable[static_cast<int>(target)][i] = s_lookupTable[static_cast<int>(source)][i];
            }
        };

        copyRow(CellType::LINE2,  CellType::LINE3);
        copyRow(CellType::TRI3,   CellType::TRI6);
        copyRow(CellType::QUAD4,  CellType::QUAD9);
        copyRow(CellType::TET4,   CellType::TET10);
        copyRow(CellType::HEX8,   CellType::HEX27);
        copyRow(CellType::PRISM6, CellType::PRISM18);

        s_isInitialized = true;
    }

    // 3. Fast Lookup
    const int typeIdx = static_cast<int>(type);
    const int ordIdx = orderToIndex(order);

    if (typeIdx >= 0 && typeIdx < constants::MAX_CELL_TYPES && ordIdx >= 0) {
        const ReferenceElement* element = s_lookupTable[typeIdx][ordIdx];
        if (element) {
            return *element;
        }
    }

    // 4. Error handling
    CFEM_ERROR << "ReferenceElementFactory: Unsupported combination [Type: " 
               << typeIdx << ", Order: " << order << "]";
               
    throw std::invalid_argument("Requested element type or order is not implemented in the Factory.");
}

} // namespace cfem