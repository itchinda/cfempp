#include "VTKExporter.h"
#include <algorithm>
#include <iostream>

namespace cfem {

void VTKExporter::write(const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open VTK file for writing: " + filename);
    }

    // --- 1. HEADER ---
    out << "# vtk DataFile Version 3.0\n";
    out << "CFEM Simulation Export\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n\n";

    // --- 2. POINTS (Vertices) ---
    out << "POINTS " << m_mesh.getNumVertices() << " double\n";
    for (size_t i = 0; i < m_mesh.getNumVertices(); ++i) {
        const auto& coords = m_mesh.getVertex(i).coords;
        out << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
    }
    out << "\n";

    // --- 3. CELLS ---
    size_t totalCellSize = 0;
    for (size_t i = 0; i < m_mesh.getNumCells(); ++i) {
        totalCellSize += 1 + m_mesh.getCell(i).vertex_ids.size();
    }

    out << "CELLS " << m_mesh.getNumCells() << " " << totalCellSize << "\n";
    for (size_t i = 0; i < m_mesh.getNumCells(); ++i) {
        const auto& cell = m_mesh.getCell(i);
        out << cell.vertex_ids.size();
        for (auto vid : cell.vertex_ids) out << " " << vid;
        out << "\n";
    }

    // --- 4. CELL TYPES ---
    out << "\nCELL_TYPES " << m_mesh.getNumCells() << "\n";
    for (size_t i = 0; i < m_mesh.getNumCells(); ++i) {
        out << getVTKCellType(m_mesh.getCell(i).type) << "\n";
    }

    // --- 5. POINT DATA ---
    bool pointDataHeaderWritten = false;
    auto ensurePointHeader = [&]() {
        if (!pointDataHeaderWritten) {
            out << "\nPOINT_DATA " << m_mesh.getNumVertices() << "\n";
            pointDataHeaderWritten = true;
        }
    };

    for (auto& field : m_fields) {
        ensurePointHeader();
        auto space = field->getSpace();
        const auto& data = field->data();
        if (space->getDimension() == 1) {
            out << "SCALARS " << field->getName() << " double 1\nLOOKUP_TABLE default\n";
            for (auto val : data) out << val << "\n";
        } 
        else if (space->getDimension() == 3) {
            out << "VECTORS " << field->getName() << " double\n";
            for (size_t i = 0; i < data.size(); i += 3) {
                out << data[i] << " " << data[i+1] << " " << data[i+2] << "\n";
            }
        }
    }

    for (const auto& entityName : m_exportedEntityNames) {
        MeshEntity* entity = const_cast<Mesh&>(m_mesh).getEntity(entityName);
        if (!entity || entity->getVertexIds().empty()) continue;
        ensurePointHeader();
        std::string safeName = "NodeMask_" + entityName;
        std::replace(safeName.begin(), safeName.end(), ' ', '_');
        out << "SCALARS " << safeName << " float 1\nLOOKUP_TABLE default\n";
        const auto& vids = entity->getVertexIds();
        for (CfemInt i = 0; i < (CfemInt)m_mesh.getNumVertices(); ++i) {
            out << (std::binary_search(vids.begin(), vids.end(), i) ? "1.0" : "0.0") << "\n";
        }
    }

    // --- 6. CELL DATA ---
    bool cellDataHeaderWritten = false;
    auto ensureCellHeader = [&]() {
        if (!cellDataHeaderWritten) {
            out << "\nCELL_DATA " << m_mesh.getNumCells() << "\n";
            cellDataHeaderWritten = true;
        }
    };

    // 6a. Entity cell masks
    for (const auto& entityName : m_exportedEntityNames) {
        MeshEntity* entity = const_cast<Mesh&>(m_mesh).getEntity(entityName);
        if (!entity || entity->getCellIds().empty()) continue;
        ensureCellHeader();
        std::string safeName = "CellMask_" + entityName;
        std::replace(safeName.begin(), safeName.end(), ' ', '_');
        out << "SCALARS " << safeName << " float 1\nLOOKUP_TABLE default\n";
        const auto& cids = entity->getCellIds();
        for (CfemInt i = 0; i < (CfemInt)m_mesh.getNumCells(); ++i) {
            out << (std::binary_search(cids.begin(), cids.end(), i) ? "1.0" : "0.0") << "\n";
        }
    }

    // 6b. Mesh Quality (New)
    if (m_exportQuality) {
        ensureCellHeader();
        out << "SCALARS Mesh_Quality float 1\n";
        out << "LOOKUP_TABLE default\n";
        for (size_t i = 0; i < m_mesh.getNumCells(); ++i) {
            out << m_mesh.computeCellQuality(i) << "\n";
        }
    }

    out.close();
}

int VTKExporter::getVTKCellType(CellType type) const {
    switch (type) {
        case CellType::LINE2:   return 3;  
        case CellType::LINE3:   return 21; 
        case CellType::TRI3:    return 5;  
        case CellType::TRI6:    return 22; 
        case CellType::QUAD4:   return 9;  
        case CellType::QUAD9:   return 28; 
        case CellType::TET4:    return 10; 
        case CellType::TET10:   return 24; 
        case CellType::HEX8:    return 12; 
        case CellType::HEX27:   return 29; 
        case CellType::PRISM6:  return 13; 
        case CellType::PRISM18: return 32; 
        default: return 0;
    }
}

}