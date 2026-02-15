/**
 * @file VTKExporter.cpp
 * @brief Implementation of the VTK Legacy format exporter for CFEM simulations.
 */

#include "VTKExporter.h"
#include "LagrangeSpace.h"
#include "ReferenceElementFactory.h"
#include <algorithm>
#include <iostream>
#include <iomanip>

namespace cfem {

/**
 * @brief Exports the mesh and associated fields to a VTK Legacy (.vtk) file.
 * * This method follows the VTK Unstructured Grid specification. It handles:
 * 1. Geometry (Points and Cells)
 * 2. Field Interpolation (Mapping P1/P2 physics to geometric vertices)
 * 3. Metadata (Mesh quality and Entity masks)
 * * @param filename Path to the output file.
 * @throws std::runtime_error If the file cannot be opened for writing.
 */
void VTKExporter::write(const std::string& filename) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open VTK file for writing: " + filename);
    }

    // High precision is required to avoid artifacts in derivative calculations in ParaView
    out << std::scientific << std::setprecision(16);

    // --- 1. HEADER ---
    out << "# vtk DataFile Version 3.0\n";
    out << "CFEM Simulation Export\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n\n";

    // --- 2. POINTS ---
    out << "POINTS " << m_mesh.getNumVertices() << " double\n";
    for (size_t i = 0; i < m_mesh.getNumVertices(); ++i) {
        const auto& coords = m_mesh.getVertex(i).coords;
        out << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
    }

    // --- 3. CELLS ---
    // Calculate total size: sum of (1 + number of nodes per cell)
    size_t totalCellSize = 0;
    for (size_t i = 0; i < m_mesh.getNumCells(); ++i) {
        totalCellSize += 1 + m_mesh.getCellConnectivity(i).size();
    }

    out << "\nCELLS " << m_mesh.getNumCells() << " " << totalCellSize << "\n";
    for (size_t i = 0; i < m_mesh.getNumCells(); ++i) {
        const auto& cell_ids = m_mesh.getCellConnectivity(i);
        out << cell_ids.size();
        for (auto vid : cell_ids) out << " " << vid;
        out << "\n";
    }

    // --- 4. CELL TYPES ---
    out << "\nCELL_TYPES " << m_mesh.getNumCells() << "\n";
    for (size_t i = 0; i < m_mesh.getNumCells(); ++i) {
        out << getVTKCellType(m_mesh.getCell(i).type) << "\n";
    }

    // --- 5. POINT DATA (Nodal Interpolation) ---
    if (!m_fields.empty() || !m_exportedEntityNames.empty()) {
        out << "\nPOINT_DATA " << m_mesh.getNumVertices() << "\n";

        for (const auto* field : m_fields) {
            /** @note Sync PETSc global data to CPU local buffer before interpolation */
            field->updateLocalBuffer(); 
            
            auto space = std::dynamic_pointer_cast<LagrangeSpace>(field->getSpace());
            int dim = field->getDimension(); 
            
            // Temporary buffers for interpolated values
            std::vector<double> interpolatedValues(m_mesh.getNumVertices() * dim, 0.0);
            std::vector<bool> nodeComputed(m_mesh.getNumVertices(), false);

            /** * Loop over cells to interpolate the physical field at geometric node locations.
             * This ensures that even if the field is P1 and the mesh is P2, 
             * every geometric vertex (including mid-edge nodes) receives a valid value.
             */
            for (CfemInt c = 0; c < (CfemInt)m_mesh.getNumCells(); ++c) {
                const auto& cell = m_mesh.getCell(c);
                const auto& cell_ids =  m_mesh.getCellConnectivity(c);
                const auto& refEl = space->getRefElement(cell.type);
                const auto& xiNodes = refEl.referenceNodes(); 

                for (size_t n = 0; n < cell_ids.size(); ++n) {
                    CfemInt vid = cell_ids[n];
                    if (nodeComputed[vid]) continue;

                    // Evaluate the field at the local xi coordinate of the mesh node
                    std::vector<double> val;
                    field->evaluateAt(c, xiNodes[n], val);

                    for (int d = 0; d < dim; ++d) {
                        interpolatedValues[vid * dim + d] = val[d];
                    }
                    nodeComputed[vid] = true;
                }
            }

            // Write Field Scalars or Vectors
            if (dim == 1) {
                out << "SCALARS " << field->getName() << " double 1\nLOOKUP_TABLE default\n";
                for (size_t i = 0; i < m_mesh.getNumVertices(); ++i) {
                    out << interpolatedValues[i] << "\n";
                }
            } else {
                out << "VECTORS " << field->getName() << " double\n";
                for (size_t i = 0; i < m_mesh.getNumVertices(); ++i) {
                    double vx = interpolatedValues[i * dim + 0];
                    double vy = (dim >= 2) ? interpolatedValues[i * dim + 1] : 0.0;
                    double vz = (dim >= 3) ? interpolatedValues[i * dim + 2] : 0.0;
                    out << vx << " " << vy << " " << vz << "\n";
                }
            }
        }

        // --- Node Masks ---
        for (const auto& entityName : m_exportedEntityNames) {
            MeshEntity* entity = const_cast<Mesh&>(m_mesh).getEntity(entityName);
            if (!entity) continue;
            
            std::string safeName = "NodeMask_" + entityName;
            std::replace(safeName.begin(), safeName.end(), ' ', '_');
            out << "SCALARS " << safeName << " float 1\nLOOKUP_TABLE default\n";
            
            const auto& vids = entity->getVertexIds();
            for (CfemInt i = 0; i < (CfemInt)m_mesh.getNumVertices(); ++i) {
                // Use binary search on sorted vertex IDs for O(log n) lookup
                out << (std::binary_search(vids.begin(), vids.end(), i) ? "1.0" : "0.0") << "\n";
            }
        }
    }

    // --- 6. CELL DATA ---
    if (m_exportQuality || !m_exportedEntityNames.empty()) {
        out << "\nCELL_DATA " << m_mesh.getNumCells() << "\n";

        if (m_exportQuality) {
            out << "SCALARS Mesh_Quality float 1\nLOOKUP_TABLE default\n";
            for (size_t i = 0; i < m_mesh.getNumCells(); ++i) {
                out << m_mesh.computeCellQuality(i) << "\n";
            }
        }

        for (const auto& entityName : m_exportedEntityNames) {
            MeshEntity* entity = const_cast<Mesh&>(m_mesh).getEntity(entityName);
            if (!entity) continue;
            
            std::string safeName = "CellMask_" + entityName;
            std::replace(safeName.begin(), safeName.end(), ' ', '_');
            out << "SCALARS " << safeName << " float 1\nLOOKUP_TABLE default\n";
            
            const auto& cids = entity->getCellIds();
            for (CfemInt i = 0; i < (CfemInt)m_mesh.getNumCells(); ++i) {
                out << (std::binary_search(cids.begin(), cids.end(), i) ? "1.0" : "0.0") << "\n";
            }
        }
    }

    out.close();
}

/**
 * @brief Maps internal CellType to VTK legacy cell identifiers.
 * @param type The internal CFEM CellType.
 * @return The corresponding VTK integer ID.
 */
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
        default: CFEM_ERROR("Unsupported cell type."); // VTK_EMPTY_CELL case 13: return CellType::PRISM6; case 32: return CellType::PRISM18;
    }
}

} // namespace cfem