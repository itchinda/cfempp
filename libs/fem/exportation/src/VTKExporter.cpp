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

#include "VTKExporter.h"
#include "CfemExpression.h"
#include <algorithm>
#include <iomanip>
#include <fstream>

namespace cfem
{
    VTKExporter::VTKExporter(const Mesh &mesh) : VTKFamilyExporter(mesh)
    {  
        m_outputFileExtension = "vtk";
    }

    // ========================================================================
    // --- FORMATTER 1: ASCII CONTINUOUS
    // ========================================================================
    void VTKExporter::writeASCIIContinuous(const std::string &filename, 
                                           const std::vector<ExportEntry>& fields, 
                                           CfemInt entityId) const
    {
        std::ofstream out(filename);
        if (!out.is_open()) CFEM_ERROR << "Cannot open VTK file: " << filename;

        auto geo = buildContinuousGeometry(entityId);
        auto data = evaluateFields(geo, fields, entityId, false);

        out << "# vtk DataFile Version 3.0\nCFEM++ Continuous Export\nASCII\nDATASET UNSTRUCTURED_GRID\n";
        
        // --- POINTS ---
        out << "POINTS " << geo.numActiveVerts << " double\n";
        for (size_t i = 0; i < geo.numActiveVerts; ++i) {
            out << std::setprecision(16) << geo.points[3*i] << " " << geo.points[3*i+1] << " " << geo.points[3*i+2] << "\n";
        }

        // --- CELLS ---
        int totalConnSize = geo.connectivity.size() + geo.numActiveCells;
        out << "\nCELLS " << geo.numActiveCells << " " << totalConnSize << "\n";
        for (size_t i = 0; i < geo.numActiveCells; ++i) {
            int start = (i == 0) ? 0 : geo.offsets[i-1];
            int end = geo.offsets[i];
            out << (end - start);
            for (int j = start; j < end; ++j) out << " " << geo.connectivity[j];
            out << "\n";
        }

        out << "\nCELL_TYPES " << geo.numActiveCells << "\n";
        for (auto t : geo.types) out << static_cast<int>(t) << "\n";

        // --- POINT DATA (Fields & Node Masks) ---
        if (!fields.empty() || !data.nodeMasks.empty()) {
            out << "\nPOINT_DATA " << geo.numActiveVerts << "\n";
            for (size_t i = 0; i < fields.size(); ++i) {
                int dim = fields[i].expr->getDimension();
                if (dim == 1) {
                    out << "SCALARS " << fields[i].name << " double 1\nLOOKUP_TABLE default\n";
                    for (size_t v = 0; v < geo.numActiveVerts; ++v) out << data.fields[i][v] << "\n";
                } else {
                    out << "VECTORS " << fields[i].name << " double\n";
                    for (size_t v = 0; v < geo.numActiveVerts; ++v) {
                        out << data.fields[i][v * dim + 0] << " " << data.fields[i][v * dim + 1] << " " 
                            << (dim >= 3 ? data.fields[i][v * dim + 2] : 0.0) << "\n";
                    }
                }
            }
            // Node Masks
            for (size_t i = 0; i < data.nodeMasks.size(); ++i) {
                out << "SCALARS NodeMask_" << m_exportedEntityNames[i] << " float 1\nLOOKUP_TABLE default\n";
                for (size_t v = 0; v < geo.numActiveVerts; ++v) out << data.nodeMasks[i][v] << "\n";
            }
        }

        // --- CELL DATA (Quality & Cell Masks) ---
        if (m_exportQuality || !data.cellMasks.empty()) {
            out << "\nCELL_DATA " << geo.numActiveCells << "\n";
            if (m_exportQuality) {
                out << "SCALARS Mesh_Quality float 1\nLOOKUP_TABLE default\n";
                for (auto q : data.cellQuality) out << q << "\n";
            }
            for (size_t i = 0; i < data.cellMasks.size(); ++i) {
                out << "SCALARS CellMask_" << m_exportedEntityNames[i] << " float 1\nLOOKUP_TABLE default\n";
                for (auto m : data.cellMasks[i]) out << m << "\n";
            }
        }
    }

    // ========================================================================
    // --- FORMATTER 2: BINARY CONTINUOUS
    // ========================================================================
    void VTKExporter::writeBINARYContinuous(const std::string &filename, 
                                            const std::vector<ExportEntry>& fields, 
                                            CfemInt entityId) const
    {
        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open()) CFEM_ERROR << "Cannot open VTK file: " << filename;

        // Big-Endian serializer constraint for VTK legacy
        auto writeBigEndian = [&](const auto *d, std::size_t count) {
            using T = std::decay_t<decltype(d[0])>;
            for (std::size_t i = 0; i < count; ++i) {
                T val = d[i];
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
                char *p = reinterpret_cast<char *>(&val);
                std::reverse(p, p + sizeof(T));
#endif
                out.write(reinterpret_cast<const char *>(&val), sizeof(T));
            }
        };

        auto geo = buildContinuousGeometry(entityId);
        auto data = evaluateFields(geo, fields, entityId, false);

        out << "# vtk DataFile Version 3.0\nCFEM++ Continuous Export\nBINARY\nDATASET UNSTRUCTURED_GRID\n";
        
        // --- POINTS ---
        out << "POINTS " << geo.numActiveVerts << " double\n";
        writeBigEndian(geo.points.data(), geo.points.size());
        out << "\n";

        // --- CELLS ---
        int totalConnSize = geo.connectivity.size() + geo.numActiveCells;
        out << "CELLS " << geo.numActiveCells << " " << totalConnSize << "\n";
        std::vector<int> cellBuffer;
        cellBuffer.reserve(totalConnSize);
        for (size_t i = 0; i < geo.numActiveCells; ++i) {
            int start = (i == 0) ? 0 : geo.offsets[i-1];
            int end = geo.offsets[i];
            cellBuffer.push_back(end - start);
            for (int j = start; j < end; ++j) cellBuffer.push_back(geo.connectivity[j]);
        }
        writeBigEndian(cellBuffer.data(), cellBuffer.size());
        out << "\n";

        out << "CELL_TYPES " << geo.numActiveCells << "\n";
        std::vector<int> typeBuffer(geo.types.begin(), geo.types.end());
        writeBigEndian(typeBuffer.data(), typeBuffer.size());
        out << "\n";

        // --- POINT DATA ---
        if (!fields.empty() || !data.nodeMasks.empty()) {
            out << "POINT_DATA " << geo.numActiveVerts << "\n";
            for (size_t i = 0; i < fields.size(); ++i) {
                int dim = fields[i].expr->getDimension();
                if (dim == 1) {
                    out << "SCALARS " << fields[i].name << " double 1\nLOOKUP_TABLE default\n";
                    writeBigEndian(data.fields[i].data(), data.fields[i].size());
                    out << "\n";
                } else {
                    out << "VECTORS " << fields[i].name << " double\n";
                    std::vector<double> vecBuffer(3 * geo.numActiveVerts, 0.0);
                    for (size_t v = 0; v < geo.numActiveVerts; ++v) {
                        vecBuffer[3*v+0] = data.fields[i][v*dim+0];
                        vecBuffer[3*v+1] = data.fields[i][v*dim+1];
                        if (dim >= 3) vecBuffer[3*v+2] = data.fields[i][v*dim+2];
                    }
                    writeBigEndian(vecBuffer.data(), vecBuffer.size());
                    out << "\n";
                }
            }
            // Node Masks
            for (size_t i = 0; i < data.nodeMasks.size(); ++i) {
                out << "SCALARS NodeMask_" << m_exportedEntityNames[i] << " float 1\nLOOKUP_TABLE default\n";
                writeBigEndian(data.nodeMasks[i].data(), data.nodeMasks[i].size());
                out << "\n";
            }
        }

        // --- CELL DATA ---
        if (m_exportQuality || !data.cellMasks.empty()) {
            out << "CELL_DATA " << geo.numActiveCells << "\n";
            if (m_exportQuality) {
                out << "SCALARS Mesh_Quality float 1\nLOOKUP_TABLE default\n";
                writeBigEndian(data.cellQuality.data(), data.cellQuality.size());
                out << "\n";
            }
            for (size_t i = 0; i < data.cellMasks.size(); ++i) {
                out << "SCALARS CellMask_" << m_exportedEntityNames[i] << " float 1\nLOOKUP_TABLE default\n";
                writeBigEndian(data.cellMasks[i].data(), data.cellMasks[i].size());
                out << "\n";
            }
        }
    }

    // ========================================================================
    // --- FORMATTER 3 & 4: DISCONTINUOUS
    // ========================================================================
    // Grâce à l'abstraction, la logique de formatage est rigoureusement identique.
    // On rappelle les noyaux en activant le mode "discontinu".
    
    void VTKExporter::writeASCIIDiscontinuous(const std::string &filename, 
                                              const std::vector<ExportEntry>& fields, 
                                              CfemInt entityId) const
    {
        std::ofstream out(filename);
        if (!out.is_open()) CFEM_ERROR << "Cannot open VTK file: " << filename;

        auto geo = buildDiscontinuousGeometry(entityId);
        auto data = evaluateFields(geo, fields, entityId, true);

        out << "# vtk DataFile Version 3.0\nCFEM++ Discontinuous Export\nASCII\nDATASET UNSTRUCTURED_GRID\n";
        
        // --- POINTS ---
        out << "POINTS " << geo.numActiveVerts << " double\n";
        for (size_t i = 0; i < geo.numActiveVerts; ++i) {
            out << std::setprecision(16) << geo.points[3*i] << " " << geo.points[3*i+1] << " " << geo.points[3*i+2] << "\n";
        }

        // --- CELLS ---
        int totalConnSize = geo.connectivity.size() + geo.numActiveCells;
        out << "\nCELLS " << geo.numActiveCells << " " << totalConnSize << "\n";
        for (size_t i = 0; i < geo.numActiveCells; ++i) {
            int start = (i == 0) ? 0 : geo.offsets[i-1];
            int end = geo.offsets[i];
            out << (end - start);
            for (int j = start; j < end; ++j) out << " " << geo.connectivity[j];
            out << "\n";
        }

        out << "\nCELL_TYPES " << geo.numActiveCells << "\n";
        for (auto t : geo.types) out << static_cast<int>(t) << "\n";

        // --- POINT DATA ---
        if (!fields.empty() || !data.nodeMasks.empty()) {
            out << "\nPOINT_DATA " << geo.numActiveVerts << "\n";
            for (size_t i = 0; i < fields.size(); ++i) {
                int dim = fields[i].expr->getDimension();
                if (dim == 1) {
                    out << "SCALARS " << fields[i].name << " double 1\nLOOKUP_TABLE default\n";
                    for (size_t v = 0; v < geo.numActiveVerts; ++v) out << data.fields[i][v] << "\n";
                } else {
                    out << "VECTORS " << fields[i].name << " double\n";
                    for (size_t v = 0; v < geo.numActiveVerts; ++v) {
                        out << data.fields[i][v * dim + 0] << " " << data.fields[i][v * dim + 1] << " " 
                            << (dim >= 3 ? data.fields[i][v * dim + 2] : 0.0) << "\n";
                    }
                }
            }
            for (size_t i = 0; i < data.nodeMasks.size(); ++i) {
                out << "SCALARS NodeMask_" << m_exportedEntityNames[i] << " float 1\nLOOKUP_TABLE default\n";
                for (size_t v = 0; v < geo.numActiveVerts; ++v) out << data.nodeMasks[i][v] << "\n";
            }
        }

        // --- CELL DATA ---
        if (m_exportQuality || !data.cellMasks.empty()) {
            out << "\nCELL_DATA " << geo.numActiveCells << "\n";
            if (m_exportQuality) {
                out << "SCALARS Mesh_Quality float 1\nLOOKUP_TABLE default\n";
                for (auto q : data.cellQuality) out << q << "\n";
            }
            for (size_t i = 0; i < data.cellMasks.size(); ++i) {
                out << "SCALARS CellMask_" << m_exportedEntityNames[i] << " float 1\nLOOKUP_TABLE default\n";
                for (auto m : data.cellMasks[i]) out << m << "\n";
            }
        }
    }

    void VTKExporter::writeBINARYDiscontinuous(const std::string &filename, 
                                               const std::vector<ExportEntry>& fields, 
                                               CfemInt entityId) const
    {
        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open()) CFEM_ERROR << "Cannot open VTK file: " << filename;

        auto writeBigEndian = [&](const auto *d, std::size_t count) {
            using T = std::decay_t<decltype(d[0])>;
            for (std::size_t i = 0; i < count; ++i) {
                T val = d[i];
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
                char *p = reinterpret_cast<char *>(&val);
                std::reverse(p, p + sizeof(T));
#endif
                out.write(reinterpret_cast<const char *>(&val), sizeof(T));
            }
        };

        auto geo = buildDiscontinuousGeometry(entityId);
        auto data = evaluateFields(geo, fields, entityId, true);

        out << "# vtk DataFile Version 3.0\nCFEM++ Discontinuous Export\nBINARY\nDATASET UNSTRUCTURED_GRID\n";
        
        // --- POINTS ---
        out << "POINTS " << geo.numActiveVerts << " double\n";
        writeBigEndian(geo.points.data(), geo.points.size());
        out << "\n";

        // --- CELLS ---
        int totalConnSize = geo.connectivity.size() + geo.numActiveCells;
        out << "CELLS " << geo.numActiveCells << " " << totalConnSize << "\n";
        std::vector<int> cellBuffer;
        cellBuffer.reserve(totalConnSize);
        for (size_t i = 0; i < geo.numActiveCells; ++i) {
            int start = (i == 0) ? 0 : geo.offsets[i-1];
            int end = geo.offsets[i];
            cellBuffer.push_back(end - start);
            for (int j = start; j < end; ++j) cellBuffer.push_back(geo.connectivity[j]);
        }
        writeBigEndian(cellBuffer.data(), cellBuffer.size());
        out << "\n";

        out << "CELL_TYPES " << geo.numActiveCells << "\n";
        std::vector<int> typeBuffer(geo.types.begin(), geo.types.end());
        writeBigEndian(typeBuffer.data(), typeBuffer.size());
        out << "\n";

        // --- POINT DATA ---
        if (!fields.empty() || !data.nodeMasks.empty()) {
            out << "POINT_DATA " << geo.numActiveVerts << "\n";
            for (size_t i = 0; i < fields.size(); ++i) {
                int dim = fields[i].expr->getDimension();
                if (dim == 1) {
                    out << "SCALARS " << fields[i].name << " double 1\nLOOKUP_TABLE default\n";
                    writeBigEndian(data.fields[i].data(), data.fields[i].size());
                    out << "\n";
                } else {
                    out << "VECTORS " << fields[i].name << " double\n";
                    std::vector<double> vecBuffer(3 * geo.numActiveVerts, 0.0);
                    for (size_t v = 0; v < geo.numActiveVerts; ++v) {
                        vecBuffer[3*v+0] = data.fields[i][v*dim+0];
                        vecBuffer[3*v+1] = data.fields[i][v*dim+1];
                        if (dim >= 3) vecBuffer[3*v+2] = data.fields[i][v*dim+2];
                    }
                    writeBigEndian(vecBuffer.data(), vecBuffer.size());
                    out << "\n";
                }
            }
            for (size_t i = 0; i < data.nodeMasks.size(); ++i) {
                out << "SCALARS NodeMask_" << m_exportedEntityNames[i] << " float 1\nLOOKUP_TABLE default\n";
                writeBigEndian(data.nodeMasks[i].data(), data.nodeMasks[i].size());
                out << "\n";
            }
        }

        // --- CELL DATA ---
        if (m_exportQuality || !data.cellMasks.empty()) {
            out << "CELL_DATA " << geo.numActiveCells << "\n";
            if (m_exportQuality) {
                out << "SCALARS Mesh_Quality float 1\nLOOKUP_TABLE default\n";
                writeBigEndian(data.cellQuality.data(), data.cellQuality.size());
                out << "\n";
            }
            for (size_t i = 0; i < data.cellMasks.size(); ++i) {
                out << "SCALARS CellMask_" << m_exportedEntityNames[i] << " float 1\nLOOKUP_TABLE default\n";
                writeBigEndian(data.cellMasks[i].data(), data.cellMasks[i].size());
                out << "\n";
            }
        }
    }

} // namespace cfem