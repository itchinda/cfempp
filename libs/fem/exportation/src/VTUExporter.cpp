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

#include "VTUExporter.h"
#include "CfemExpression.h"

#include <algorithm>
#include <iomanip>
#include <fstream>
#include <unordered_map>
#include <cstdint>

namespace cfem
{

    VTUExporter::VTUExporter(const Mesh &mesh) : VTKFamilyExporter(mesh)
    { 
        m_outputFileExtension = "vtu";
    }

    // ========================================================================
    // --- FORMATTER 1: ASCII CONTINUOUS (I/O ONLY)
    // ========================================================================
    void VTUExporter::writeASCIIContinuous(const std::string &filename, const std::vector<ExportEntry>& fields, CfemInt entityId) const
    {
        std::ofstream out(filename);
        if (!out.is_open()) CFEM_ERROR << "Cannot open VTU file: " << filename;

        auto geo = buildContinuousGeometry(entityId);
        auto data = evaluateFields(geo, fields, entityId, false);

        out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n";
        out << "<Piece NumberOfPoints=\"" << geo.numActiveVerts << "\" NumberOfCells=\"" << geo.numActiveCells << "\">\n";

        out << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (size_t i = 0; i < geo.numActiveVerts; ++i) {
            out << std::setprecision(16) << geo.points[3*i] << " " << geo.points[3*i+1] << " " << geo.points[3*i+2] << "\n";
        }
        out << "</DataArray>\n</Points>\n";

        out << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (CfemInt c : geo.connectivity) out << c << " ";
        out << "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        for (CfemInt o : geo.offsets) out << o << "\n";
        out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (auto t : geo.types) out << static_cast<int>(t) << "\n";
        out << "</DataArray>\n</Cells>\n";

        if (!fields.empty() || !data.nodeMasks.empty()) {
            out << "<PointData>\n";
            for (size_t i = 0; i < fields.size(); ++i) {
                CfemInt dim = fields[i].expr->getNumComponents();
                out << "<DataArray type=\"Float64\" Name=\"" << fields[i].name << "\" NumberOfComponents=\"" << dim << "\" format=\"ascii\">\n";
                for (size_t v = 0; v < geo.numActiveVerts; ++v) {
                    for (CfemInt d = 0; d < dim; ++d) out << data.fields[i][v * dim + d] << " ";
                    out << "\n";
                }
                out << "</DataArray>\n";
            }
            for (size_t i = 0; i < data.nodeMasks.size(); ++i) {
                out << "<DataArray type=\"Float32\" Name=\"NodeMask_" << m_exportedEntityNames[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
                for (size_t v = 0; v < geo.numActiveVerts; ++v) out << data.nodeMasks[i][v] << "\n";
                out << "</DataArray>\n";
            }
            out << "</PointData>\n";
        }
        
        if (m_exportQuality || !data.cellMasks.empty()) {
            out << "<CellData>\n";
            if (m_exportQuality) {
                out << "<DataArray type=\"Float32\" Name=\"Mesh_Quality\" NumberOfComponents=\"1\" format=\"ascii\">\n";
                for (auto q : data.cellQuality) out << q << "\n";
                out << "</DataArray>\n";
            }
            for (size_t i = 0; i < data.cellMasks.size(); ++i) {
                out << "<DataArray type=\"Float32\" Name=\"CellMask_" << m_exportedEntityNames[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
                for (auto m : data.cellMasks[i]) out << m << "\n";
                out << "</DataArray>\n";
            }
            out << "</CellData>\n";
        }
        
        out << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    }

    // ========================================================================
    // --- FORMATTER 2: BINARY CONTINUOUS (I/O ONLY)
    // ========================================================================
    void VTUExporter::writeBINARYContinuous(const std::string &filename, const std::vector<ExportEntry>& fields, CfemInt entityId) const
    {
        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open()) CFEM_ERROR << "Cannot open VTU file: " << filename;

        auto geo = buildContinuousGeometry(entityId);
        auto data = evaluateFields(geo, fields, entityId, false);

        std::vector<char> appended;
        auto appendBlock = [&](const auto &vec) -> size_t {
            if (vec.empty()) return 0;
            size_t offset = appended.size();
            uint32_t bytes = static_cast<uint32_t>(vec.size() * sizeof(typename std::decay_t<decltype(vec)>::value_type));
            appended.insert(appended.end(), reinterpret_cast<const char *>(&bytes), reinterpret_cast<const char *>(&bytes) + sizeof(uint32_t));
            appended.insert(appended.end(), reinterpret_cast<const char *>(vec.data()), reinterpret_cast<const char *>(vec.data()) + bytes);
            return offset;
        };

        size_t ptsOff = appendBlock(geo.points);
        size_t connOff = appendBlock(geo.connectivity);
        size_t offOff = appendBlock(geo.offsets);
        size_t typeOff = appendBlock(geo.types);
        
        std::vector<size_t> fieldOffs;
        for (const auto& f : data.fields) fieldOffs.push_back(appendBlock(f));
        std::vector<size_t> nodeMaskOffs;
        for (const auto& m : data.nodeMasks) nodeMaskOffs.push_back(appendBlock(m));
        
        size_t qOff = m_exportQuality ? appendBlock(data.cellQuality) : 0;
        std::vector<size_t> cellMaskOffs;
        for (const auto& m : data.cellMasks) cellMaskOffs.push_back(appendBlock(m));

        out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n";
        out << "<Piece NumberOfPoints=\"" << geo.numActiveVerts << "\" NumberOfCells=\"" << geo.numActiveCells << "\">\n";
        
        out << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << ptsOff << "\"/>\n</Points>\n";
        out << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"" << connOff << "\"/>\n";
        out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"" << offOff << "\"/>\n";
        out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << typeOff << "\"/>\n</Cells>\n";

        if (!fields.empty() || !data.nodeMasks.empty()) {
            out << "<PointData>\n";
            for (size_t i = 0; i < fields.size(); ++i) {
                out << "<DataArray type=\"Float64\" Name=\"" << fields[i].name 
                    << "\" NumberOfComponents=\"" << fields[i].expr->getNumComponents() 
                    << "\" format=\"appended\" offset=\"" << fieldOffs[i] << "\"/>\n";
            }
            for (size_t i = 0; i < data.nodeMasks.size(); ++i) {
                out << "<DataArray type=\"Float32\" Name=\"NodeMask_" << m_exportedEntityNames[i] 
                    << "\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << nodeMaskOffs[i] 
                    << "\"/>\n";
            }
            out << "</PointData>\n";
        }
        if (m_exportQuality || !data.cellMasks.empty()) {
            out << "<CellData>\n";
            if (m_exportQuality) {
                out << "<DataArray type=\"Float32\" Name=\"Mesh_Quality\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << qOff << "\"/>\n";
            }
            for (size_t i = 0; i < data.cellMasks.size(); ++i) {
                out << "<DataArray type=\"Float32\" Name=\"CellMask_" << m_exportedEntityNames[i] 
                    << "\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << cellMaskOffs[i] 
                    << "\"/>\n";
            }
            out << "</CellData>\n";
        }

        out << "</Piece>\n</UnstructuredGrid>\n<AppendedData encoding=\"raw\">\n_";
        out.write(appended.data(), appended.size());
        out << "\n</AppendedData>\n</VTKFile>\n";
    }

    // ========================================================================
    // --- FORMATTER 3: ASCII DISCONTINUOUS (I/O ONLY)
    // ========================================================================
    void VTUExporter::writeASCIIDiscontinuous(const std::string &filename, const std::vector<ExportEntry>& fields, CfemInt entityId) const
    {
        std::ofstream out(filename);
        if (!out.is_open()) CFEM_ERROR << "Cannot open VTU file: " << filename;

        auto geo = buildDiscontinuousGeometry(entityId);
        auto data = evaluateFields(geo, fields, entityId, true);

        out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n";
        out << "<Piece NumberOfPoints=\"" << geo.numActiveVerts << "\" NumberOfCells=\"" << geo.numActiveCells << "\">\n";

        out << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (size_t i = 0; i < geo.numActiveVerts; ++i) {
            out << std::setprecision(16) << geo.points[3*i] << " " << geo.points[3*i+1] << " " << geo.points[3*i+2] << "\n";
        }
        out << "</DataArray>\n</Points>\n";

        out << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (CfemInt c : geo.connectivity) out << c << " ";
        out << "\n</DataArray>\n<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        for (CfemInt o : geo.offsets) out << o << "\n";
        out << "</DataArray>\n<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (auto t : geo.types) out << static_cast<int>(t) << "\n";
        out << "</DataArray>\n</Cells>\n";

        // --- POINT DATA (Fields & Node Masks) ---
        if (!fields.empty() || !data.nodeMasks.empty()) {
            out << "<PointData>\n";
            for (size_t i = 0; i < fields.size(); ++i) {
                CfemInt dim = fields[i].expr->getNumComponents();
                out << "<DataArray type=\"Float64\" Name=\"" << fields[i].name << "\" NumberOfComponents=\"" << dim << "\" format=\"ascii\">\n";
                for (size_t v = 0; v < geo.numActiveVerts; ++v) {
                    for (CfemInt d = 0; d < dim; ++d) out << data.fields[i][v * dim + d] << " ";
                    out << "\n";
                }
                out << "</DataArray>\n";
            }
            for (size_t i = 0; i < data.nodeMasks.size(); ++i) {
                out << "<DataArray type=\"Float32\" Name=\"NodeMask_" << m_exportedEntityNames[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
                for (size_t v = 0; v < geo.numActiveVerts; ++v) out << data.nodeMasks[i][v] << "\n";
                out << "</DataArray>\n";
            }
            out << "</PointData>\n";
        }
        
        // --- CELL DATA (Quality & Cell Masks) ---
        if (m_exportQuality || !data.cellMasks.empty()) {
            out << "<CellData>\n";
            if (m_exportQuality) {
                out << "<DataArray type=\"Float32\" Name=\"Mesh_Quality\" NumberOfComponents=\"1\" format=\"ascii\">\n";
                for (auto q : data.cellQuality) out << q << "\n";
                out << "</DataArray>\n";
            }
            for (size_t i = 0; i < data.cellMasks.size(); ++i) {
                out << "<DataArray type=\"Float32\" Name=\"CellMask_" << m_exportedEntityNames[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
                for (auto m : data.cellMasks[i]) out << m << "\n";
                out << "</DataArray>\n";
            }
            out << "</CellData>\n";
        }
        
        out << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    }

    // ========================================================================
    // FORMATTER 4: BINARY DISCONTINUOUS (I/O ONLY)
    // ========================================================================
    void VTUExporter::writeBINARYDiscontinuous(const std::string &filename, const std::vector<ExportEntry>& fields, CfemInt entityId) const
    {
        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open()) CFEM_ERROR << "Cannot open VTU file: " << filename;

        auto geo = buildDiscontinuousGeometry(entityId);
        auto data = evaluateFields(geo, fields, entityId, true);

        std::vector<char> appended;
        auto appendBlock = [&](const auto &vec) -> size_t {
            if (vec.empty()) return 0;
            size_t offset = appended.size();
            uint32_t bytes = static_cast<uint32_t>(vec.size() * sizeof(typename std::decay_t<decltype(vec)>::value_type));
            appended.insert(appended.end(), reinterpret_cast<const char *>(&bytes), reinterpret_cast<const char *>(&bytes) + sizeof(uint32_t));
            appended.insert(appended.end(), reinterpret_cast<const char *>(vec.data()), reinterpret_cast<const char *>(vec.data()) + bytes);
            return offset;
        };

        size_t ptsOff = appendBlock(geo.points);
        size_t connOff = appendBlock(geo.connectivity);
        size_t offOff = appendBlock(geo.offsets);
        size_t typeOff = appendBlock(geo.types);
        
        // --- OFFSETS COMPUTATION ---
        std::vector<size_t> fieldOffs;
        for (const auto& f : data.fields) fieldOffs.push_back(appendBlock(f));
        std::vector<size_t> nodeMaskOffs;
        for (const auto& m : data.nodeMasks) nodeMaskOffs.push_back(appendBlock(m));
        
        size_t qOff = m_exportQuality ? appendBlock(data.cellQuality) : 0;
        std::vector<size_t> cellMaskOffs;
        for (const auto& m : data.cellMasks) cellMaskOffs.push_back(appendBlock(m));

        // --- XML HEADER & GEOMETRY ---
        out << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n<UnstructuredGrid>\n";
        out << "<Piece NumberOfPoints=\"" << geo.numActiveVerts << "\" NumberOfCells=\"" << geo.numActiveCells << "\">\n";
        
        out << "<Points>\n<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << ptsOff << "\"/>\n</Points>\n";
        out << "<Cells>\n<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"" << connOff << "\"/>\n";
        out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"" << offOff << "\"/>\n";
        out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << typeOff << "\"/>\n</Cells>\n";

        // --- POINT DATA ---
        if (!fields.empty() || !data.nodeMasks.empty()) {
            out << "<PointData>\n";
            for (size_t i = 0; i < fields.size(); ++i){
                out << "<DataArray type=\"Float64\" Name=\"" << fields[i].name 
                    << "\" NumberOfComponents=\"" << fields[i].expr->getNumComponents() 
                    << "\" format=\"appended\" offset=\"" << fieldOffs[i] 
                    << "\"/>\n";
            }
            for (size_t i = 0; i < data.nodeMasks.size(); ++i) {
                out << "<DataArray type=\"Float32\" Name=\"NodeMask_" << m_exportedEntityNames[i] 
                    << "\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << nodeMaskOffs[i] 
                    << "\"/>\n";
            }
            out << "</PointData>\n";
        }
        
        // --- CELL DATA ---
        if (m_exportQuality || !data.cellMasks.empty()) {
            out << "<CellData>\n";
            if (m_exportQuality) {
                out << "<DataArray type=\"Float32\" Name=\"Mesh_Quality\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << qOff
                    << "\"/>\n";
            }
            for (size_t i = 0; i < data.cellMasks.size(); ++i) {
                out << "<DataArray type=\"Float32\" Name=\"CellMask_" << m_exportedEntityNames[i] 
                    << "\" NumberOfComponents=\"1\" format=\"appended\" offset=\"" << cellMaskOffs[i] 
                    << "\"/>\n";
                }
            out << "</CellData>\n";
        }

        // --- APPENDED RAW DATA ---
        out << "</Piece>\n</UnstructuredGrid>\n<AppendedData encoding=\"raw\">\n_";
        out.write(appended.data(), appended.size());
        out << "\n</AppendedData>\n</VTKFile>\n";
    }

} // namespace cfem