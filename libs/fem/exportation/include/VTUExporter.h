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

/**
 * @file VTUExporter.h
 * @brief VTK XML (.vtu) mesh exporter for CFEM.
 *
 * This class exports mesh geometry, connectivity, and attached
 * expressions (scalar, vector, tensor) to the VTK XML Unstructured
 * Grid format (.vtu). The VTU format is the modern VTK standard and
 * is recommended for visualization with tools such as ParaView.
 *
 * @ingroup MeshExporter
 *
 * # Overview
 * VTUExporter allows:
 * - Exporting vertex coordinates
 * - Exporting cell connectivity, offsets, and VTK cell types
 * - Attaching and exporting symbolic expressions or FE fields
 * - Optional filtering of mesh entities by name
 * - Optional cell quality metrics
 * - Writing either ASCII or binary-encoded data
 *
 * Expressions added through `addExpression()` or `addField()`
 * are evaluated during export and written as point-based data
 * fields in the output file.
 *
 * # Usage Example
 * ```cpp
 * std::shared_ptr<Mesh> mesh = ...;  // Create or load a mesh
 *
 * VTUExporter exporter(*mesh);
 * exporter.addExpression("Pressure", pressureExpr);
 * exporter.addExpression("Velocity", velocityExpr);
 *
 * exporter.setExportedEntities({"xmin", "xmax"});
 * exporter.exportCellQuality(true);
 *
 * exporter.setFormat(VTUExporter::Format::BINARY);
 * exporter.write("cube"); // Produces cube.vtu
 * ```
 *
 * # Output Structure
 * The generated VTU file contains the following sections:
 *
 * - **Points**
 *   Vertex coordinates of the mesh.
 *
 * - **Cells**
 *   Cell connectivity, offsets, and VTK cell type identifiers.
 *
 * - **PointData**
 *   Evaluated expressions and attached FE fields.
 *
 * - **CellData**
 *   Optional cell quality metrics and entity masks.
 *
 * # ASCII vs Binary
 *
 * - **ASCII**
 *   Human-readable text format. Useful for debugging but produces
 *   significantly larger files.
 *
 * - **Binary**
 *   Base64-encoded binary blocks embedded in the XML file. This is
 *   the recommended format for large meshes because it provides
 *   better performance and smaller file sizes.
 *
 * # Notes
 *
 * - Expressions must be prepared before export.
 * - Expressions are evaluated at mesh nodes.
 * - The exporter does not modify the mesh.
 * - Physical entity filters must correspond to valid mesh entities.
 *
 * # References
 * - VTK XML File Format: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 */


#ifndef CFEM_VTUEXPORTER_H
#define CFEM_VTUEXPORTER_H

#include <vector>
#include <string>
#include <memory>
#include "MeshExporter.h"

namespace cfem
{
    class Mash;
    /**
     * @brief Export CFEM meshes and associated expressions to VTK XML (.vtu) format.
     */
    class VTUExporter : public MeshExporter
    {
    public:
        enum class Format
        {
            ASCII,
            BINARY
        };

        VTUExporter(const Mesh &mesh);

        /**
        * @brief Write mesh and attached fields to file.
        *
        * @param filename Output file name.
        */
        void write(const std::string &filename) const override;

        /**
         * @brief Set the export format: either ASCII or binary-encoded data
         * 
         * @param format: The export format
         */
        void setFormat(Format format) {m_format = format;}

    private:
        Format m_format{Format::BINARY};
        void writeASCIIContinuous(const std::string &filename, const std::vector<ExportEntry>& contField) const;
        void writeBINARYContinuous(const std::string &filename, const std::vector<ExportEntry>& contField) const;
        void writeASCIIDiscontinuous(const std::string &filename, const std::vector<ExportEntry>& discontField) const;
        void writeBINARYDiscontinuous(const std::string &filename, const std::vector<ExportEntry>& discontField) const;
    };

}
#endif