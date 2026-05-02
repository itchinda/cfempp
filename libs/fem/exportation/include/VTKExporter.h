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
 * @file VTKExporter.h
 * @brief Legacy VTK (.vtk) mesh exporter for CFEM++.
 */

#ifndef CFEM_VTKEXPORTER_H
#define CFEM_VTKEXPORTER_H

#include <vector>
#include <string>
#include <memory>
#include "VTKFamilyExporter.h"

namespace cfem
{
    class Mesh;

    /**
     * @class VTKExporter
     * @brief Legacy VTK (.vtk) unstructured grid exporter for CFEM++.
     *
     * @ingroup MeshExporter
     * 
     * @details
     * This class exports mesh geometry and evaluations to the classic VTK DataFile Version 3.0 format.
     * While VTU (XML) is the modern standard, the legacy format remains exceptionally useful for 
     * academic purposes, custom parsing scripts, and rapid visual debugging.
     *
     * **Key Features:**
     * - **Architectural Parity:** Maintains strict mathematical parity with `VTUExporter`, including 
     *   Entity Filtering, Discontinuous Topology (DG), and the Double Mapping Strategy for boundaries.
     * - **Human-Readable ASCII:** The ASCII mode produces a very straightforward text structure 
     *   (`POINTS`, `CELLS`, `POINT_DATA`), making it the ultimate tool for debugging nodal values 
     *   or geometric coordinates in a simple text editor.
     * - **Strict Big-Endian Binary:** The binary mode automatically enforces Big-Endian byte order 
     *   during serialization, ensuring strict compliance with the legacy VTK specification regardless 
     *   of the host machine's architecture.
     *
     * # Usage Example
     * ```cpp
     * std::shared_ptr<Mesh> mesh = ...;  // Create or load a mesh
     *
     * VTKExporter exporter(*mesh);
     * exporter.addExpression("Pressure", pressureExpr);
     * exporter.addExpression("Velocity", velocityExpr);
     *
     * exporter.setFormat(VTKExporter::Format::BINARY);
     * exporter.write("cube"); // Produces cube.cont.vtk
     * ```
     * @warning The VTK Legacy format does not natively support advanced parallel/distributed data 
     *          (.pvtu) or multi-block structures. For HPC production runs, use `VTUExporter`.
     * @see VTUExporter
     */
    class VTKExporter : public VTKFamilyExporter
    {
    public:
        /**
         * @brief Constructs a VTKExporter attached to a specific mesh.
         * @param mesh The reference mesh to be exported.
         */
        explicit VTKExporter(const Mesh &mesh);

    private:

        /**
         * @brief Writes continuous fields to an ASCII formatted file.
         * 
         * @details This method iterates over the active elements (filtered by `entityId`).
         * It uses a Global-to-Local mapping to ensure only the necessary vertices are written.
         * When exporting a boundary (facet), it applies a Double Mapping strategy:
         * the physical position and Jacobians are computed from the 3D Parent Cell, while 
         * geometric normals are extracted from the 2D Facet mapping.
         * 
         * @param filename The path and name of the output file.
         * @param contFields A vector of `ExportEntry` containing continuous expressions/fields.
         * @param entityId The ID of the specific mesh entity to export.
         */
        void writeASCIIContinuous(const std::string&              filename, 
                                  const std::vector<ExportEntry>& contFields, 
                                  CfemInt                         entityId) const override;

        /**
         * @brief Writes continuous fields to a Binary formatted file.
         * 
         * @details Implements the same robust mathematical logic (Global-to-Local mapping and 
         * Double Mapping) as `writeASCIIContinuous`, but serializes the data into raw binary 
         * blocks for enhanced performance and reduced file size.
         * - For VTU: Data is written as appended Base64/Raw blocks.
         * - For VTK Legacy: Data is forced into Big-Endian byte order.
         * 
         * @param filename The path and name of the output file.
         * @param contFields A vector of `ExportEntry` containing continuous expressions/fields.
         * @param entityId The ID of the specific mesh entity to export.
         */
        void writeBINARYContinuous(const std::string&              filename, 
                                   const std::vector<ExportEntry>& contFields,
                                   CfemInt                         entityId) const override;

        /**
         * @brief Writes discontinuous (DG) fields to an ASCII formatted file.
         * 
         * @details For Discontinuous Galerkin (DG) fields or expressions that exhibit jumps 
         * across element interfaces, this method duplicates the mesh vertices for each cell.
         * This creates a purely disconnected topology (`totalDgPoints`), allowing the 
         * visualization software to render sharp gradients and interfacial jumps accurately.
         * 
         * @param filename The path and name of the output file.
         * @param discontFields A vector of `ExportEntry` containing discontinuous expressions/fields.
         * @param entityId The ID of the specific mesh entity to export.
         */
        void writeASCIIDiscontinuous(const std::string&              filename, 
                                     const std::vector<ExportEntry>& discontFields, 
                                     CfemInt                         entityId) const override;

        /**
         * @brief Writes discontinuous (DG) fields to a Binary formatted file.
         * 
         * @details Implements the discontinuous topology generation (vertex duplication per cell) 
         * and the Double Mapping strategy, but outputs the arrays as binary blocks. This is the 
         * standard and most efficient method for exporting large-scale DG simulations.
         * 
         * @param filename The path and name of the output file.
         * @param discontFields A vector of `ExportEntry` containing discontinuous expressions/fields.
         * @param entityId The ID of the specific mesh entity to export.
         */
        void writeBINARYDiscontinuous(const std::string&              filename, 
                                      const std::vector<ExportEntry>& discontFields, 
                                      CfemInt                         entityId) const override;
    };

}
#endif