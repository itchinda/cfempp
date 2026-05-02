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
 * @brief Modern VTK XML (.vtu) unstructured grid exporter for CFEM++.
 */

#ifndef CFEM_VTUEXPORTER_H
#define CFEM_VTUEXPORTER_H

#include <vector>
#include <string>
#include <memory>
#include "VTKFamilyExporter.h"

namespace cfem
{
    class Mesh;
    /**
     * @class VTUExporter
     * @brief Modern VTK XML (.vtu) unstructured grid exporter for CFEM++.
     * 
     * @ingroup MeshExporter
     *
     * @details
     * This class exports mesh geometry, connectivity, and attached symbolic expressions
     * to the VTK XML Unstructured Grid format (`.vtu`). It represents the modern standard
     * for Paraview visualization, supporting highly optimized binary data blocks.
     *
     * **Key Features:**
     * - **Entity Filtering:** Capable of exporting the entire universe or isolating specific 
     *   sub-domains and boundary manifolds (e.g., extracting a single 2D boundary from a 3D mesh).
     * - **Continuous & Discontinuous Galerkin (DG):** Automatically detects field continuity. 
     *   For DG fields, it duplicates nodes to allow exact representation of jumps at cell interfaces.
     * - **Double Mapping Strategy:** When evaluating fields on boundaries (facets), it seamlessly 
     *   combines the 3D Parent Cell mapping (for exact physical gradients via the inverse Jacobian) 
     *   with the 2D Facet mapping (for exact geometric normals).
     * - **Appended Raw Data:** In binary mode, uses the `<AppendedData>` section to write raw, 
     *   contiguous memory blocks, drastically reducing file size and loading times in ParaView.
     *
     * # Usage Example
     * ```cpp
     * std::shared_ptr<Mesh> mesh = ...;  // Create or load a mesh
     *
     * VTUExporter exporter(*mesh);
     * exporter.addExpression("Pressure", pressureExpr);
     * exporter.addExpression("Velocity", velocityExpr);
     *
     * exporter.setFormat(VTUExporter::Format::BINARY);
     * exporter.write("cube"); // Produces cube.cont.vtu
     * ```
     * @note For huge production simulations, prefer `Format::BINARY`.
     * @see VTKExporter for the legacy ASCII/Binary format (useful for text-editor debugging).
     * 
     * @ref
     * - VTK XML File Format: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
     */
    class VTUExporter : public VTKFamilyExporter
    {
    public:
        /**
         * @brief Constructs a VTUExporter attached to a specific mesh.
         * @param mesh The reference mesh to be exported.
         */
        VTUExporter(const Mesh &mesh);

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
        void writeASCIIContinuous(const std::string &filename, 
                                const std::vector<ExportEntry>& contFields, 
                                CfemInt entityId) const override;

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
        void writeBINARYContinuous(const std::string &filename, 
                                const std::vector<ExportEntry>& contFields, 
                                CfemInt entityId) const override;

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
        void writeASCIIDiscontinuous(const std::string &filename, 
                                    const std::vector<ExportEntry>& discontFields, 
                                    CfemInt entityId) const override;

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
        void writeBINARYDiscontinuous(const std::string &filename, 
                                    const std::vector<ExportEntry>& discontFields, 
                                    CfemInt entityId) const override;
        };

}
#endif