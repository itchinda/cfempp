/**
 * @file VTKExporter.h
 * @brief Legacy VTK (.vtk) mesh exporter for CFEM.
 * 
 * This class exports mesh geometry, connectivity, and attached
 * expressions (scalar, vector, tensor) to the legacy VTK file format.
 * It supports both ASCII and binary output.
 * 
 * @ingroup MeshExporter
 *
 * # Overview
 * VTKExporter allows:
 * - Exporting vertex coordinates
 * - Exporting cell connectivity and VTK cell types
 * - Attaching and exporting expressions defined on nodes
 * - Optional point/cell masks based on named mesh entities
 * - Optional cell quality metrics
 *
 * # Usage Example
 * ```cpp
 * std::shared_ptr<Mesh> mesh = ...;  // Create or load a mesh
 * VTKExporter exporter(*mesh);
 * exporter.addExpression("Pressure", pressureExpr);
 * exporter.addExpression("Velocity", velocityExpr);
 * exporter.setExportedEntities({"xmin", "xmax"});
 * exporter.exportCellQuality(true);
 * exporter.write("cube.vtk");  // Chooses ASCII or binary internally
 * ```
 *
 * # Binary vs ASCII
 * - **ASCII**: Human-readable, larger file size, slower to write.
 * - **Binary**: Raw binary blocks, faster to write, smaller file.
 *
 * # Data Layout
 * - `POINTS`: 3D vertex coordinates
 * - `CELLS`: connectivity, offsets, and VTK cell types
 * - `POINT_DATA`: expressions and optional point masks
 * - `CELL_DATA`: optional cell quality and cell masks
 *
 * # Notes
 * - Expressions must be prepared via `CfemExpression::prepare()`.
 * - Only node-based expressions are supported in legacy VTK.
 * - Named entities must exist in the mesh for masks to be generated.
 * - Binary files use Little Endian encoding.
 */


#ifndef CFEM_VTKEXPORTER_H
#define CFEM_VTKEXPORTER_H

#include <vector>
#include <string>
#include <memory>
#include "MeshExporter.h"

namespace cfem
{
    class Mesh;

    /**
     * @brief Export CFEM meshes and associated expressions to VTK  (.vtk) format.
     */
    class VTKExporter : public MeshExporter
    {
    public:
        enum class Format
        {
            ASCII,
            BINARY
        };

        VTKExporter(const Mesh &mesh);

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
        void writeASCII(const std::string &filename) const;
        void writeBINARY(const std::string &filename) const;
    };

}
#endif