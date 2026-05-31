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
 * @file VTKFamilyExporter.h
 * @brief Intermediate base class for all unstructured grid exporters.
 */
#ifndef CFEM_VTKFAMILYEXPORTER_H
#define CFEM_VTKFAMILYEXPORTER_H

#include "MeshExporter.h"
#include <vector>
#include <cstdint>
#include <string>

namespace cfem
{
    /**
     * @brief Output encoding formats for VTK-family exporters.
     */
    enum class VTKExportFormat
    {
        ASCII,  ///< Plain text format. The best choice for visual inspection and debugging.
        BINARY  ///< Big-Endian or Base64 binary format. Smaller size and faster I/O.
    };

    /**
     * @class VTKFamilyExporter
     * @brief Intermediate base class for all unstructured grid exporters.
     * * @ingroup MeshExporter
     * * @details
     * This class centralizes all the heavy mathematical and topological logic required 
     * for exporting unstructured meshes, including:
     * - **Entity Filtering:** Extraction of specific sub-domains or boundary manifolds.
     * - **Discontinuous Topology:** Generation of Discontinuous Galerkin (DG) topology via node duplication.
     * - **Field Evaluation:** Evaluating expressions using the "Double Mapping" strategy (Volume + Facet) to properly handle exact physical gradients and geometric normals.
     * * Derived classes (such as `VTUExporter` and `VTKExporter`) only need to implement 
     * the low-level formatting and I/O logic (ASCII, Binary, XML, etc.) by overriding 
     * the specialized protected `write...` methods.
     */
    class VTKFamilyExporter : public MeshExporter
    {
    public:

        virtual ~VTKFamilyExporter() = default;

        using MeshExporter::write;

        /**
         * @brief Evaluates all attached expressions and writes the mesh file.
         * * @param filename Base name of the output file. Extensions (e.g., `.cont.vtu`, `.disc.vtk`) 
         * are automatically appended based on field continuity.
         * @param entityId Optional ID of the specific mesh entity to export. Defaults to the 
         * entire universe (`UNIVERSE_ID`).
         * @throws CFEM_ERROR If the file cannot be created or opened.
         */
        void write(const std::string &filename, CfemInt entityId = UNIVERSE_ID) const override;

        /**
         * @brief Sets the data encoding format.
         * @param format The desired export format (ASCII or BINARY).
         */
        void setFormat(VTKExportFormat format) { m_format = format; }

    protected:
        
        std::string m_outputFileExtension {""};     ///< File extension defined by derived classes (e.g., "vtu", "vtk").
        VTKExportFormat m_format{VTKExportFormat::BINARY};  ///< Active export format.

        /**
         * @brief Constructs a VTKFamilyExporter tied to a specific mesh.
         * @param mesh The CFEM++ mesh object to be exported.
         */
        explicit VTKFamilyExporter(const Mesh &mesh) : MeshExporter(mesh) {}

        /**
         * @brief Internal geometry representation decoupled from the CFEM++ core mesh.
         * @details This structure holds the flattened arrays strictly required by VTK/VTU formats,
         * filtering out inactive elements if a specific sub-domain is being exported.
         */
        struct ExtractedGeometry {
            CfemInt numActiveCells = 0;       ///< Number of cells actually exported.
            CfemInt numActiveVerts = 0;       ///< Number of unique vertices actually exported.
            
            std::vector<double> points;       ///< Flattened 3D coordinates [x0, y0, z0, x1, y1, z1, ...].
            std::vector<int> connectivity;    ///< Flattened node indices defining each cell's connectivity.
            std::vector<int> offsets;         ///< Array indicating the end index of each cell within the connectivity array.
            std::vector<unsigned char> types; ///< VTK cell type identifiers (e.g., VTK_TRIANGLE, VTK_TETRA).
            
            std::vector<CfemInt> globalToLocal; ///< Mapping from CFEM++ global node IDs to local export node IDs.
            std::vector<CfemInt> localToGlobal; ///< Mapping from local export node IDs to CFEM++ global node IDs.
            
            const std::vector<CfemInt>* pActiveElements = nullptr; ///< Pointer to the list of active cell IDs in the mesh.
        };

        /**
         * @brief Container for evaluated nodal fields and element metadata.
         */
        struct EvaluatedData {
            std::vector<std::vector<CfemReal>> fields;    ///< Evaluated scalar/vector field values. Outer size: number of fields. Inner size: number of nodes.
            std::vector<std::vector<float>> nodeMasks;    ///< Node validity masks (e.g., for thresholding or clipping in ParaView).
            std::vector<float> cellQuality;               ///< Quality metrics computed per cell.
            std::vector<std::vector<float>> cellMasks;    ///< Cell validity masks.
        };

        // --- Core Mathematical Kernels ---
        
        /**
         * @brief Builds the topological and geometric data for a standard continuous export.
         * @details Nodes shared between adjacent cells are written once, maintaining standard connectivity.
         * @param entityId The ID of the mesh entity to extract.
         * @return The populated ExtractedGeometry structure.
         */
        ExtractedGeometry buildContinuousGeometry(CfemInt entityId) const;

        /**
         * @brief Builds the topological and geometric data for a discontinuous (DG) export.
         * @details Nodes shared between adjacent cells are intentionally duplicated. Each cell receives 
         * its own unique set of nodes. This prevents ParaView/VisIt from applying nodal interpolation 
         * smoothing across material interfaces or elements evaluating discontinuous gradients.
         * @param entityId The ID of the mesh entity to extract.
         * @return The populated ExtractedGeometry structure.
         */
        ExtractedGeometry buildDiscontinuousGeometry(CfemInt entityId) const;
        
        /**
         * @brief Evaluates all requested fields (expressions) on the extracted geometry.
         * @details Utilizes the CFEM++ functional evaluator to compute field values at nodes. 
         * It intelligently handles both standard volume evaluations and boundary evaluations 
         * (which require projecting facet coordinates into the parent volume).
         * @param geo The previously extracted geometry to evaluate fields on.
         * @param entries The list of expressions/fields to evaluate.
         * @param entityId The ID of the targeted mesh entity.
         * @param isDiscontinuous Flag indicating if the geometry is DG (affects node-to-cell mapping logic).
         * @return An EvaluatedData structure containing the flattened results ready for VTK serialization.
         */
        EvaluatedData evaluateFields(const ExtractedGeometry&        geo, 
                                     const std::vector<ExportEntry>& entries, 
                                     CfemInt                         entityId, 
                                     bool                            isDiscontinuous) const;

        // --- Abstract Formatters to be implemented by derived classes ---

        /**
         * @brief Writes continuous fields to an ASCII formatted file.
         * @param filename The path and name of the output file.
         * @param contFields A vector of `ExportEntry` containing continuous expressions/fields.
         * @param entityId The ID of the specific mesh entity to export.
         */
        virtual void writeASCIIContinuous(const std::string&              filename, 
                                          const std::vector<ExportEntry>& contFields, 
                                          CfemInt                         entityId) const {};

        /**
         * @brief Writes continuous fields to a Binary formatted file.
         * @param filename The path and name of the output file.
         * @param contFields A vector of `ExportEntry` containing continuous expressions/fields.
         * @param entityId The ID of the specific mesh entity to export.
         */
        virtual void writeBINARYContinuous(const std::string&              filename, 
                                           const std::vector<ExportEntry>& contFields,
                                           CfemInt                         entityId) const {};

        /**
         * @brief Writes discontinuous (DG) fields to an ASCII formatted file.
         * @param filename The path and name of the output file.
         * @param discontFields A vector of `ExportEntry` containing discontinuous expressions/fields.
         * @param entityId The ID of the specific mesh entity to export.
         */
        virtual void writeASCIIDiscontinuous(const std::string&              filename, 
                                             const std::vector<ExportEntry>& discontFields, 
                                             CfemInt                         entityId) const {};

        /**
         * @brief Writes discontinuous (DG) fields to a Binary formatted file.
         * @param filename The path and name of the output file.
         * @param discontFields A vector of `ExportEntry` containing discontinuous expressions/fields.
         * @param entityId The ID of the specific mesh entity to export.
         */
        virtual void writeBINARYDiscontinuous(const std::string&              filename, 
                                              const std::vector<ExportEntry>& discontFields, 
                                              CfemInt                         entityId) const {};
        
    };

} // namespace cfem
#endif