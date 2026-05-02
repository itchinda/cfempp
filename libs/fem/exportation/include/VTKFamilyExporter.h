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
     * 
     * @ingroup MeshExporter
     * 
     * @details
     * This class centralizes all the heavy mathematical and topological logic required 
     * for exporting unstructured meshes, including:
     * - **Entity Filtering:** Extraction of specific sub-domains or boundary manifolds.
     * - **Discontinuous Topology:** Generation of Discontinuous Galerkin (DG) topology via node duplication.
     * - **Field Evaluation:** Evaluating expressions using the "Double Mapping" strategy (Volume + Facet) to properly handle exact physical gradients and geometric normals.
     * 
     * Derived classes (such as `VTUExporter` and `VTKExporter`) only need to implement 
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
         * 
         * @param filename Base name of the output file. Extensions (e.g., `.cont.vtu`, `.disc.vtk`) 
         *                 are automatically appended based on field continuity.
         * @param entityId Optional ID of the specific mesh entity to export. Defaults to the 
         *                 entire universe (`UNIVERSE_ID`).
         * @throws CFEM_ERROR If the file cannot be created or opened.
         */
        void write(const std::string &filename, CfemInt entityId = UNIVERSE_ID) const override;

        /**
         * @brief Sets the data encoding format.
         * @param format The desired export format (ASCII or BINARY).
         */
        void setFormat(VTKExportFormat format) { m_format = format; }

    protected:
        
        std::string m_outputFileExtension {""};             ///< File extension defined by derived classes (e.g., "vtu", "vtk")
        VTKExportFormat m_format{VTKExportFormat::BINARY};  ///< Active export format

        explicit VTKFamilyExporter(const Mesh &mesh) : MeshExporter(mesh) {}

        /**
         * @brief Internal geometry representation decoupled from the CFEM++ core mesh.
         */
        struct ExtractedGeometry {
            CfemInt numActiveCells = 0;
            CfemInt numActiveVerts = 0;
            std::vector<double> points;
            std::vector<int> connectivity;
            std::vector<int> offsets;
            std::vector<unsigned char> types;
            
            std::vector<CfemInt> globalToLocal; 
            std::vector<CfemInt> localToGlobal;
            
            const std::vector<CfemInt>* pActiveElements = nullptr;
        };

        /**
         * @brief Container for evaluated nodal fields and element metadata.
         */
        struct EvaluatedData {
            std::vector<std::vector<CfemReal>> fields;
            std::vector<std::vector<float>> nodeMasks;
            std::vector<float> cellQuality;
            std::vector<std::vector<float>> cellMasks;
        };

        // --- Core Mathematical Kernels ---
        
        ExtractedGeometry buildContinuousGeometry(CfemInt entityId) const;
        ExtractedGeometry buildDiscontinuousGeometry(CfemInt entityId) const;
        
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