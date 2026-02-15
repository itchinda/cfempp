#ifndef CFEM_MESH_EXPORTER_H
#define CFEM_MESH_EXPORTER_H

#include "cfemutils.h"
#include "FEField.h"
#include "CfemExpression.h"
#include "cfem_topology.h"

#include <memory>
#include <string>

namespace cfem {

    class Mesh;
    //class CfemExpression;

    /**
     * @defgroup MeshExporter Mesh Export Infrastructure
     * @brief Classes responsible for exporting meshes and associated fields.
     *
     * This module provides the abstraction layer used to export CFEM meshes
     * and associated scalar fields to external visualization or data formats
     * such as VTK, VTU, XDMF(planned), etc.
     *
     * Exporters operate on a given Mesh instance and allow attaching
     * expressions or finite element fields that are evaluated during output.
     *
     * @{
     */

    /**
     * @brief Abstract base class for mesh exporters.
     *
     * Provides a polymorphic interface for exporting mesh geometry
     * and associated scalar fields to various file formats.
     *
     * Derived classes implement the write() function to support
     * specific output formats.
     *
     * The exporter does not modify the mesh. It only reads geometry,
     * connectivity, and evaluates attached expressions.
     */
    class MeshExporter
    {
    public:

        /**
         * @brief Construct an exporter for a given mesh.
         *
         * @param mesh Reference to the mesh to export.
         */
        explicit MeshExporter(const Mesh& mesh)
            : m_mesh(mesh)
        {}

        /// Virtual destructor.
        virtual ~MeshExporter() = default;

        /**
         * @brief Attach an expression to be exported.
         *
         * The expression is evaluated at nodes or cells depending
         * on the exporter implementation.
         *
         * The exporter takes shared ownership of the expression.
         *
         * @param name Name under which the field will appear in the output file.
         * @param expr Shared pointer to the expression to export.
         *
         * @throws Runtime error if a field with the same name already exists.
         */
        void addExpression(const std::string& name,
                           std::shared_ptr<const sym::CfemExpression> expr)
        {
            for (const auto& [name_used, stored_expr] : m_entries) {
                if (name_used == name) {
                    CFEM_ERROR("A field with name '" + name + "' is already present.");
                }
            }
            m_entries.push_back({name, std::move(expr)});
        }

        /**
         * @brief Attach a finite element field for export.
         *
         * Convenience wrapper that forwards to addExpression().
         *
         * @param field Shared pointer to the FEField.
         */
        void addField(std::shared_ptr<const FEField> field)
        {
            // addExpression(field->getName(), std::move(field));
            addExpression(field->getName(), field);
        }

        /**
         * @brief Restrict export to specific physical entities.
         *
         * If empty, all entities are exported.
         *
         * @param entityNames List of physical entity names to export.
         */
        void setExportedEntities(const std::vector<std::string>& entityNames)
        {
            m_exportedEntityNames = entityNames;
        }

        /**
         * @brief Enable or disable cell quality export.
         *
         * When enabled, the exporter may compute and output
         * element quality metrics.
         *
         * @param enable True to export quality metrics.
         */
        void exportCellQuality(bool enable)
        {
            m_exportQuality = enable;
        }

        /**
         * @brief Write mesh and attached fields to file.
         *
         * Must be implemented by derived classes.
         *
         * @param filename Output file name.
         */
        virtual void write(const std::string& filename) const = 0;

    protected:

        /**
         * @brief Internal representation of an exported field.
         */
        struct ExportEntry {
            std::string name;                                   ///< Field name
            std::shared_ptr<const sym::CfemExpression> expr;         ///< Expression to evaluate
        };

        const Mesh& m_mesh;                                     ///< Mesh being exported
        bool m_exportQuality {false};                           ///< Flag for quality export
        std::vector<ExportEntry> m_entries;                     ///< Attached fields
        std::vector<std::string> m_exportedEntityNames;         ///< Filtered entities

        /**
         * @brief Maps internal CellType to VTK legacy cell identifiers.
         * @param type The internal CFEM CellType.
         * @return The corresponding VTK integer ID.
         */
        static int getVTKCellType(CellType type)
        {
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

        std::string validateFilename(std::string filename,
                                     const std::string& requiredExtension) const;
    
    };

    /** @} */

} // namespace cfem

#endif