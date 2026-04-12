#include "MeshExporter.h"
#include <regex>
#include <filesystem>
#include "Mesh.h"
#include "FEField.h"
#include "CfemExpression.h"

namespace fs = std::filesystem;

namespace cfem
{
    MeshExporter::MeshExporter(const Mesh& mesh)
            : m_mesh(mesh)
        {}

    void MeshExporter::addExpression(const std::string& name,
                           std::shared_ptr<const sym::CfemExpression> expr)
    {
        for (const auto& [name_used, stored_expr] : m_entries) {
            if (name_used == name) {
                CFEM_ERROR << "A field with name '" << name << "' is already present.";
            }
        }
        m_entries.push_back({name, std::move(expr)});
    }

    void MeshExporter::addField(std::shared_ptr<const FEField> field)
    {
        addExpression(field->getName(), field);
    }

    void MeshExporter::setExportedEntities(const std::vector<std::string>& entityNames)
    {
        m_exportedEntityNames = entityNames;
    }

    void MeshExporter::exportCellQuality(bool enable)
    {
        m_exportQuality = enable;
    }

    std::string MeshExporter::validateFilename(std::string filename,
                                               const std::string &requiredExtension) const
    {
        if (filename.empty())
            CFEM_ERROR << "Filename cannot be empty.";

        // Remove trailing spaces
        while (!filename.empty() && std::isspace(filename.back()))
            filename.pop_back();

        if (filename.empty())
            CFEM_ERROR << "Filename cannot be empty after trimming.";

        // Reject dangerous characters
        static const std::regex allowedPattern(
            R"(^[a-zA-Z0-9_\-./\\]+$)");

        if (!std::regex_match(filename, allowedPattern))
            CFEM_ERROR << "Filename contains invalid characters: " << filename;

        // Normalize extension
        std::filesystem::path p(filename);

        if (!p.has_extension())
        {
            p.replace_extension(requiredExtension);
        }
        else if (p.extension() != requiredExtension)
        {
            CFEM_ERROR << "Invalid file extension. Expected " << requiredExtension;
        }

        return p.string();
    }

    CfemInt MeshExporter::getVTKCellType(CellType type)
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
            default: 
                CFEM_ERROR << "Unsupported cell type."; // VTK_EMPTY_CELL case 13: return CellType::PRISM6; case 32: return CellType::PRISM18;
                return 0;
        }
    }
}
