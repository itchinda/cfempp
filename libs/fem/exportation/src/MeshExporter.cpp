#include "MeshExporter.h"
#include <regex>
#include <filesystem>

namespace fs = std::filesystem;

namespace cfem
{
    std::string MeshExporter::validateFilename(std::string filename,
                                               const std::string &requiredExtension) const
    {
        if (filename.empty())
            CFEM_ERROR("Filename cannot be empty.");

        // Remove trailing spaces
        while (!filename.empty() && std::isspace(filename.back()))
            filename.pop_back();

        if (filename.empty())
            CFEM_ERROR("Filename cannot be empty after trimming.");

        // Reject dangerous characters
        static const std::regex allowedPattern(
            R"(^[a-zA-Z0-9_\-./\\]+$)");

        if (!std::regex_match(filename, allowedPattern))
            CFEM_ERROR("Filename contains invalid characters: " + filename);

        // Normalize extension
        std::filesystem::path p(filename);

        if (!p.has_extension())
        {
            p.replace_extension(requiredExtension);
        }
        else if (p.extension() != requiredExtension)
        {
            CFEM_ERROR("Invalid file extension. Expected " + requiredExtension);
        }

        return p.string();
    }
}
