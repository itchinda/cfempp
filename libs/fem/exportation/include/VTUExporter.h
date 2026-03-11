#ifndef CFEM_VTUEXPORTER_H
#define CFEM_VTUEXPORTER_H

#include <vector>
#include <string>
#include <memory>
#include "Mesh.h"
#include "MeshExporter.h"

namespace cfem
{

    class VTUExporter : public MeshExporter
    {
    public:
        enum class Format
        {
            ASCII,
            BINARY
        };

        VTUExporter(const Mesh &mesh);
        void write(const std::string &filename) const override;
        void setFormat(Format format) {m_format = format;}

    private:
        Format m_format{Format::BINARY};
        void writeASCII(const std::string &filename) const;
        void writeBINARY(const std::string &filename) const;
    };

}
#endif