#ifndef CFEM_VTKEXPORTER_H
#define CFEM_VTKEXPORTER_H

#include <fstream>
#include <vector>
#include "Mesh.h"
#include "Field.h" // Assuming this is where your Field base class lives

namespace cfem {

class VTKExporter {
public:
    VTKExporter(const Mesh& mesh) : m_mesh(mesh) {}

    void addField(std::shared_ptr<Field> field) {
        m_fields.push_back(field);
    }

    void setExportedEntities(const std::vector<std::string>& entityNames) {
        m_exportedEntityNames = entityNames;
    }
    void setExportQuality(bool enable) { m_exportQuality = enable; }

    void write(const std::string& filename);

private:
    const Mesh& m_mesh;
    bool m_exportQuality;
    std::vector<std::shared_ptr<Field>> m_fields;
    std::vector<std::string> m_exportedEntityNames;

    // Helper to map CellType to VTK cell identifiers
    int getVTKCellType(CellType type) const;
};

}

#endif