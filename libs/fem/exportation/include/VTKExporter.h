#ifndef CFEM_VTKEXPORTER_H
#define CFEM_VTKEXPORTER_H

#include <fstream>
#include <vector>
#include <string>
#include <memory>
#include "Mesh.h"
#include "FEField.h"

namespace cfem {

class VTKExporter {
public:
    VTKExporter(const Mesh& mesh) : m_mesh(mesh), m_exportQuality(false) {}

     void addField(const FEField& field) {
        m_fields.push_back(&field);
    }

    void setExportedEntities(const std::vector<std::string>& entityNames) {
        m_exportedEntityNames = entityNames;
    }

    void setExportQuality(bool enable) { m_exportQuality = enable; }

    void write(const std::string& filename);

private:
    const Mesh& m_mesh;
    bool m_exportQuality;
    std::vector<const FEField*> m_fields;
    std::vector<std::string> m_exportedEntityNames;

    int getVTKCellType(CellType type) const;
};

} // namespace cfem

#endif