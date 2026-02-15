#include "VTUExporter.h"
#include "CfemExpression.h"
#include "FEField.h"
#include "LagrangeSpace.h"
#include "ReferenceElementFactory.h"

#include <algorithm>
#include <iomanip>

namespace cfem
{

    VTUExporter::VTUExporter(const Mesh &mesh) : MeshExporter(mesh)
    {  }

    void VTUExporter::write(const std::string &filename) const
    {
        const std::string validated = validateFilename(filename, ".vtu");
        if (m_format == Format::ASCII)
        {
            writeASCII(validated);
        }
        else
        {
            writeBINARY(validated);
        }
    }
    /**
     * @brief Exports the mesh and associated fields to a VTU Legacy (.vtu) ASCII file.
     *      This method follows the VTK Unstructured Grid specification. It handles:
     * 1. Geometry (Points and Cells)
     * 2. Field Interpolation (Mapping P1/P2 physics to geometric vertices)
     * 3. Metadata (Mesh quality and Entity masks)
     * @param filename Path to the output file.
     * @throws CFEM_ERROR If the file cannot be opened for writing.
     */
    void VTUExporter::writeASCII(const std::string &filename) const
    {
        std::ofstream out(filename);
        if (!out.is_open())
            CFEM_ERROR("Cannot open VTU file for writing: " + filename);

        const size_t nVerts = m_mesh.getNumVertices();
        const size_t nCells = m_mesh.getNumCells();

        // ============================================================
        // 1) WRITE XML HEADER
        // ============================================================

        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        out << "<UnstructuredGrid>\n";
        out << "<Piece NumberOfPoints=\"" << nVerts
            << "\" NumberOfCells=\"" << nCells << "\">\n";

        // ============================================================
        // 2) POINTS
        // ============================================================

        out << "<Points>\n";
        out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

        for (size_t i = 0; i < nVerts; ++i)
        {
            const auto &c = m_mesh.getVertex(i).coords;
            out << std::setprecision(16) << c[0] << " "
                << c[1] << " " << c[2] << "\n";
        }

        out << "</DataArray>\n";
        out << "</Points>\n";

        // ============================================================
        // 3) CELLS
        // ============================================================

        out << "<Cells>\n";

        // Connectivity
        out << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (size_t i = 0; i < nCells; ++i)
        {
            const auto &conn = m_mesh.getCellConnectivity(i);
            for (auto vid : conn)
                out << vid << " ";
            out << "\n";
        }
        out << "</DataArray>\n";

        // Offsets
        out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        int runningOffset = 0;
        for (size_t i = 0; i < nCells; ++i)
        {
            const auto &conn = m_mesh.getCellConnectivity(i);
            runningOffset += conn.size();
            out << runningOffset << "\n";
        }
        out << "</DataArray>\n";

        // Types
        out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (size_t i = 0; i < nCells; ++i)
        {
            out << static_cast<int>(getVTKCellType(m_mesh.getCell(i).type)) << "\n";
        }
        out << "</DataArray>\n";

        out << "</Cells>\n";

        // ============================================================
        // 4) POINT DATA
        // ============================================================

        if (!m_entries.empty())
        {
            out << "<PointData>\n";

            // Evaluate expressions at each vertex
            std::vector<std::vector<CfemReal>> allResults(m_entries.size());
            for (size_t i = 0; i < m_entries.size(); ++i)
            {
                m_entries[i].expr->prepare();
                allResults[i].resize(nVerts * m_entries[i].expr->getDimension(), 0.0);
            }

            std::vector<bool> computed(nVerts, false);
            EvalContext ctx;
            DynamicVector<Vector3D> physCoords;

            struct TypeCache
            {
                std::vector<ShapeInfo> nodeShapes;
                const ReferenceElement *refEl;
            };
            std::map<CellType, TypeCache> geoCache;

            for (CellType type : m_mesh.getUniqueCellTypes())
            {
                const auto &refEl = ReferenceElementFactory::getReferenceElement(type);
                auto &cache = geoCache[type];
                cache.refEl = &refEl;
                for (const auto &xi : refEl.referenceNodes())
                {
                    ShapeInfo info;
                    refEl.evaluate(xi, ShapeEvaluationFlags::Value, info);
                    cache.nodeShapes.push_back(info);
                }
            }

            for (CfemInt c = 0; c < (CfemInt)nCells; ++c)
            {
                const auto &cell = m_mesh.getCell(c);
                const auto &connectivityCell = m_mesh.getCellConnectivity(c);
                auto &cache = geoCache[cell.type];

                m_mesh.getCellVerticesCoords(c, physCoords);
                const auto &refEl = *cache.refEl;
                const auto &xiNodes = refEl.referenceNodes();

                for (size_t n = 0; n < connectivityCell.size(); ++n)
                {
                    CfemInt vid = connectivityCell[n];
                    if (computed[vid])
                        continue;

                    MappingInfo map;
                    refEl.computeMapping(xiNodes[n], physCoords, ShapeEvaluationFlags::None, cache.nodeShapes[n], map);
                    ctx.physPos = map.physicalPosition;

                    for (size_t i = 0; i < m_entries.size(); ++i)
                    {
                        int dim = m_entries[i].expr->getDimension();
                        std::span<CfemReal> outSpan(&allResults[i][vid * dim], dim);
                        m_entries[i].expr->evaluate(c, xiNodes[n], outSpan, ctx);
                    }

                    computed[vid] = true;
                }
            }

            // Write evaluated expressions
            for (size_t i = 0; i < m_entries.size(); ++i)
            {
                int dim = m_entries[i].expr->getDimension();
                out << "<DataArray type=\"Float64\" Name=\"" << m_entries[i].name
                    << "\" NumberOfComponents=\"" << dim << "\" format=\"ascii\">\n";

                for (size_t v = 0; v < nVerts; ++v)
                {
                    for (int d = 0; d < dim; ++d)
                        out << allResults[i][v * dim + d] << " ";
                    out << "\n";
                }

                out << "</DataArray>\n";
            }

            out << "</PointData>\n";
        }

        // ============================================================
        // 5) CELL DATA (optional quality, masks)
        // ============================================================

        if (m_exportQuality)
        {
            out << "<CellData>\n";
            out << "<DataArray type=\"Float32\" Name=\"Mesh_Quality\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            for (size_t i = 0; i < nCells; ++i)
                out << m_mesh.computeCellQuality(i) << "\n";
            out << "</DataArray>\n";
            out << "</CellData>\n";
        }

        // ============================================================
        // 6) CLOSE XML
        // ============================================================

        out << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    }

    /**
     * @brief Exports the mesh and associated fields to a VTU Legacy (.vtu) ASCII file.
     *      This method follows the VTK Unstructured Grid specification. It handles:
     * 1. Geometry (Points and Cells)
     * 2. Field Interpolation (Mapping P1/P2 physics to geometric vertices)
     * 3. Metadata (Mesh quality and Entity masks)
     * @param filename Path to the output file.
     * @throws CFEM_ERROR If the file cannot be opened for writing.
     */
    void VTUExporter::writeBINARY(const std::string &filename) const
    {
        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open())
            CFEM_ERROR("Cannot open VTU file for writing: " + filename);

        const size_t nVerts = m_mesh.getNumVertices();
        const size_t nCells = m_mesh.getNumCells();

        // ============================================================
        // 1) BUILD GEOMETRY BUFFERS
        // ============================================================

        std::vector<double> points(3 * nVerts);
        for (size_t i = 0; i < nVerts; ++i)
        {
            const auto &c = m_mesh.getVertex(i).coords;
            points[3 * i + 0] = c[0];
            points[3 * i + 1] = c[1];
            points[3 * i + 2] = c[2];
        }

        std::vector<int> connectivity;
        std::vector<int> offsets(nCells);
        std::vector<unsigned char> types(nCells);

        int runningOffset = 0;
        for (size_t i = 0; i < nCells; ++i)
        {
            const auto &conn = m_mesh.getCellConnectivity(i);
            for (auto vid : conn)
                connectivity.push_back((int)vid);

            runningOffset += conn.size();
            offsets[i] = runningOffset;
            types[i] = (unsigned char)getVTKCellType(m_mesh.getCell(i).type);
        }

        // ============================================================
        // 2) EVALUATE EXPRESSIONS (same logic as ASCII version)
        // ============================================================

        size_t numEntries = m_entries.size();
        std::vector<std::vector<CfemReal>> allResults(numEntries);

        for (size_t i = 0; i < numEntries; ++i)
        {
            m_entries[i].expr->prepare();
            allResults[i].resize(nVerts * m_entries[i].expr->getDimension(), 0.0);
        }

        std::vector<bool> computed(nVerts, false);
        EvalContext ctx;
        DynamicVector<Vector3D> physCoords;

        struct TypeCache
        {
            std::vector<ShapeInfo> nodeShapes;
            const ReferenceElement *refEl;
        };
        std::map<CellType, TypeCache> geoCache;

        for (CellType type : m_mesh.getUniqueCellTypes())
        {
            const auto &refEl = ReferenceElementFactory::getReferenceElement(type);
            auto &cache = geoCache[type];
            cache.refEl = &refEl;

            for (const auto &xi : refEl.referenceNodes())
            {
                ShapeInfo info;
                refEl.evaluate(xi, ShapeEvaluationFlags::Value, info);
                cache.nodeShapes.push_back(info);
            }
        }

        for (CfemInt c = 0; c < (CfemInt)nCells; ++c)
        {
            const auto &cell = m_mesh.getCell(c);
            const auto &connectivityCell = m_mesh.getCellConnectivity(c);
            auto &cache = geoCache[cell.type];

            m_mesh.getCellVerticesCoords(c, physCoords);

            const auto &refEl = *cache.refEl;
            const auto &xiNodes = refEl.referenceNodes();

            for (size_t n = 0; n < connectivityCell.size(); ++n)
            {
                CfemInt vid = connectivityCell[n];
                if (computed[vid])
                    continue;

                MappingInfo map;
                refEl.computeMapping(
                    xiNodes[n],
                    physCoords,
                    ShapeEvaluationFlags::None,
                    cache.nodeShapes[n],
                    map);

                ctx.physPos = map.physicalPosition;

                for (size_t i = 0; i < numEntries; ++i)
                {
                    int dim = m_entries[i].expr->getDimension();
                    std::span<CfemReal> outSpan(
                        &allResults[i][vid * dim], dim);

                    m_entries[i].expr->evaluate(c, xiNodes[n], outSpan, ctx);
                }

                computed[vid] = true;
            }
        }

        // ============================================================
        // 3) BUILD NODE MASKS
        // ============================================================

        std::vector<std::vector<float>> nodeMasks;

        for (const auto &entityName : m_exportedEntityNames)
        {
            MeshEntity *entity =
                const_cast<Mesh &>(m_mesh).getEntity(entityName);
            if (!entity)
                continue;

            std::vector<float> mask(nVerts, 0.0f);
            const auto &vids = entity->getVertexIds();

            for (CfemInt i = 0; i < (CfemInt)nVerts; ++i)
                if (std::binary_search(vids.begin(), vids.end(), i))
                    mask[i] = 1.0f;

            nodeMasks.push_back(std::move(mask));
        }

        // ============================================================
        // 4) BUILD CELL DATA (Quality + Masks)
        // ============================================================

        std::vector<float> quality;
        if (m_exportQuality)
        {
            quality.resize(nCells);
            for (size_t i = 0; i < nCells; ++i)
                quality[i] = m_mesh.computeCellQuality(i);
        }

        std::vector<std::vector<float>> cellMasks;

        for (const auto &entityName : m_exportedEntityNames)
        {
            MeshEntity *entity =
                const_cast<Mesh &>(m_mesh).getEntity(entityName);
            if (!entity)
                continue;

            std::vector<float> mask(nCells, 0.0f);
            const auto &cids = entity->getCellIds();

            for (CfemInt i = 0; i < (CfemInt)nCells; ++i)
                if (std::binary_search(cids.begin(), cids.end(), i))
                    mask[i] = 1.0f;

            cellMasks.push_back(std::move(mask));
        }

        // ============================================================
        // 5) APPENDED BINARY SECTION
        // ============================================================

        std::vector<char> appended;

        auto appendBlock = [&](auto &vec) -> size_t
        {
            size_t offset = appended.size();
            uint32_t bytes = static_cast<uint32_t>(
                vec.size() * sizeof(typename std::decay_t<decltype(vec)>::value_type));

            appended.insert(appended.end(),
                            reinterpret_cast<char *>(&bytes),
                            reinterpret_cast<char *>(&bytes) + sizeof(uint32_t));

            appended.insert(appended.end(),
                            reinterpret_cast<char *>(vec.data()),
                            reinterpret_cast<char *>(vec.data()) + bytes);

            return offset;
        };

        size_t pointsOffset = appendBlock(points);
        size_t connOffset = appendBlock(connectivity);
        size_t offOffset = appendBlock(offsets);
        size_t typeOffset = appendBlock(types);

        std::vector<size_t> exprOffsets;
        for (auto &r : allResults)
            exprOffsets.push_back(appendBlock(r));

        std::vector<size_t> nodeMaskOffsets;
        for (auto &m : nodeMasks)
            nodeMaskOffsets.push_back(appendBlock(m));

        size_t qualityOffset = 0;
        if (m_exportQuality)
            qualityOffset = appendBlock(quality);

        std::vector<size_t> cellMaskOffsets;
        for (auto &m : cellMasks)
            cellMaskOffsets.push_back(appendBlock(m));

        // ============================================================
        // 6) WRITE XML HEADER
        // ============================================================

        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
            << "byte_order=\"LittleEndian\">\n";
        out << "<UnstructuredGrid>\n";
        out << "<Piece NumberOfPoints=\"" << nVerts
            << "\" NumberOfCells=\"" << nCells << "\">\n";

        // Points
        out << "<Points>\n";
        out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" "
            << "format=\"appended\" offset=\"" << pointsOffset << "\"/>\n";
        out << "</Points>\n";

        // Cells
        out << "<Cells>\n";
        out << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"" << connOffset << "\"/>\n";
        out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"" << offOffset << "\"/>\n";
        out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << typeOffset << "\"/>\n";
        out << "</Cells>\n";

        // PointData
        if (!m_entries.empty() || !nodeMasks.empty())
        {
            out << "<PointData>\n";

            for (size_t i = 0; i < m_entries.size(); ++i)
            {
                int dim = m_entries[i].expr->getDimension();
                out << "<DataArray type=\"Float64\" Name=\"" << m_entries[i].name
                    << "\" NumberOfComponents=\"" << dim
                    << "\" format=\"appended\" offset=\""
                    << exprOffsets[i] << "\"/>\n";
            }

            for (size_t i = 0; i < nodeMasks.size(); ++i)
            {
                out << "<DataArray type=\"Float32\" Name=\"NodeMask_"
                    << m_exportedEntityNames[i]
                    << "\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                    << nodeMaskOffsets[i] << "\"/>\n";
            }

            out << "</PointData>\n";
        }

        // CellData
        if (m_exportQuality || !cellMasks.empty())
        {
            out << "<CellData>\n";

            if (m_exportQuality)
                out << "<DataArray type=\"Float32\" Name=\"Mesh_Quality\" "
                    << "NumberOfComponents=\"1\" format=\"appended\" offset=\""
                    << qualityOffset << "\"/>\n";

            for (size_t i = 0; i < cellMasks.size(); ++i)
                out << "<DataArray type=\"Float32\" Name=\"CellMask_"
                    << m_exportedEntityNames[i]
                    << "\" NumberOfComponents=\"1\" format=\"appended\" offset=\""
                    << cellMaskOffsets[i] << "\"/>\n";

            out << "</CellData>\n";
        }

        out << "</Piece>\n</UnstructuredGrid>\n";

        out << "<AppendedData encoding=\"raw\">\n_";
        out.write(appended.data(), appended.size());
        out << "\n</AppendedData>\n</VTKFile>\n";
    }

}
