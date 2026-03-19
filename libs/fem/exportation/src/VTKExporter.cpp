/**
 * @file VTKExporter.cpp
 * @brief Implementation of the VTK Legacy format exporter for CFEM simulations.
 */

#include "VTKExporter.h"
#include "CfemExpression.h"
#include "FEField.h"
#include "LagrangeSpace.h"
#include "Mesh.h"
#include "ReferenceElementFactory.h"
#include <algorithm>
#include <iomanip>

namespace cfem
{
    VTKExporter::VTKExporter(const Mesh &mesh) : MeshExporter(mesh)
    { }

    void VTKExporter::write(const std::string &filename) const
    {
        const std::string validated = validateFilename(filename, ".vtk");
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
     * @brief Exports the mesh and associated fields to a VTK Legacy (.vtk) ASCII file.
     * * This method follows the VTK Unstructured Grid specification. It handles:
     * 1. Geometry (Points and Cells)
     * 2. Field Interpolation (Mapping P1/P2 physics to geometric vertices)
     * 3. Metadata (Mesh quality and Entity masks)
     * * @param filename Path to the output file.
     * @throws CFEM_ERROR If the file cannot be opened for writing.
     */
    void VTKExporter::writeASCII(const std::string &filename) const
    {

        std::ofstream out(filename);
        if (!out.is_open())
        {
            CFEM_ERROR << "Cannot open VTK file for writing: " << filename;
        }

        // High precision is required to avoid artifacts in derivative calculations in ParaView
        out << std::scientific << std::setprecision(12);

        // --- 1. HEADER ---
        out << "# vtk DataFile Version 3.0\n";
        out << "CFEM Simulation Export\n";
        out << "ASCII\n";
        out << "DATASET UNSTRUCTURED_GRID\n\n";

        // --- 2. POINTS ---
        out << "POINTS " << m_mesh.getNumVertices() << " double\n";
        for (size_t i = 0; i < m_mesh.getNumVertices(); ++i)
        {
            const auto &coords = m_mesh.getVertex(i).coords;
            out << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
        }

        // --- 3. CELLS ---
        // Calculate total size: sum of (1 + number of nodes per cell)
        size_t totalCellSize = 0;
        for (size_t i = 0; i < m_mesh.getNumCells(); ++i)
        {
            totalCellSize += 1 + m_mesh.getCellConnectivity(i).size();
        }

        out << "\nCELLS " << m_mesh.getNumCells() << " " << totalCellSize << "\n";
        for (size_t i = 0; i < m_mesh.getNumCells(); ++i)
        {
            const auto &cell_ids = m_mesh.getCellConnectivity(i);
            out << cell_ids.size();
            for (auto vid : cell_ids)
                out << " " << vid;
            out << "\n";
        }

        // --- 4. CELL TYPES ---
        out << "\nCELL_TYPES " << m_mesh.getNumCells() << "\n";
        for (size_t i = 0; i < m_mesh.getNumCells(); ++i)
        {
            out << getVTKCellType(m_mesh.getCell(i).type) << "\n";
        }

        // --- 2. POINT DATA OPTIMIZATION ---
        if (!m_entries.empty() || !m_exportedEntityNames.empty())
        {
            out << "\nPOINT_DATA " << m_mesh.getNumVertices() << "\n";

            // Pre-prepare all expressions and allocate all result buffers
            size_t numEntries = m_entries.size();
            std::vector<std::vector<CfemReal>> allResults(numEntries);
            for (size_t i = 0; i < numEntries; ++i)
            {
                m_entries[i].expr->prepare();
                allResults[i].resize(m_mesh.getNumVertices() * m_entries[i].expr->getDimension(), 0.0);
            }

            std::vector<bool> computed(m_mesh.getNumVertices(), false);
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

            // SINGLE PASS OVER THE MESH
            for (CfemInt c = 0; c < (CfemInt)m_mesh.getNumCells(); ++c)
            {
                const auto &cell = m_mesh.getCell(c);
                const auto &connectivity = m_mesh.getCellConnectivity(c);
                auto &cache = geoCache[cell.type];

                m_mesh.getCellVerticesCoords(c, physCoords);

                const auto &refEl = ReferenceElementFactory::getReferenceElement(cell.type);
                const auto &xiNodes = refEl.referenceNodes();

                for (size_t n = 0; n < connectivity.size(); ++n)
                {
                    CfemInt vid = connectivity[n];
                    if (computed[vid])
                        continue;

                    ctx.resetBuffers(); 

                    refEl.computeMapping(xiNodes[n],
                                        physCoords,
                                        ShapeEvaluationFlags::None,
                                        cache.nodeShapes[n], ctx.scratchMapping);

                    ctx.physPos = ctx.scratchMapping.physicalPosition;

                    for (size_t i = 0; i < numEntries; ++i)
                    {
                        int dim = m_entries[i].expr->getDimension();
                        std::span<CfemReal> outSpan(&allResults[i][vid * dim], dim);
                        m_entries[i].expr->evaluate(c, xiNodes[n], outSpan, ctx);
                    }
                    computed[vid] = true;
                }
            }

            // WRITE ALL EVALUATED DATA TO FILE
            for (size_t i = 0; i < numEntries; ++i)
            {
                int dim = m_entries[i].expr->getDimension();
                const std::string &name = m_entries[i].name;

                if (dim == 1)
                {
                    out << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
                    for (size_t j = 0; j < m_mesh.getNumVertices(); ++j)
                        out << allResults[i][j] << "\n";
                }
                else
                {
                    out << "VECTORS " << name << " double\n";
                    for (size_t j = 0; j < m_mesh.getNumVertices(); ++j)
                    {
                        out << allResults[i][j * dim] << " "
                            << (dim >= 2 ? allResults[i][j * dim + 1] : 0.0) << " "
                            << (dim >= 3 ? allResults[i][j * dim + 2] : 0.0) << "\n";
                    }
                }
            }

            // B. Node Masks (Boundary Conditions / Entity Sets)
            for (const auto &entityName : m_exportedEntityNames)
            {
                MeshEntity *entity = const_cast<Mesh &>(m_mesh).getEntity(entityName);
                if (!entity)
                    continue;

                std::string safeName = "NodeMask_" + entityName;
                std::replace(safeName.begin(), safeName.end(), ' ', '_');
                out << "SCALARS " << safeName << " float 1\nLOOKUP_TABLE default\n";

                const auto &vids = entity->getVertexIds();
                for (CfemInt i = 0; i < (CfemInt)m_mesh.getNumVertices(); ++i)
                {
                    out << (std::binary_search(vids.begin(), vids.end(), i) ? "1.0" : "0.0") << "\n";
                }
            }
        }

        // --- 3. CELL DATA (Quality + Cell Masks) ---
        bool hasCellData = m_exportQuality || !m_exportedEntityNames.empty();
        if (hasCellData)
        {
            out << "\nCELL_DATA " << m_mesh.getNumCells() << "\n";

            if (m_exportQuality)
            {
                out << "SCALARS Mesh_Quality float 1\nLOOKUP_TABLE default\n";
                for (size_t i = 0; i < m_mesh.getNumCells(); ++i)
                {
                    out << m_mesh.computeCellQuality(i) << "\n";
                }
            }

            for (const auto &entityName : m_exportedEntityNames)
            {
                MeshEntity *entity = const_cast<Mesh &>(m_mesh).getEntity(entityName);
                if (!entity)
                    continue;

                std::string safeName = "CellMask_" + entityName;
                std::replace(safeName.begin(), safeName.end(), ' ', '_');
                out << "SCALARS " << safeName << " float 1\nLOOKUP_TABLE default\n";

                const auto &cids = entity->getCellIds();
                for (CfemInt i = 0; i < (CfemInt)m_mesh.getNumCells(); ++i)
                {
                    out << (std::binary_search(cids.begin(), cids.end(), i) ? "1.0" : "0.0") << "\n";
                }
            }
        }

        out.close();
    }

    /**
     * @brief Exports the mesh and associated fields to a VTK Legacy (.vtk) binary file.
     * * This method follows the VTK Unstructured Grid specification. It handles:
     * 1. Geometry (Points and Cells)
     * 2. Field Interpolation (Mapping P1/P2 physics to geometric vertices)
     * 3. Metadata (Mesh quality and Entity masks)
     * * @param filename Path to the output file.
     * @throws CFEM_ERROR If the file cannot be opened for writing.
     */
    void VTKExporter::writeBINARY(const std::string &filename) const
    {
        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open())
            CFEM_ERROR << "Cannot open VTK file for writing: " << filename;

        auto writeBigEndian = [&](auto *data, std::size_t count)
        {
            using T = std::remove_pointer_t<decltype(data)>;
            for (std::size_t i = 0; i < count; ++i)
            {
                T value = data[i];

#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
                char *p = reinterpret_cast<char *>(&value);
                std::reverse(p, p + sizeof(T));
#endif
                out.write(reinterpret_cast<char *>(&value), sizeof(T));
            }
        };

        // --------------------------------------------------
        // HEADER (ASCII)
        // --------------------------------------------------

        out << "# vtk DataFile Version 3.0\n";
        out << "CFEM Simulation Export\n";
        out << "BINARY\n";
        out << "DATASET UNSTRUCTURED_GRID\n";

        // --------------------------------------------------
        // POINTS
        // --------------------------------------------------

        std::size_t nVerts = m_mesh.getNumVertices();
        out << "POINTS " << nVerts << " double\n";

        std::vector<double> pointBuffer(3 * nVerts);

        for (std::size_t i = 0; i < nVerts; ++i)
        {
            const auto &coords = m_mesh.getVertex(i).coords;
            pointBuffer[3 * i + 0] = coords[0];
            pointBuffer[3 * i + 1] = coords[1];
            pointBuffer[3 * i + 2] = coords[2];
        }

        writeBigEndian(pointBuffer.data(), pointBuffer.size());
        out << "\n";

        // --------------------------------------------------
        // CELLS
        // --------------------------------------------------

        std::size_t nCells = m_mesh.getNumCells();
        std::size_t totalCellSize = 0;

        for (std::size_t i = 0; i < nCells; ++i)
            totalCellSize += 1 + m_mesh.getCellConnectivity(i).size();

        out << "CELLS " << nCells << " " << totalCellSize << "\n";

        std::vector<int> cellBuffer;
        cellBuffer.reserve(totalCellSize);

        for (std::size_t i = 0; i < nCells; ++i)
        {
            const auto &conn = m_mesh.getCellConnectivity(i);
            cellBuffer.push_back((int)conn.size());
            for (auto vid : conn)
                cellBuffer.push_back((int)vid);
        }

        writeBigEndian(cellBuffer.data(), cellBuffer.size());
        out << "\n";

        // --------------------------------------------------
        // CELL TYPES
        // --------------------------------------------------

        out << "CELL_TYPES " << nCells << "\n";

        std::vector<int> typeBuffer(nCells);

        for (std::size_t i = 0; i < nCells; ++i)
            typeBuffer[i] = getVTKCellType(m_mesh.getCell(i).type);

        writeBigEndian(typeBuffer.data(), typeBuffer.size());
        out << "\n";

        // --------------------------------------------------
        // POINT DATA (UNCHANGED EVALUATION LOGIC)
        // --------------------------------------------------

        if (!m_entries.empty() || !m_exportedEntityNames.empty())
        {
            out << "POINT_DATA " << nVerts << "\n";

            std::size_t numEntries = m_entries.size();
            std::vector<std::vector<double>> allResults(numEntries);

            for (std::size_t i = 0; i < numEntries; ++i)
            {
                m_entries[i].expr->prepare();
                allResults[i].resize(nVerts * m_entries[i].expr->getDimension(), 0.0);
            }

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

            std::vector<bool> computed(nVerts, false);
            EvalContext ctx;
            DynamicVector<Vector3D> physCoords;

            for (CfemInt c = 0; c < (CfemInt)nCells; ++c)
            {
                const auto &cell = m_mesh.getCell(c);
                const auto &connectivity = m_mesh.getCellConnectivity(c);
                const auto &refEl = ReferenceElementFactory::getReferenceElement(cell.type);
                const auto &xiNodes = refEl.referenceNodes();
                auto &cache = geoCache[cell.type];

                m_mesh.getCellVerticesCoords(c, physCoords);

                for (size_t n = 0; n < connectivity.size(); ++n)
                {
                    CfemInt vid = connectivity[n];
                    if (computed[vid])
                        continue;

                    ctx.resetBuffers(); 

                    refEl.computeMapping(xiNodes[n],
                                        physCoords,
                                        ShapeEvaluationFlags::None,
                                        cache.nodeShapes[n], ctx.scratchMapping);

                    ctx.physPos = ctx.scratchMapping.physicalPosition;

                    for (size_t i = 0; i < numEntries; ++i)
                    {
                        int dim = m_entries[i].expr->getDimension();
                        std::span<CfemReal> outSpan(&allResults[i][vid * dim], dim);
                        m_entries[i].expr->evaluate(c, xiNodes[n], outSpan, ctx);
                    }
                    computed[vid] = true;
                }
            }

            // ---- Write evaluated fields in binary ----

            for (std::size_t i = 0; i < numEntries; ++i)
            {
                int dim = m_entries[i].expr->getDimension();
                const std::string &name = m_entries[i].name;

                if (dim == 1)
                {
                    out << "SCALARS " << name << " double 1\n";
                    out << "LOOKUP_TABLE default\n";
                    writeBigEndian(allResults[i].data(), allResults[i].size());
                    out << "\n";
                }
                else
                {
                    out << "VECTORS " << name << " double\n";

                    std::vector<double> vecBuffer(3 * nVerts, 0.0);

                    for (std::size_t j = 0; j < nVerts; ++j)
                    {
                        for (int d = 0; d < dim && d < 3; ++d)
                            vecBuffer[3 * j + d] = allResults[i][j * dim + d];
                    }

                    writeBigEndian(vecBuffer.data(), vecBuffer.size());
                    out << "\n";
                }
            }
        }

        // --------------------------------------------------
        // CELL DATA
        // --------------------------------------------------

        if (m_exportQuality)
        {
            out << "CELL_DATA " << nCells << "\n";
            out << "SCALARS Mesh_Quality float 1\n";
            out << "LOOKUP_TABLE default\n";

            std::vector<float> quality(nCells);
            for (std::size_t i = 0; i < nCells; ++i)
                quality[i] = (float)m_mesh.computeCellQuality(i);

            writeBigEndian(quality.data(), quality.size());
            out << "\n";
        }

        out.close();
    }

} // namespace cfem