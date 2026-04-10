#include "VTUExporter.h"
#include "CfemExpression.h"
#include "FEField.h"
#include "LagrangeSpace.h"
#include "Mesh.h"
#include "ReferenceElementFactory.h"

#include <algorithm>
#include <iomanip>

namespace cfem
{

    VTUExporter::VTUExporter(const Mesh &mesh) : MeshExporter(mesh)
    {  }

    // void VTUExporter::write(const std::string &filename) const
    // {
    //     const std::string validated = validateFilename(filename, ".vtu");
    //     if (m_format == Format::ASCII)
    //     {
    //         writeASCII(validated);
    //     }
    //     else
    //     {
    //         writeBINARY(validated);
    //     }
    // }


    void VTUExporter::write(const std::string &filename) const
    {
        std::vector<ExportEntry> continuousFields;
        std::vector<ExportEntry> discontinuousFields;

        for (const auto& entry : m_entries) {
            // On suppose que l'expression peut nous dire si elle est continue
            // ou on vérifie le type de FEField/FESpace associé
            if (entry.expr->isContinuous()) {
                continuousFields.push_back(entry);
            } else {
                CFEM_INFO << "Discontinuous field detected";
                discontinuousFields.push_back(entry);
            }
        }

        // Export Continu
        if (!continuousFields.empty()) {
            std::string contFile = validateFilename(filename, ".cont.vtu");

            if (m_format == Format::ASCII){
                writeASCIIContinuous(contFile, continuousFields);
            }
            else{
                writeBINARYContinuous(contFile, continuousFields);
            }

        }

        // Export Discontinu
        if (!discontinuousFields.empty()) {
            std::string discFile = validateFilename(filename, ".disc.vtu");
            if (m_format == Format::ASCII){
                writeASCIIDiscontinuous(discFile, discontinuousFields);
            }
            else{
                writeBINARYDiscontinuous(discFile, discontinuousFields);
            }
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
    void VTUExporter::writeASCIIContinuous(const std::string &filename,
                                           const std::vector<ExportEntry>& continuousFields) const
    {
        std::ofstream out(filename);
        if (!out.is_open())
            CFEM_ERROR << "Cannot open VTU file for writing: " << filename;

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

        if (!continuousFields.empty())
        {
            out << "<PointData>\n";

            // Evaluate expressions at each vertex
            std::vector<std::vector<CfemReal>> allResults(continuousFields.size());
            for (size_t i = 0; i < continuousFields.size(); ++i)
            {
                continuousFields[i].expr->prepare();
                allResults[i].resize(nVerts * continuousFields[i].expr->getDimension(), 0.0);
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
            CfemInt geoOrder = m_mesh.getOrderInt();

            for (CellType type : m_mesh.getCellTypesPresent())
            {
                const auto &refEl = ReferenceElementFactory::getReferenceElement(type, geoOrder);
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

                    ctx.resetBuffers(); 

                    refEl.computeMapping(xiNodes[n],
                                        physCoords,
                                        ShapeEvaluationFlags::PhysicalPosition,
                                        cache.nodeShapes[n], ctx.scratchMapping);

                    ctx.physPos = ctx.scratchMapping.physicalPosition;

                    for (size_t i = 0; i < continuousFields.size(); ++i)
                    {
                        int dim = continuousFields[i].expr->getDimension();
                        std::span<CfemReal> outSpan(&allResults[i][vid * dim], dim);
                        continuousFields[i].expr->evaluate(c, xiNodes[n], outSpan, ctx);
                    }
                    computed[vid] = true;
                }
            }

            // Write evaluated expressions
            for (size_t i = 0; i < continuousFields.size(); ++i)
            {
                int dim = continuousFields[i].expr->getDimension();
                out << "<DataArray type=\"Float64\" Name=\"" << continuousFields[i].name
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
    void VTUExporter::writeBINARYContinuous(const std::string &filename,
                                            const std::vector<ExportEntry>& continuousFields) const
    {
        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open())
            CFEM_ERROR << "Cannot open VTU file for writing: " << filename;

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

        size_t numEntries = continuousFields.size();
        std::vector<std::vector<CfemReal>> allResults(numEntries);

        for (size_t i = 0; i < numEntries; ++i)
        {
            continuousFields[i].expr->prepare();
            allResults[i].resize(nVerts * continuousFields[i].expr->getDimension(), 0.0);
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
        CfemInt geoOrder = m_mesh.getOrderInt();

        for (CellType type : m_mesh.getCellTypesPresent())
        {
            const auto &refEl = ReferenceElementFactory::getReferenceElement(type, geoOrder);
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

                ctx.resetBuffers(); 

                refEl.computeMapping(xiNodes[n],
                                     physCoords,
                                     ShapeEvaluationFlags::PhysicalPosition,
                                     cache.nodeShapes[n], ctx.scratchMapping);

                ctx.physPos = ctx.scratchMapping.physicalPosition;

                for (size_t i = 0; i < continuousFields.size(); ++i)
                {
                    int dim = continuousFields[i].expr->getDimension();
                    std::span<CfemReal> outSpan(&allResults[i][vid * dim], dim);
                    continuousFields[i].expr->evaluate(c, xiNodes[n], outSpan, ctx);
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
        if (!continuousFields.empty() || !nodeMasks.empty())
        {
            out << "<PointData>\n";

            for (size_t i = 0; i < continuousFields.size(); ++i)
            {
                int dim = continuousFields[i].expr->getDimension();
                out << "<DataArray type=\"Float64\" Name=\"" << continuousFields[i].name
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

    /**
     * @brief Exports the discontinuous fields to an ASCII VTU file.
     * Nodes are duplicated for each cell to allow discontinuous visualization.
     */
    void VTUExporter::writeASCIIDiscontinuous(const std::string &filename,
                                              const std::vector<ExportEntry>& discontFields) const
    {
        CFEM_INFO << "Exporting discontinuous fields";
        std::ofstream out(filename);
        if (!out.is_open())
            CFEM_ERROR << "Cannot open VTU file for writing: " << filename;

        const size_t nCells = m_mesh.getNumCells();
        
        // 1) CALCULATE TOTAL DUPLICATED POINTS
        size_t totalDgPoints = 0;
        for (size_t i = 0; i < nCells; ++i) {
            totalDgPoints += m_mesh.getCellConnectivity(i).size();
        }

        // 2) WRITE XML HEADER
        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        out << "<UnstructuredGrid>\n";
        out << "<Piece NumberOfPoints=\"" << totalDgPoints << "\" NumberOfCells=\"" << nCells << "\">\n";

        // 3) POINTS (Duplicated)
        out << "<Points>\n";
        out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (size_t c = 0; c < nCells; ++c) {
            const auto& connectivity = m_mesh.getCellConnectivity(c);
            for (CfemInt vId : connectivity) {
                const auto& coords = m_mesh.getVertex(vId).coords;
                out << std::setprecision(16) << coords[0] << " " << coords[1] << " " << coords[2] << "\n";
            }
        }
        out << "</DataArray>\n";
        out << "</Points>\n";

        // 4) CELLS
        out << "<Cells>\n";
        
        // Connectivity (Linear mapping)
        out << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (size_t i = 0; i < totalDgPoints; ++i) {
            out << i << " ";
        }
        out << "\n</DataArray>\n";

        // Offsets
        out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        int runningOffset = 0;
        for (size_t c = 0; c < nCells; ++c) {
            runningOffset += m_mesh.getCellConnectivity(c).size();
            out << runningOffset << "\n";
        }
        out << "</DataArray>\n";

        // Types
        out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
        for (size_t c = 0; c < nCells; ++c) {
            out << static_cast<int>(getVTKCellType(m_mesh.getCell(c).type)) << "\n";
        }
        out << "</DataArray>\n";
        out << "</Cells>\n";

        // 5) POINT DATA (Evaluations)
        if (!discontFields.empty()) {
            out << "<PointData>\n";
            
            EvalContext ctx;
            DynamicVector<Vector3D> physCoords;
            CfemInt geoOrder = m_mesh.getOrderInt();

            for (const auto& entry : discontFields) {
                entry.expr->prepare();
                int dim = entry.expr->getDimension();
                
                out << "<DataArray type=\"Float64\" Name=\"" << entry.name 
                    << "\" NumberOfComponents=\"" << dim << "\" format=\"ascii\">\n";

                for (CfemInt c = 0; c < (CfemInt)nCells; ++c) {
                    const auto& cell = m_mesh.getCell(c);
                    const auto& refEl = ReferenceElementFactory::getReferenceElement(cell.type, geoOrder);
                    const auto& xiNodes = refEl.referenceNodes();
                    
                    m_mesh.getCellVerticesCoords(c, physCoords);

                    for (size_t n = 0; n < xiNodes.size(); ++n) {
                        ctx.resetBuffers();
                        ShapeInfo info; // We need basic shape info to evaluate properly mapping
                        refEl.evaluate(xiNodes[n], ShapeEvaluationFlags::Value, info);
                        refEl.computeMapping(xiNodes[n], physCoords, ShapeEvaluationFlags::PhysicalPosition, info, ctx.scratchMapping);
                        ctx.physPos = ctx.scratchMapping.physicalPosition;

                        std::vector<CfemReal> values(dim, 0.0);
                        std::span<CfemReal> outSpan(values.data(), dim);
                        entry.expr->evaluate(c, xiNodes[n], outSpan, ctx);

                        for (int d = 0; d < dim; ++d) {
                            out << values[d] << " ";
                        }
                        out << "\n";
                    }
                }
                out << "</DataArray>\n";
            }
            out << "</PointData>\n";
        }

        // 6) CLOSE XML
        out << "</Piece>\n</UnstructuredGrid>\n</VTKFile>\n";
    }
    
    /**
     * @brief Exports the discontinuous fields to a Binary VTU file.
     */
    void VTUExporter::writeBINARYDiscontinuous(const std::string &filename,
                                               const std::vector<ExportEntry>& discontFields) const
    {
        CFEM_INFO << "Exporting discontinuous fields";
        std::ofstream out(filename, std::ios::binary);
        if (!out.is_open())
            CFEM_ERROR << "Cannot open VTU file for writing: " << filename;

        const size_t nCells = m_mesh.getNumCells();

        // 1) CALCULATE TOTAL DUPLICATED POINTS
        size_t totalDgPoints = 0;
        for (size_t i = 0; i < nCells; ++i) {
            totalDgPoints += m_mesh.getCellConnectivity(i).size();
        }

        // 2) BUILD GEOMETRY BUFFERS
        std::vector<double> points(3 * totalDgPoints);
        std::vector<int> connectivity(totalDgPoints);
        std::vector<int> offsets(nCells);
        std::vector<unsigned char> types(nCells);

        int runningOffset = 0;
        int ptIdx = 0;

        for (size_t c = 0; c < nCells; ++c) {
            const auto &conn = m_mesh.getCellConnectivity(c);
            for (auto vid : conn) {
                const auto &coord = m_mesh.getVertex(vid).coords;
                points[3 * ptIdx + 0] = coord[0];
                points[3 * ptIdx + 1] = coord[1];
                points[3 * ptIdx + 2] = coord[2];
                
                connectivity[ptIdx] = ptIdx; // Linear mapping for duplicated points
                ptIdx++;
            }
            runningOffset += conn.size();
            offsets[c] = runningOffset;
            types[c] = (unsigned char)getVTKCellType(m_mesh.getCell(c).type);
        }

        // 3) EVALUATE EXPRESSIONS
        size_t numEntries = discontFields.size();
        std::vector<std::vector<CfemReal>> allResults(numEntries);

        for (size_t i = 0; i < numEntries; ++i) {
            discontFields[i].expr->prepare();
            allResults[i].resize(totalDgPoints * discontFields[i].expr->getDimension(), 0.0);
        }

        EvalContext ctx;
        DynamicVector<Vector3D> physCoords;
        CfemInt geoOrder = m_mesh.getOrderInt();
        int globalPtIdx = 0;

        for (CfemInt c = 0; c < (CfemInt)nCells; ++c) {
            const auto &cell = m_mesh.getCell(c);
            const auto &refEl = ReferenceElementFactory::getReferenceElement(cell.type, geoOrder);
            const auto &xiNodes = refEl.referenceNodes();

            m_mesh.getCellVerticesCoords(c, physCoords);

            for (size_t n = 0; n < xiNodes.size(); ++n) {
                ctx.resetBuffers();
                ShapeInfo info;
                refEl.evaluate(xiNodes[n], ShapeEvaluationFlags::Value, info);
                refEl.computeMapping(xiNodes[n], physCoords, ShapeEvaluationFlags::PhysicalPosition, info, ctx.scratchMapping);
                ctx.physPos = ctx.scratchMapping.physicalPosition;

                for (size_t i = 0; i < numEntries; ++i) {
                    int dim = discontFields[i].expr->getDimension();
                    std::span<CfemReal> outSpan(&allResults[i][globalPtIdx * dim], dim);
                    discontFields[i].expr->evaluate(c, xiNodes[n], outSpan, ctx);
                }
                globalPtIdx++;
            }
        }

        // 4) APPENDED BINARY SECTION
        std::vector<char> appended;

        auto appendBlock = [&](auto &vec) -> size_t {
            size_t offset = appended.size();
            uint32_t bytes = static_cast<uint32_t>(vec.size() * sizeof(typename std::decay_t<decltype(vec)>::value_type));

            appended.insert(appended.end(),
                            reinterpret_cast<char *>(&bytes),
                            reinterpret_cast<char *>(&bytes) + sizeof(uint32_t));

            appended.insert(appended.end(),
                            reinterpret_cast<char *>(vec.data()),
                            reinterpret_cast<char *>(vec.data()) + bytes);

            return offset;
        };

        size_t pointsOffset = appendBlock(points);
        size_t connOffset   = appendBlock(connectivity);
        size_t offOffset    = appendBlock(offsets);
        size_t typeOffset   = appendBlock(types);

        std::vector<size_t> exprOffsets;
        for (auto &r : allResults) {
            exprOffsets.push_back(appendBlock(r));
        }

        // 5) WRITE XML HEADER
        out << "<?xml version=\"1.0\"?>\n";
        out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        out << "<UnstructuredGrid>\n";
        out << "<Piece NumberOfPoints=\"" << totalDgPoints << "\" NumberOfCells=\"" << nCells << "\">\n";

        // Points
        out << "<Points>\n";
        out << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"appended\" offset=\"" << pointsOffset << "\"/>\n";
        out << "</Points>\n";

        // Cells
        out << "<Cells>\n";
        out << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\"" << connOffset << "\"/>\n";
        out << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\"" << offOffset << "\"/>\n";
        out << "<DataArray type=\"UInt8\" Name=\"types\" format=\"appended\" offset=\"" << typeOffset << "\"/>\n";
        out << "</Cells>\n";

        // PointData
        if (!discontFields.empty()) {
            out << "<PointData>\n";
            for (size_t i = 0; i < numEntries; ++i) {
                int dim = discontFields[i].expr->getDimension();
                out << "<DataArray type=\"Float64\" Name=\"" << discontFields[i].name
                    << "\" NumberOfComponents=\"" << dim
                    << "\" format=\"appended\" offset=\"" << exprOffsets[i] << "\"/>\n";
            }
            out << "</PointData>\n";
        }

        out << "</Piece>\n</UnstructuredGrid>\n";

        // 6) APPEND RAW DATA
        out << "<AppendedData encoding=\"raw\">\n_";
        out.write(appended.data(), appended.size());
        out << "\n</AppendedData>\n</VTKFile>\n";
    }


}