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

#include "cfemutils.h"

#include "VTKFamilyExporter.h"
#include "CfemExpression.h"
#include "FEField.h"
#include "LagrangeSpace.h"
#include "Mesh.h"
#include "ReferenceElementFactory.h"


namespace cfem
{
    void VTKFamilyExporter::write(const std::string &filename, CfemInt entityId) const
    {
        std::vector<ExportEntry> continuousFields;
        std::vector<ExportEntry> discontinuousFields;

        for (const auto& entry : m_entries) {
            if (entry.expr->isContinuous()) {
                continuousFields.push_back(entry);
            } else {
                CFEM_INFO << "Discontinuous field detected: " << entry.name;
                discontinuousFields.push_back(entry);
            }
        }

        if (!continuousFields.empty()) {
            std::string contFile = validateFilename(filename, ".cont." + m_outputFileExtension);
            if (m_format == VTKExportFormat::ASCII) writeASCIIContinuous(contFile, continuousFields, entityId);
            else                                    writeBINARYContinuous(contFile, continuousFields, entityId);
        }

        if (!discontinuousFields.empty()) {
            std::string discFile = validateFilename(filename, ".disc." + m_outputFileExtension);
            if (m_format == VTKExportFormat::ASCII) writeASCIIDiscontinuous(discFile, discontinuousFields, entityId);
            else                                    writeBINARYDiscontinuous(discFile, discontinuousFields, entityId);
        }
    }

    // ========================================================================
    // --- KERNEL 1: GEOMETRY EXTRACTION (CONTINUOUS)
    // ========================================================================
    VTKFamilyExporter::ExtractedGeometry VTKFamilyExporter::buildContinuousGeometry(CfemInt entityId) const
    {
        ExtractedGeometry geo;
        const MeshEntity& entity = m_mesh.getEntityConst(entityId);
        bool isVolume = (entity.getDimension() == m_mesh.getDimension());
        
        geo.pActiveElements = isVolume ? &entity.getCellIds() : &entity.getFacetIds();
        geo.numActiveCells = geo.pActiveElements->size();

        geo.globalToLocal.assign(m_mesh.getNumVertices(), -1);
        geo.localToGlobal.reserve(geo.numActiveCells * 4); 

        int runningOffset = 0;

        for (CfemInt elemId : *geo.pActiveElements) {
            const auto& conn = isVolume ? m_mesh.getCellConnectivity(elemId) : m_mesh.getFacetConnectivity(elemId);
            for (auto gId : conn) {
                if (geo.globalToLocal[gId] == -1) {
                    geo.globalToLocal[gId] = geo.localToGlobal.size();
                    geo.localToGlobal.push_back(gId);
                }
                geo.connectivity.push_back(geo.globalToLocal[gId]);
            }
            runningOffset += conn.size();
            geo.offsets.push_back(runningOffset);
            
            auto elemType = isVolume ? m_mesh.getCell(elemId).type : m_mesh.getFacet(elemId).type;
            geo.types.push_back((unsigned char)getVTKCellType(elemType));
        }

        geo.numActiveVerts = geo.localToGlobal.size();
        geo.points.resize(3 * geo.numActiveVerts);
        for (size_t i = 0; i < geo.numActiveVerts; ++i) {
            const auto& c = m_mesh.getVertex(geo.localToGlobal[i]).coords;
            geo.points[3*i+0] = c[0]; geo.points[3*i+1] = c[1]; geo.points[3*i+2] = c[2];
        }

        return geo;
    }

    // ========================================================================
    // --- KERNEL 2: GEOMETRY EXTRACTION (DISCONTINUOUS)
    // ========================================================================
    VTKFamilyExporter::ExtractedGeometry VTKFamilyExporter::buildDiscontinuousGeometry(CfemInt entityId) const
    {
        ExtractedGeometry geo;
        const MeshEntity& entity = m_mesh.getEntityConst(entityId);
        bool isVolume = (entity.getDimension() == m_mesh.getDimension());
        
        geo.pActiveElements = isVolume ? &entity.getCellIds() : &entity.getFacetIds();
        geo.numActiveCells = geo.pActiveElements->size();

        size_t totalDgPoints = 0;
        for (CfemInt elemId : *geo.pActiveElements) {
            totalDgPoints += isVolume ? m_mesh.getCellConnectivity(elemId).size() : m_mesh.getFacetConnectivity(elemId).size();
        }

        geo.numActiveVerts = totalDgPoints;
        geo.points.resize(3 * totalDgPoints);
        geo.connectivity.resize(totalDgPoints);
        geo.localToGlobal.resize(totalDgPoints); // We map each duplicated node back to its global ID

        int ptIdx = 0;
        for (CfemInt elemId : *geo.pActiveElements) {
            const auto& conn = isVolume ? m_mesh.getCellConnectivity(elemId) : m_mesh.getFacetConnectivity(elemId);
            for (auto gId : conn) {
                const auto& c = m_mesh.getVertex(gId).coords;
                geo.points[3*ptIdx+0] = c[0]; geo.points[3*ptIdx+1] = c[1]; geo.points[3*ptIdx+2] = c[2];
                geo.connectivity[ptIdx] = ptIdx;
                geo.localToGlobal[ptIdx] = gId; 
                ptIdx++;
            }
            geo.offsets.push_back(ptIdx);
            auto elemType = isVolume ? m_mesh.getCell(elemId).type : m_mesh.getFacet(elemId).type;
            geo.types.push_back((unsigned char)getVTKCellType(elemType));
        }

        return geo;
    }

    VTKFamilyExporter::EvaluatedData VTKFamilyExporter::evaluateFields(const ExtractedGeometry& geo, 
                                                                       const std::vector<ExportEntry>& entries, 
                                                                       CfemInt entityId, 
                                                                       bool isDiscontinuous) const
    {
        EvaluatedData data;
        const MeshEntity& entity = m_mesh.getEntityConst(entityId);
        bool isVolume = (entity.getDimension() == m_mesh.getDimension());
        bool isFacet  = (entity.getDimension() == m_mesh.getDimension() - 1);

        size_t numEntries = entries.size();
        ShapeEvaluationFlags fieldFlags = ShapeEvaluationFlags::None;
        if (isDiscontinuous && isFacet) fieldFlags |= ShapeEvaluationFlags::Normal;

        data.fields.resize(numEntries);
        for (size_t i = 0; i < numEntries; ++i) {
            entries[i].expr->prepare();
            fieldFlags |= entries[i].expr->requiredShapeFlags();
            data.fields[i].resize(geo.numActiveVerts * entries[i].expr->getDimension(), 0.0);
        }

        // FIX: O(1) Cache using unordered_map
        struct TypeCache { std::vector<ShapeInfo> nodeShapes; const ReferenceElement *refEl; };
        std::unordered_map<CellType, TypeCache> geoCache; 
        CfemInt geoOrder = m_mesh.getOrderInt();

        auto getOrCacheRefEl = [&](CellType type) -> TypeCache& {
            auto it = geoCache.find(type);
            if (it == geoCache.end()) {
                auto &cache = geoCache[type];
                cache.refEl = &ReferenceElementFactory::getReferenceElement(type, geoOrder);
                for (const auto &xi : cache.refEl->referenceNodes()) {
                    ShapeInfo info; cache.refEl->evaluate(xi, ShapeEvaluationFlags::Value | fieldFlags, info);
                    info.shapeMustBeUpdated = false;
                    cache.nodeShapes.push_back(info);
                }
                return cache;
            }
            return it->second;
        };

        EvalContext ctx;
        MappingInfo scratchFacetMapping;
        DynamicVector<Vector3D> physCoordsElem, physCoordsParent;
        std::vector<uint8_t> computed(geo.numActiveVerts, 0);

        int ptIdx = 0;
        for (CfemInt elemId : *geo.pActiveElements) {
            CfemInt parentCellId = elemId;
            const std::vector<CfemInt>* localFaceConnectivity = nullptr;

            if (isFacet) {
                const auto& face = m_mesh.getFacet(elemId);
                parentCellId = face.left_cell;
                localFaceConnectivity = &Cell::getLocalFacesIndices(m_mesh.getCell(parentCellId).type)[face.left_cell_facet_id];
            }

            auto elemType = isVolume ? m_mesh.getCell(elemId).type : m_mesh.getFacet(elemId).type;
            auto parentCellType = m_mesh.getCell(parentCellId).type;

            auto &elemCache = getOrCacheRefEl(elemType);         
            auto &parentCache = getOrCacheRefEl(parentCellType); 

            if (isFacet) m_mesh.getFacetVerticesCoords(elemId, physCoordsElem);
            m_mesh.getCellVerticesCoords(parentCellId, physCoordsParent);

            const auto &conn = isVolume ? m_mesh.getCellConnectivity(elemId) : m_mesh.getFacetConnectivity(elemId);
            
            for (size_t n = 0; n < conn.size(); ++n) {
                // FIX: Unified index resolution. O(1) for both Continuous and DG
                CfemInt evalIdx = isDiscontinuous ? ptIdx : geo.globalToLocal[conn[n]];
                
                if (!isDiscontinuous && computed[evalIdx]) continue;

                CfemInt parentNodeIdx = isVolume ? n : (*localFaceConnectivity)[n];
                ctx.resetBuffers(); 
                MappingInfo& mapParent = ctx.scratchMapping;

                parentCache.refEl->computeMapping(parentCache.refEl->referenceNodes()[parentNodeIdx], physCoordsParent,
                                           ShapeEvaluationFlags::PhysicalPosition | fieldFlags,
                                           parentCache.nodeShapes[parentNodeIdx], mapParent);

                if (isFacet && hasFlag(fieldFlags, ShapeEvaluationFlags::Normal)) {
                    elemCache.refEl->computeMapping(elemCache.refEl->referenceNodes()[n], physCoordsElem,
                                             ShapeEvaluationFlags::Normal, elemCache.nodeShapes[n], scratchFacetMapping);
                    mapParent.normal = scratchFacetMapping.normal;
                }

                ctx.physPos = mapParent.physicalPosition;
                ctx.geometricMapping = &mapParent;

                for (size_t i = 0; i < numEntries; ++i) {
                    int dim = entries[i].expr->getDimension();
                    std::span<CfemReal> outSpan(&data.fields[i][evalIdx * dim], dim);
                    entries[i].expr->evaluate(parentCellId, parentCache.refEl->referenceNodes()[parentNodeIdx], outSpan, ctx);
                }
                
                if (!isDiscontinuous) computed[evalIdx] = 1;
                if (isDiscontinuous) ptIdx++;
            }
        }

        for (const auto &name : m_exportedEntityNames) {
            
            if ( !m_mesh.isPartOfTheMeshEntities(name) ) continue;

            const MeshEntity& ent = m_mesh.getEntityConst(name);

            std::vector<uint8_t> isInEntity(m_mesh.getNumVertices(), 0);
            for (CfemInt vId : ent.getVertexIds()) isInEntity[vId] = 1;

            std::vector<float> mask(geo.numActiveVerts, 0.0f);
            for (size_t i = 0; i < geo.numActiveVerts; ++i) {
                if (isInEntity[geo.localToGlobal[i]]) mask[i] = 1.0f;
            }
            data.nodeMasks.push_back(std::move(mask));
        }

        if (m_exportQuality) {
            data.cellQuality.reserve(geo.numActiveCells);
            for (CfemInt elemId : *geo.pActiveElements) {
                data.cellQuality.push_back(isVolume ? m_mesh.computeCellQuality(elemId) : m_mesh.computeFacetQuality(elemId));
            }
        }
        
        for (const auto &name : m_exportedEntityNames) {
            if ( !m_mesh.isPartOfTheMeshEntities(name) ) continue;
            
            const MeshEntity& ent = m_mesh.getEntityConst(name);

            std::vector<uint8_t> isEntityCell(m_mesh.getNumCells(), 0);
            std::vector<uint8_t> isEntityFacet(m_mesh.getNumFacets(), 0);
            
            if (isVolume) for (CfemInt cId : ent.getCellIds()) isEntityCell[cId] = 1;
            else          for (CfemInt fId : ent.getFacetIds()) isEntityFacet[fId] = 1;

            std::vector<float> mask(geo.numActiveCells, 0.0f);
            int idx = 0;
            for (CfemInt elemId : *geo.pActiveElements) {
                if ((isVolume && isEntityCell[elemId]) || (!isVolume && isEntityFacet[elemId])) mask[idx] = 1.0f;
                idx++;
            }
            data.cellMasks.push_back(std::move(mask));
        }

        return data;
    }
}