#include "MeshQualityTool.h"
#include "Mesh.h"

namespace cfem
{
    QualityReport MeshQualityTool::analyze(const Mesh& mesh) {
        QualityReport report;

        for (size_t i = 0; i < mesh.getNumCells(); ++i) {
            const auto& cell = mesh.getCell(i);
            
            // Signed Volume (Jacobian)
            CfemReal vol = mesh.computeCellMeasure(cell.type, mesh.getCellConnectivity(cell.id) );
            
            if (vol < -cfem::ZERO_TOLERANCE) {
                report.invertedCount++;
            } else if (std::abs(vol) <= cfem::ZERO_TOLERANCE) {
                report.degenerateCount++;
            }

            report.minVolume = std::min(report.minVolume, vol);
            report.maxVolume = std::max(report.maxVolume, vol);
        }

        return report;
    }

    void MeshQualityTool::printSummary(const QualityReport& report) {
        std::cout << "\n--- CFEM MESH QUALITY REPORT ---" << std::endl;
        std::cout << "Min Volume: " << report.minVolume << std::endl;
        std::cout << "Max Volume: " << report.maxVolume << std::endl;
        
        if (report.invertedCount > 0) {
            CFEM_WARN << "Mesh contains " << std::to_string(report.invertedCount) << " inverted elements!";
        }
        if (report.degenerateCount > 0) {
            CFEM_WARN << "Mesh contains " << std::to_string(report.degenerateCount) << " degenerate (zero-volume) elements!";
        }
        
        if (report.invertedCount == 0 && report.degenerateCount == 0) {
            std::cout << "Result: Mesh is valid for simulation." << std::endl;
        }
        std::cout << "--------------------------------\n" << std::endl;
    }

    
}