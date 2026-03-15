#pragma once

#include "cfemutils.h"
#include <limits>

namespace cfem
{
    class Mesh;

struct QualityReport {
    CfemReal minVolume = std::numeric_limits<CfemReal>::max();
    CfemReal maxVolume = -std::numeric_limits<CfemReal>::max();
    int invertedCount = 0;
    int degenerateCount = 0;
};

class MeshQualityTool {
public:
    static QualityReport analyze(const Mesh& mesh);
    static void printSummary(const QualityReport& report);
};

}