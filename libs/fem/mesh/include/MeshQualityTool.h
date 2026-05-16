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

#pragma once

#include "cfemutils.h"
#include <limits>

namespace cfem
{
    class Mesh;

struct QualityReport {
    CfemReal minVolume = std::numeric_limits<CfemReal>::max();
    CfemReal maxVolume = -std::numeric_limits<CfemReal>::max();
    CfemInt invertedCount = 0;
    CfemInt degenerateCount = 0;
};

class MeshQualityTool {
public:
    static QualityReport analyze(const Mesh& mesh);
    static void printSummary(const QualityReport& report);
};

}