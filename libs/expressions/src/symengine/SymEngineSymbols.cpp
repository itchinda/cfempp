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

#include "symengine/SymEngineSymbols.h"

namespace cfem::sym
{
    inline std::vector<SymEngine::RCP<const SymEngine::Symbol>> extract_symbols(const SymEngine::Expression& ex)
    {
        auto s = SymEngine::free_symbols(*ex.get_basic()); 
        std::vector<SymEngine::RCP<const SymEngine::Symbol>> symbols;
        symbols.reserve(s.size());
        for (auto &b : s) {
            if (SymEngine::is_a<SymEngine::Symbol>(*b)) {
                symbols.push_back(SymEngine::rcp_static_cast<const SymEngine::Symbol>(b));
            }
        }
        return symbols;
    }

}