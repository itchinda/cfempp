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