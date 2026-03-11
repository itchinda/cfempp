#include "symengine/SymEngineSymbols.h"

namespace cfem::sym
{
    inline std::vector<RCP<const Symbol>> extract_symbols(const Expression& ex)
    {
        auto s = SymEngine::free_symbols(*ex.get_basic()); 
        std::vector<RCP<const Symbol>> symbols;
        symbols.reserve(s.size());
        for (auto &b : s) {
            if (SymEngine::is_a<Symbol>(*b)) {
                symbols.push_back(SymEngine::rcp_static_cast<const Symbol>(b));
            }
        }
        return symbols;
    }

}