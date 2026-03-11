#ifndef CFEM_REFERENCE_ELEMENT_FACTORY_H
#define CFEM_REFERENCE_ELEMENT_FACTORY_H

#include "ReferenceElement.h"
#include "cfem_topology.h" // Pour CellType
#include <map>
#include <memory>

namespace cfem {

class ReferenceElementFactory {
public:
    /**
     * @brief Retourne l'élément de référence correspondant au type de cellule.
     * Cette instance contient les coordonnées xi de TOUS les noeuds du maillage.
     */
    static const ReferenceElement& getReferenceElement(CellType type);

private:
    // On garde les instances en mémoire pour éviter de les recréer
    static std::map<CellType, std::unique_ptr<ReferenceElement>> m_instances;
};




} // namespace cfem

#endif