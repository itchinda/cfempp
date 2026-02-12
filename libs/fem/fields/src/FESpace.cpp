#include "FESpace.h"
#include "FEField.h"

namespace cfem
{

    void FESpace::registerField(std::shared_ptr<FEField> field)
    {
        m_fields.push_back(field);
    }

    /**
     * @brief Point d'entrée pour la cascade de mise à jour.
     */
    void FESpace::updateAfterMeshChange() {
        this->validate();
        this->enumerateDofs();

        // Notifier les champs (température, pression, etc.)
        for (auto it = m_fields.begin(); it != m_fields.end(); ) {
            if (auto field = it->lock()) {
                field->reallocate(); // Redimensionner les vecteurs PETSc
                ++it;
            } else {
                it = m_fields.erase(it);
            }
        }
    }

} // namespace