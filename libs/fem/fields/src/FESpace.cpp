#include "FESpace.h"
#include "FEField.h"

namespace cfem
{

    void FESpace::registerField(std::shared_ptr<FEField> field)
    {
        NamingPolicy::validate(field->getName());

        if (m_mesh->nameAlreadyUsed(field->getName()))
        {
            CFEM_ERROR << " Global naming collision: Field with name '" << field->getName() << "' already exists.";
        }
        
        m_fields.push_back(field);
        m_mesh->registerFieldName(field->getName());
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