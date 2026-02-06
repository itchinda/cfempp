#ifndef CFEM_REFERENCE_ELEMENT_
#define CFEM_REFERENCE_ELEMENT_

#include "cfemutils.h"
#include "Vector3D.h"
#include "DynamicVector.h"
#include "DynamicMatrix.h"
#include "Tensor.h"

namespace cfem {

/**
 * @struct ShapeInfo
 * @brief Buffer de stockage pour les évaluations des fonctions de forme.
 * @note Réutiliser cet objet pour tous les éléments d'un même type afin d'éviter les allocations.
 */
struct ShapeInfo {
    DynamicVector<CfemReal> values;                // Stocke les N_i
    DynamicVector<Vector3D> gradients;             // Stocke les dN_i/dxi
    DynamicVector<SymTensorO2> hessians; // Stocke les d2N_i/dxi2

    CfemIntFlag validFlags = cfem::ShapeActiveField::None;

    /** @brief Prépare la taille du buffer pour éviter les réallocations. */
    void setup(int numNodes) {
        values.resize(numNodes);
        gradients.resize(numNodes);
        // On n'alloue les hessiennes que si nécessaire
        validFlags = cfem::ShapeActiveField::None;
    }
    
    bool has(cfem::ShapeActiveField field) const {
        return (validFlags & field) != 0;
    }
};

struct MappingInfo {
    TensorO2 jacobian;  
    TensorO2 invJacobian;
    CfemReal detJ;
    DynamicVector<Vector3D> physicalShapeGradients;
    DynamicVector<SymTensorO2> physicalShapeHessians;

    inline CfemReal dNdx(int nodeIdx) const { return physicalShapeGradients[nodeIdx].x; }
    inline CfemReal dNdy(int nodeIdx) const { return physicalShapeGradients[nodeIdx].y; }
    inline CfemReal dNdz(int nodeIdx) const { return physicalShapeGradients[nodeIdx].z; }

    inline CfemReal d2Ndxdx(int nodeIdx) const { return physicalShapeHessians[nodeIdx].xx(); }
    inline CfemReal d2Ndydy(int nodeIdx) const { return physicalShapeHessians[nodeIdx].yy(); }
    inline CfemReal d2Ndzdz(int nodeIdx) const { return physicalShapeHessians[nodeIdx].zz(); }
    inline CfemReal d2Ndydz(int nodeIdx) const { return physicalShapeHessians[nodeIdx].yz(); }
    inline CfemReal d2Ndxdz(int nodeIdx) const { return physicalShapeHessians[nodeIdx].xz(); }
    inline CfemReal d2Ndxdy(int nodeIdx) const { return physicalShapeHessians[nodeIdx].xy(); }

};


/**
 * @class ReferenceElement
 * @brief Interface pour les éléments finis de référence (coordonnées locales xi).
 */
class ReferenceElement {
protected:
    DynamicVector<Vector3D> m_nodes;

public:
    virtual ~ReferenceElement() = default;

    virtual int dimension() const = 0;
    virtual int numNodes() const = 0;

    /**
     * @brief Évalue les fonctions de forme au point xi.
     * @param xi Coordonnées locales dans l'élément parent.
     * @param requested Masque des champs à calculer.
     * @param out Buffer ShapeInfo à remplir.
     */
    virtual void evaluate(const Vector3D& xi, 
                          CfemIntFlag requested, 
                          ShapeInfo& out) const = 0;

    const DynamicVector<Vector3D>& referenceNodes() const { return m_nodes; }


    /**
     * @brief Computes the transformation from reference to physical space
     * @param xi Point in reference space
     * @param physicalNodesCoords Physical coordinates of the element's nodes
     * @param out Calculated Jacobian and metrics
     */

    void computeMapping(const Vector3D& xi, 
                        const DynamicVector<Vector3D>& physicalNodesCoords, 
                        CfemIntFlag requestedFields, // Utilisation du bitmask
                        MappingInfo& out) const {
        
        ShapeInfo shape;
        shape.setup(numNodes());
        
        // On propage les flags à l'évaluation des fonctions de forme
        // On s'assure que les gradients sont TOUJOURS calculés, car le Jacobien en dépend.
        CfemIntFlag internalFields = requestedFields | ShapeActiveField::Gradient;
        evaluate(xi, internalFields, shape);
        
        // 1. Calcul du Jacobien (toujours nécessaire pour le mapping de base)
        computeJacobian(shape, physicalNodesCoords, out); 
        
        // 2. Gradients Physiques (si demandés)
        if (requestedFields & ShapeActiveField::Gradient) {
            computePhysicalGradients(shape, out);
        }

        // 3. Hessiennes Physiques (si demandées)
        if (requestedFields & ShapeActiveField::Hessian) {
            computePhysicalHessians(physicalNodesCoords, shape, out);
        }
    }

    void computeJacobian(const ShapeInfo& shape, 
                        const DynamicVector<Vector3D>& physicalNodesCoords, 
                        MappingInfo& out) const {
        
        out.jacobian.setZero();
        for (int i = 0; i < numNodes(); ++i) {
            // J_00 = sum( x_i * dNi/dxi )
            out.jacobian(0,0) += physicalNodesCoords[i].x * shape.gradients[i].x;
            out.jacobian(0,1) += physicalNodesCoords[i].x * shape.gradients[i].y;
            out.jacobian(0,2) += physicalNodesCoords[i].x * shape.gradients[i].z;
            
            out.jacobian(1,0) += physicalNodesCoords[i].y * shape.gradients[i].x;
            out.jacobian(1,1) += physicalNodesCoords[i].y * shape.gradients[i].y;
            out.jacobian(1,2) += physicalNodesCoords[i].y * shape.gradients[i].z;
            
            out.jacobian(2,0) += physicalNodesCoords[i].z * shape.gradients[i].x;
            out.jacobian(2,1) += physicalNodesCoords[i].z * shape.gradients[i].y;
            out.jacobian(2,2) += physicalNodesCoords[i].z * shape.gradients[i].z;
        }

        if (dimension() == 2) {
            out.jacobian(2,2) = 1.0; 
        }
        
        out.detJ = det(out.jacobian);
        out.invJacobian = inverse(out.jacobian, out.detJ);
    }

    void computePhysicalGradients(const ShapeInfo& shape, 
                                  MappingInfo& out) const {
        
        // Transform gradients: Grad_xyz = Grad_xi * J^-1
        out.physicalShapeGradients.resize(numNodes());
        
        for (int a = 0; a < out.physicalShapeGradients.size(); ++a) {
            const Vector3D& g = shape.gradients[a];
            out.physicalShapeGradients[a].x = out.invJacobian(0,0)*g.x + out.invJacobian(1,0)*g.y + out.invJacobian(2,0)*g.z;
            out.physicalShapeGradients[a].y = out.invJacobian(0,1)*g.x + out.invJacobian(1,1)*g.y + out.invJacobian(2,1)*g.z;
            out.physicalShapeGradients[a].z = out.invJacobian(0,2)*g.x + out.invJacobian(1,2)*g.y + out.invJacobian(2,2)*g.z;
        }
    }

    /**
     * @brief Calcule les Hessiennes physiques pour chaque noeud.
     * @param physicalNodesCoords Coordonnées réelles des noeuds de l'élément.
     * @param shape Contient les gradients et hessiennes de référence.
     * @param out MappingInfo à remplir (doit avoir jacobian et invJacobian à jour).
     */
    void computePhysicalHessians(const DynamicVector<Vector3D>& physicalNodesCoords, 
                                const ShapeInfo& shape, 
                                MappingInfo& out) const {
        
        int nNodes = numNodes();
        out.physicalShapeHessians.resize(nNodes);

        // 1. Calculer la courbure de la transformation : H_xi(x_k)
        // C'est la dérivée seconde de la géométrie par rapport à l'espace de référence.
        SymTensorO2 H_xi_X, H_xi_Y, H_xi_Z;
        H_xi_X.setZero(); H_xi_Y.setZero(); H_xi_Z.setZero();

        for (int a = 0; a < nNodes; ++a) {
            const SymTensorO2& ha = shape.hessians[a];
            H_xi_X += ha * physicalNodesCoords[a].x;
            H_xi_Y += ha * physicalNodesCoords[a].y;
            H_xi_Z += ha * physicalNodesCoords[a].z;
        }

        // 2. Transformer chaque Hessienne de noeud
        for (int a = 0; a < nNodes; ++a) {
            // Terme correctif dû à la courbure : R = H_xi(Na) - sum( dNa/dx_k * H_xi(x_k) )
            // Note : On utilise les gradients physiques dN/dx déjà calculés
            SymTensorO2 R = shape.hessians[a];
            R -= H_xi_X * out.dNdx(a);
            R -= H_xi_Y * out.dNdy(a);
            R -= H_xi_Z * out.dNdz(a);

            // Transformation finale : H_phys = J^-T * R * J^-1
            // On peut le faire via : out.physicalShapeHessians[a] = transformHessian(R, out.invJacobian);
            out.physicalShapeHessians[a] = R.quadform (out.invJacobian); 
        }
    }
};

}

#endif