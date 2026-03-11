#pragma once

#include "FEField.h"

namespace cfem
{

    template <typename T>
    Vector3D FEField::fetchValue(
                            const T &obj, 
                            CfemInt cellId,
                            const Vector3D &xi,
                            const Vector3D &physPos,
                            std::span<CfemReal> buffer
                        )
    {
        if constexpr (std::is_convertible_v<T, const FEField &>)
        {
            obj.evaluateAt(cellId, xi, buffer);
            if (buffer.size() == 1)
                return Vector3D{buffer[0], 0.0, 0.0};
            if (buffer.size() == 2)
                return Vector3D{buffer[0], buffer[1], 0.0};
            return Vector3D{buffer[0], buffer[1], buffer[2]};
        }
        else
        {
            auto result = obj(physPos);
            if constexpr (std::is_arithmetic_v<decltype(result)>)
                return Vector3D{static_cast<CfemReal>(result), 0.0, 0.0};
            else
                return result;
        }
    }

    // template <typename Func>
    //     requires(!std::is_base_of_v<sym::CfemExpression, std::decay_t<Func>> &&
    //              !std::is_convertible_v<Func, std::shared_ptr<sym::CfemExpression>>)
    template <typename Func>
    requires (!std::is_same_v<std::decay_t<Func>, sym::CfemExpression> &&
            !std::is_same_v<std::decay_t<Func>, std::shared_ptr<sym::CfemExpression>>)
    void FEField::project(const Func& func)
    {
        auto mesh = m_space->getMesh();
        auto &vec = getData();
        CfemInt N = m_space->getNumActiveNodes();
        int dim = getDimension();
        CfemInt totalVertices = mesh->getNumVertices();

        for (CfemInt i = 0; i < totalVertices; ++i)
        {
            CfemInt nodeIdx = m_space->getLocalNodeIdx(i);
            if (nodeIdx == -1)
                continue;
            const auto &vertex = mesh->getVertex(i);
            auto result = func(vertex.coords);

            if constexpr (std::is_same_v<decltype(result), Vector3D>)
            {
                for (int d = 0; d < dim; ++d)
                    vec.setValue(nodeIdx + d * N, result[d], INSERT_VALUES);
            }
            else
            {
                vec.setValue(nodeIdx, static_cast<CfemReal>(result), INSERT_VALUES);
            }
        }
        vec.assemble();
    }

}
