/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "GenericNaturalBoundaryConditionLocalAssembler.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/Process.h"
#include "MeshLib/Elements/MapBulkElementPoint.h"

namespace ProcessLib
{
struct TransportReflectionBoundaryConditionData final {
    ParameterLib::Parameter<double> const& boundary_permeability;
    MeshLib::PropertyVector<std::size_t> const bulk_face_ids;
    MeshLib::PropertyVector<std::size_t> const bulk_element_ids;
    Process const& process;
};


template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>

class TransportReflectionBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;

protected:
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;

    struct NAndWeight
    {
        NAndWeight(typename ShapeMatricesType::ShapeMatrices::ShapeType&& N_,
                   typename ShapeMatricesType::ShapeMatrices::ShapeType&& dNdx_,
                   double const weight_)
            : N(std::move(N_)), dNdx(dNdx_), weight(weight_)
        {
        }
        typename ShapeMatricesType::ShapeMatrices::ShapeType const N;
        typename ShapeMatricesType::ShapeMatrices::ShapeType const dNdx;
        double const weight;
    };

private:
    static std::vector<NAndWeight, Eigen::aligned_allocator<NAndWeight>>
    initNsAndWeights(MeshLib::Element const& e, bool is_axially_symmetric,
                     unsigned const integration_order)
    {
        IntegrationMethod const integration_method(integration_order);
        std::vector<NAndWeight, Eigen::aligned_allocator<NAndWeight>>
            ns_and_weights;
        ns_and_weights.reserve(integration_method.getNumberOfPoints());

        auto sms = initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                     IntegrationMethod, GlobalDim>(
            e, is_axially_symmetric, integration_method);
        for (unsigned ip = 0; ip < sms.size(); ++ip)
        {
            auto& sm = sms[ip];
            double const w =
                sm.detJ * sm.integralMeasure *
                integration_method.getWeightedPoint(ip).getWeight();

            ns_and_weights.emplace_back(std::move(sm.N), std::move(sm.dNdx),w);
        }

        return ns_and_weights;
    }    
    
    
public:
    TransportReflectionBoundaryConditionLocalAssembler(MeshLib::Element const& e,
                                         std::size_t const local_matrix_size,
                                         bool is_axially_symmetric,
                                         unsigned const integration_order,
                                         TransportReflectionBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          _data(data),
          _local_rhs(local_matrix_size)
    {
    }

    // TODO also implement derivative for Jacobian in Newton scheme.
    void assemble(std::size_t const mesh_item_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& x, GlobalMatrix& /*K*/,
                  GlobalVector& b, GlobalMatrix* /*Jac*/) override
    {
        _local_rhs.setZero();

        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();
        auto const indices =
                    NumLib::getIndices(mesh_item_id, dof_table_boundary);
        std::vector<double> local_x =
            x.get(indices);
        NodalVectorType const current_variable_nodal_values =
            Eigen::Map<const NodalVectorType>(local_x.data(), local_x.size());
        std::size_t bulk_element_id = _data.bulk_element_ids[Base::_element.getID()];
        std::size_t bulk_face_id = _data.bulk_face_ids[Base::_element.getID()];
        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& ip_data = _ns_and_weights[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.weight;
            auto const& wp = Base::_integration_method.getWeightedPoint(ip);

            auto const bulk_element_point = MeshLib::getBulkElementPoint(
                _data.process.getMesh(), bulk_element_id, bulk_face_id, wp);

            // adding a alpha term to the diagonal of the stiffness matrix
            // and a alpha * u_0 term to the rhs vector
            _local_rhs.noalias() += w * N *  _data.process.getFlux(bulk_element_id, bulk_element_point, t, x).dot(dNdx.transpose() * current_variable_nodal_values);
            
        }

        b.add(indices, _local_rhs);
    }

private:
    TransportReflectionBoundaryConditionData const& _data;
    std::vector<NAndWeight, Eigen::aligned_allocator<NAndWeight>> const
        _ns_and_weights;
    typename Base::NodalVectorType _local_rhs;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib