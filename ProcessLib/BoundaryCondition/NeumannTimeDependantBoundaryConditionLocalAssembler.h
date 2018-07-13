/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/PropertyVector.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Function/Interpolation.h"

#include "GenericNonuniformNaturalBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
struct NeumannTimeDependantBoundaryConditionData
{
    /* TODO use Parameter */

    // Used for mapping boundary nodes to bulk nodes.
    std::size_t bulk_mesh_id;
    MeshLib::PropertyVector<std::size_t> const& mapping_to_bulk_nodes;
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk;
    int const variable_id_bulk;
    int const component_id_bulk;
    MeshLib::PropertyVector<double> const& gravitation_values;
    MeshLib::PropertyVector<double> const& psi_leaf;
    MeshLib::PropertyVector<double> const& resistance;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class NeumannTimeDependantBoundaryConditionLocalAssembler final
    : public GenericNonuniformNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNonuniformNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;

public:
    /// The neumann_bc_value factor is directly integrated into the local
    /// element matrix.
    NeumannTimeDependantBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool is_axially_symmetric,
        unsigned const integration_order,
        NeumannTimeDependantBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          _data(data),
          _local_rhs(local_matrix_size)
    {
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& x,
                  GlobalMatrix& /*K*/, GlobalVector& b) override
    {
        // concentration has id 0, pressure has id 1
        _local_rhs.setZero();
        auto indices = NumLib::getIndices(id, dof_table_boundary);
        std::cout << t << std::endl;
//        for(unsigned ip = 0; ip < n_integration_pointsx; ip++)
//        {
//            std::cout << x.get(ip) << "bei der id " << id << std::endl;
//        }
        // TODO lots of heap allocations.
        std::vector<double> neumann_param_nodal_values_local;
        neumann_param_nodal_values_local.reserve(indices.size());
        for (auto i : indices)
        {
            double resistance = _data.resistance.getComponent(i, 0);
            double local_neumann_gradient = resistance*_data.psi_leaf.getComponent(i, 0);
            local_neumann_gradient += resistance*_data.gravitation_values.getComponent(i, 0);

            auto const bulk_node_id =
                _data.mapping_to_bulk_nodes.getComponent(i, 0);
            
            double concentration = x[2*bulk_node_id];
            double pressure = x[2*bulk_node_id+1];
            
            local_neumann_gradient += resistance*concentration * 85000;
            local_neumann_gradient -= resistance*pressure;
            if(_data.variable_id_bulk == 0 ){local_neumann_gradient*=-concentration;}      
            std::cout << "concentration at node = " << concentration << std::endl;
            std::cout << "pressure  at node = " << x[2*bulk_node_id+1] <<std::endl;
            std::cout << "local neumann gradient = " << local_neumann_gradient <<std::endl;
            std::cout << "local neumann gradient without resistance = " << local_neumann_gradient/resistance <<std::endl;
            std::cout << "variable_id_bulk = " << _data.variable_id_bulk << std::endl;
            
            int ti(t);
            ti/= 43200;
            ti=ti%2;
            std::cout << "t_i = " << ti << "and t = " << t << std::endl;
            
            neumann_param_nodal_values_local.push_back(
                ti*local_neumann_gradient);
            std::cout <<" ---" << std::endl;
            
        }

        auto const n_integration_points = Base::_ns_and_weights.size();
        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            auto const& n_and_weight = Base::_ns_and_weights[ip];
            double flux;
            NumLib::shapeFunctionInterpolate(neumann_param_nodal_values_local,
                                             n_and_weight.N, flux);
            _local_rhs.noalias() +=
                n_and_weight.N * (n_and_weight.weight * flux);
        }

        // map boundary dof indices to bulk dof indices
        for (auto& i : indices)
        {
            std::cout << i << std::endl;
            auto const bulk_node_id =
                _data.mapping_to_bulk_nodes.getComponent(i, 0);
            std::cout << bulk_node_id << std::endl;

            std::cout << "variable id " << _data.variable_id_bulk << std::endl;

            MeshLib::Location const l{
                _data.bulk_mesh_id, MeshLib::MeshItemType::Node, bulk_node_id};

            i = _data.dof_table_bulk.getGlobalIndex(l, _data.variable_id_bulk,
                                                    _data.component_id_bulk);
            assert(i != NumLib::MeshComponentMap::nop);
        }
        b.add(indices, _local_rhs);
    }

private:
    NeumannTimeDependantBoundaryConditionData const& _data;
    typename Base::NodalVectorType _local_rhs;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib
