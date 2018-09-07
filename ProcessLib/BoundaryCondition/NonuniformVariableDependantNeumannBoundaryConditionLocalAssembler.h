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
#include "ProcessLib/Parameter/MeshNodeParameter.h"

#include "GenericNonuniformNaturalBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
struct NonuniformVariableDependantNeumannBoundaryConditionData
{
    /* TODO use Parameter */
    MeshLib::PropertyVector<double> const& constant;
    MeshLib::PropertyVector<double> const& prefac1;
    MeshLib::PropertyVector<double> const& prefac2;
    // Used for mapping boundary nodes to bulk nodes.
    std::size_t bulk_mesh_id;
    MeshLib::PropertyVector<std::size_t> const& mapping_to_bulk_nodes;
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk;
    int const variable_id_bulk;
    int const component_id_bulk;
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class NonuniformVariableDependantNeumannBoundaryConditionLocalAssembler final
    : public GenericNonuniformNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNonuniformNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;
    using NodalVectorType = typename Base::NodalVectorType;
    using NodalMatrixType = typename Base::NodalMatrixType;

public:
    /// The neumann_bc_term factor is directly integrated into the local
    /// element matrix.
    NonuniformVariableDependantNeumannBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const local_matrix_size,
        bool const is_axially_symmetric,
        unsigned const integration_order,
        NonuniformVariableDependantNeumannBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          _data(data),
          _local_rhs(local_matrix_size)
    {
    }

    void assemble(std::size_t const id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& x,
                  GlobalMatrix& /*K*/, GlobalVector& b,
                  GlobalMatrix* /*Jac*/) override
    {
        _local_rhs.setZero();

        MeshNodeParameter<double> constant_values{"ConstantValues", _data.constant};
        MeshNodeParameter<double> prefac1_values{"Prefac1Values", _data.prefac1};
        MeshNodeParameter<double> prefac2_values{"Prefac2Values", _data.prefac2};
        // Get element nodes for the interpolation from nodes to
        // integration point.
        NodalVectorType constant_node_values =
            constant_values.getNodalValuesOnElement(Base::_element, t);
        NodalVectorType prefac1_node_values =
            prefac1_values.getNodalValuesOnElement(Base::_element, t);
        NodalVectorType prefac2_node_values =
            prefac2_values.getNodalValuesOnElement(Base::_element, t);
        auto const num_nodes = x.size();
        int offset; double localv1; double localv2;
        if(id == 0){
            offset=1;}
        else{ offset=-1;}
//        std::cout << "num_nodes: " << num_nodes << std::endl;
        unsigned const n_integration_points =
            Base::_integration_method.getNumberOfPoints();
        NodalVectorType neumann_node_values;
        auto const indices = NumLib::getIndices(id, dof_table_boundary);  
        auto const local_x = x.get(indices);

       

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            localv1 = x[indices[ip]];
            localv2 = x[indices[ip]+offset];
//            std::cout << "p = " << localv1 << " c = " << localv2 << std::endl;
                    

            auto const& n_and_weight = Base::_ns_and_weights[ip];
            auto const& N = n_and_weight.N;
            auto const& w = n_and_weight.weight;
            neumann_node_values = constant_node_values + 
                    prefac1_node_values*localv1+
                    prefac2_node_values*localv2;
//            std::cout << neumann_node_values.dot(N) << "dot product" << std::endl;
//            std::cout << neumann_node_values[ip] << std::endl; 
            _local_rhs.noalias() += N * neumann_node_values.dot(N) * w;
        }

        for(unsigned int i= 0 ; i<indices.size(); i++){
//            std::cout << indices[i] << std::endl;}
//        std::cout << "__________________" << std::endl;
        b.add(indices, _local_rhs);}
//        std::cout << "__________________" << std::endl;
    }

private:
    NonuniformVariableDependantNeumannBoundaryConditionData const& _data;
    NodalVectorType _local_rhs;

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace ProcessLib
