/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "GenericNaturalBoundaryCondition.h"
#include "TransportReflectionBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
using TransportReflectionBoundaryCondition = GenericNaturalBoundaryCondition<
    TransportReflectionBoundaryConditionData,
    TransportReflectionBoundaryConditionLocalAssembler>;

/*! Creates a new uniform TransportReflection boundary condition from the given data.
 *
 * The TransportReflection boundary condition is given in the form
 * \f$ \alpha \cdot [ u_0 - u(x) ] \f$,
 * where the coefficients \f$ \alpha \f$ and \f$ u_0 \f$ are obtained from the
 * \c config, and \f$ u \f$ is the unknown to which the boundary condition is
 * applied.
 *
 * The value \f$ \alpha \cdot [ u_0 - u(x) ] \f$ is a flux. It replaces the
 * integrand in the boundary integral for the variable \f$ u \f$.
 */
std::unique_ptr<TransportReflectionBoundaryCondition> createTransportReflectionBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& bc_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim, Process const& process,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace ProcessLib