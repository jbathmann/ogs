/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "GenericNonuniformNaturalBoundaryCondition.h"
#include "MeshLib/PropertyVector.h"
#include "NeumannTimeDependantBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
using NeumannTimeDependantBoundaryCondition =
    GenericNonuniformNaturalBoundaryCondition<
        NeumannTimeDependantBoundaryConditionData,
        NeumannTimeDependantBoundaryConditionLocalAssembler>;

std::unique_ptr<NeumannTimeDependantBoundaryCondition>
createNeumannTimeDependantBoundaryCondition(
    BaseLib::ConfigTree const& config,
    NumLib::LocalToGlobalIndexMap const& dof_table, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order, const MeshLib::Mesh& bulk_mesh);
}  // ProcessLib
