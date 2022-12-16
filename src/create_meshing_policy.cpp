// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include "stdafx.h"

#include <frost/create_meshing_policy.hpp>

#include <boost/make_shared.hpp>

#include <frost/frost_parameter_interface.hpp>

meshing_policy::ptr_type create_meshing_policy( const frost_parameter_interface& params ) {
    const int meshingMethod = params.get_meshing_method();

    switch( meshingMethod ) {
    case meshing_method::union_of_spheres:
        return boost::make_shared<union_of_spheres_meshing_policy>();
    case meshing_method::metaballs:
        return boost::make_shared<metaball_meshing_policy>( params.get_metaball_radius_scale(),
                                                            params.get_metaball_isosurface_level() );
    case meshing_method::zhu_bridson:
        return boost::make_shared<zhu_bridson_meshing_policy>( params.get_zhu_bridson_blend_radius_scale(),
                                                               params.get_zhu_bridson_enable_low_density_trimming(),
                                                               params.get_zhu_bridson_low_density_trimming_threshold(),
                                                               params.get_zhu_bridson_low_density_trimming_strength() );
    case meshing_method::anisotropic:
        return boost::make_shared<anisotropic_meshing_policy>(
            params.get_anisotropic_radius_scale(), params.get_anisotropic_isosurface_level(),
            params.get_anisotropic_max_anisotropy(), params.get_anisotropic_min_neighbor_count(),
            params.get_anisotropic_position_smoothing_weight() );
    default:
        throw std::runtime_error( "Internal Error: Unrecognized implicit surface meshing method: " +
                                  boost::lexical_cast<std::string>( meshingMethod ) );
    }
}
