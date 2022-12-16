// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "stdafx.h"

#include <frost/frost.hpp>
#include <frost/frost_parameters.hpp>
#include <frost/geometry_meshing.hpp>
#include <frost/geometry_meshing_parameters_interface.hpp>

using namespace frost;
// The following class is used to mock a geometry_meshing_parameters object, which are defined
// in the FrostMY project as maya_geometry_meshing_parameters, and in FrostMax as
// max3d_geometry_meshing_parameters
class geometry_meshing_parameters : public geometry_meshing_parameters_interface {
    double get_time() { return 0; }
    geometry_type::option get_geometry_type() { return m_geometryType; }
    hard_edge_type::option get_hard_edge_type() { return hard_edge_type::vertex_normal; }
    geometry_selection_mode::option get_geometry_selection_mode() { return geometry_selection_mode::cycle; }
    boost::uint32_t get_geometry_selection_seed() { return 0; }
    geometry_sample_time_base_mode::option get_geometry_sample_time_base_mode() {
        return geometry_sample_time_base_mode::time_0;
    }
    geometry_sample_time_offset_mode::option get_geometry_sample_time_offset_mode() {
        return geometry_sample_time_offset_mode::no_offset;
    }
    double get_geometry_sample_time_max_random_offset() { return 0; }
    boost::uint32_t get_geometry_sample_time_seed() { return 0; }
    frantic::graphics::coordinate_system::option get_coordinate_system() { return m_coordinateSystem; }
    double frames_to_seconds( double frames ) { return frames * 30; }

    geometry_orientation_mode::option get_geometry_orientation_mode() { return m_geometryOrientationMode; }
    frantic::graphics::vector3f get_geometry_orientation() { return frantic::graphics::vector3f( 0, 0, 0 ); }
    frantic::graphics::vector3f get_geometry_look_at_position() { return frantic::graphics::vector3f( 0, 0, 0 ); }
    frantic::graphics::vector3f get_geometry_look_at_orientation() { return frantic::graphics::vector3f( 0, 0, 0 ); }
    frantic::tstring get_geometry_orientation_vector_channel() { return frantic::tstring( _T("Position") ); }
    float get_geometry_orientation_divergence() { return 0; }

    bool get_orientation_restrict_divergence_axis() { return false; }
    frantic::graphics::vector3f get_orientation_divergence_axis() { return frantic::graphics::vector3f( 0, 0, 0 ); }
    geometry_orientation_divergence_axis_space::option get_geometry_orientation_divergence_axis_space() {
        return geometry_orientation_divergence_axis_space::world;
    }

    frantic::graphics::coordinate_system::option m_coordinateSystem;
    geometry_type::option m_geometryType;
    geometry_orientation_mode::option m_geometryOrientationMode;

  public:
    void set_coordinate_system( frantic::graphics::coordinate_system::option coordSystem ) {
        m_coordinateSystem = coordSystem;
    }

    void set_geometry_type( const geometry_type::option type ) { m_geometryType = type; }

    void set_geometry_orientation_mode( const geometry_orientation_mode::option mode ) {
        m_geometryOrientationMode = mode;
    }
};