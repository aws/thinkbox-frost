// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/cstdint.hpp>

#include <frantic/graphics/transform4f.hpp>
#include <frantic/graphics/units.hpp>
#include <frantic/graphics/vector3f.hpp>

namespace frost {

namespace geometry_type {
enum option {
    plane,
    sprite,
    tetrahedron,
    custom_geometry,
    box,
    sphere20,
    //
    count
};
};

// How to encode hard edges in the built-in shapes?
namespace hard_edge_type {
enum option {
    // Use the "Normal" vertex channel
    vertex_normal,
    // Use the 3ds Max "SmoothingGroup" face channel
    smoothing_group,
    //
    count
};
}

namespace geometry_selection_mode {
enum option {
    cycle,
    random_by_id,
    shapeindex_channel,
    //
    count
};
};

namespace geometry_sample_time_base_mode {
enum option {
    time_0,
    current_time,
    //
    count
};
};

namespace geometry_sample_time_offset_mode {
enum option {
    no_offset,
    random_by_id,
    abstime_channel,
    timeoffset_channel,
    geomtime_channel,
    //
    count
};
};

namespace geometry_orientation_mode {
enum option {
    look_at,
    match_object,
    orientation_channel,
    vector_channel,
    specify,
    //
    count
};
};

namespace geometry_orientation_divergence_axis_space {
enum option {
    world,
    local,
    //
    count
};
};

class geometry_meshing_parameters_interface {
  public:
    virtual double get_time() = 0;

    virtual geometry_type::option get_geometry_type() = 0;
    virtual hard_edge_type::option get_hard_edge_type() = 0;
    virtual frantic::graphics::coordinate_system::option get_coordinate_system() {
        return frantic::graphics::coordinate_system::right_handed_zup;
    }

    virtual geometry_selection_mode::option get_geometry_selection_mode() = 0;
    virtual boost::uint32_t get_geometry_selection_seed() = 0;

    virtual geometry_sample_time_base_mode::option get_geometry_sample_time_base_mode() = 0;
    virtual geometry_sample_time_offset_mode::option get_geometry_sample_time_offset_mode() = 0;
    virtual double get_geometry_sample_time_max_random_offset() = 0;
    virtual boost::uint32_t get_geometry_sample_time_seed() = 0;

    virtual double frames_to_seconds( double frames ) = 0;

    virtual geometry_orientation_mode::option get_geometry_orientation_mode() = 0;
    virtual frantic::graphics::vector3f get_geometry_orientation() = 0;
    virtual frantic::graphics::vector3f get_geometry_look_at_position() = 0;
    virtual frantic::graphics::vector3f get_geometry_look_at_orientation() = 0;
    virtual frantic::tstring get_geometry_orientation_vector_channel() = 0;
    virtual float get_geometry_orientation_divergence() = 0;

    virtual bool get_orientation_restrict_divergence_axis() = 0;
    virtual frantic::graphics::vector3f get_orientation_divergence_axis() = 0;
    virtual geometry_orientation_divergence_axis_space::option get_geometry_orientation_divergence_axis_space() = 0;
};

} // namespace frost
