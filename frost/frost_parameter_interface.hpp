// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

namespace meshing_method {
enum option {
    union_of_spheres,
    metaballs,
    zhu_bridson,
    anisotropic,
    //
    count
};
};

namespace meshing_resolution_mode {
enum option {
    subdivide_max_radius,
    voxel_length,
    //
    count
};
};

class frost_parameter_interface {
  public:
    virtual int get_meshing_method() const = 0;

    virtual float get_metaball_radius_scale() const = 0;
    virtual float get_metaball_isosurface_level() const = 0;

    virtual float get_zhu_bridson_blend_radius_scale() const = 0;
    virtual bool get_zhu_bridson_enable_low_density_trimming() const = 0;
    virtual float get_zhu_bridson_low_density_trimming_threshold() const = 0;
    virtual float get_zhu_bridson_low_density_trimming_strength() const = 0;

    virtual float get_anisotropic_radius_scale() const = 0;
    virtual float get_anisotropic_isosurface_level() const = 0;
    virtual float get_anisotropic_max_anisotropy() const = 0;
    virtual int get_anisotropic_min_neighbor_count() const = 0;
    virtual float get_anisotropic_position_smoothing_weight() const = 0;

    virtual int get_vert_refinement_iterations() const = 0;

    virtual int get_meshing_resolution_mode() const = 0;

    virtual float get_meshing_resolution() const = 0;

    virtual float get_meshing_voxel_length() const = 0;
};
