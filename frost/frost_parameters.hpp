// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frost/frost_parameter_interface.hpp>

class frost_parameters : public frost_parameter_interface {
  public:
    frost_parameters();

    int get_meshing_method() const;
    void set_meshing_method( int meshingMethod );

    float get_metaball_radius_scale() const;
    void set_metaball_radius_scale( float metaballRadiusScale );
    float get_metaball_isosurface_level() const;
    void set_metaball_isosurface_level( float metaballIsosurfaceLevel );

    float get_zhu_bridson_blend_radius_scale() const;
    void set_zhu_bridson_blend_radius_scale( float zhuBridsonBlendRadiusScale );
    bool get_zhu_bridson_enable_low_density_trimming() const;
    void set_zhu_bridson_enable_low_density_trimming( bool zhuBridsonEnableLowDensityTrimming );
    float get_zhu_bridson_low_density_trimming_threshold() const;
    void set_zhu_bridson_low_density_trimming_threshold( float zhuBridsonLowDensityTrimmingThreshold );
    float get_zhu_bridson_low_density_trimming_strength() const;
    void set_zhu_bridson_low_density_trimming_strength( float zhuBridsonLowDensityTrimmingStrength );

    float get_anisotropic_radius_scale() const;
    void set_anisotropic_radius_scale( float anisotropicRadiusScale );
    float get_anisotropic_isosurface_level() const;
    void set_anisotropic_isosurface_level( float anisotropicIsosurfaceLevel );
    float get_anisotropic_max_anisotropy() const;
    void set_anisotropic_max_anisotropy( float anisotropicMaxAnisotropy );
    int get_anisotropic_min_neighbor_count() const;
    void set_anisotropic_min_neighbor_count( int anisotropicMinNeighborCount );
    float get_anisotropic_position_smoothing_weight() const;
    void set_anisotropic_position_smoothing_weight( float anisotropicPositionSmoothingWeight );

    int get_vert_refinement_iterations() const;
    void set_vert_refinement_iterations( int vertRefinementIterations );

    int get_meshing_resolution_mode() const;
    void set_meshing_resolution_mode( int meshingResolutionMode );

    float get_meshing_resolution() const;
    void set_meshing_resolution( float meshingResolution );

    float get_meshing_voxel_length() const;
    void set_meshing_voxel_length( float meshingVoxelLength );

  private:
    int m_meshingMethod;

    float m_metaballRadiusScale;
    float m_metaballIsosurfaceLevel;

    float m_zhuBridsonBlendRadiusScale;
    bool m_zhuBridsonEnableLowDensityTrimming;
    float m_zhuBridsonLowDensityTrimmingThreshold;
    float m_zhuBridsonLowDensityTrimmingStrength;

    float m_anisotropicRadiusScale;
    float m_anisotropicIsosurfaceLevel;
    float m_anisotropicMaxAnisotropy;
    int m_anisotropicMinNeighborCount;
    float m_anisotropicPositionSmoothingWeight;

    int m_vertRefinementIterations;

    int m_meshingResolutionMode;

    float m_meshingResolution;
    float m_meshingVoxelLength;
};
