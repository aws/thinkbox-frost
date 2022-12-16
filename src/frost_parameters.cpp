// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include "stdafx.h"

#include <frost/frost_parameters.hpp>

frost_parameters::frost_parameters()
    : m_meshingMethod( meshing_method::zhu_bridson )
    , m_metaballRadiusScale( 1.5f )
    , m_metaballIsosurfaceLevel( 0.3f )
    , m_zhuBridsonBlendRadiusScale( 2.2f )
    , m_zhuBridsonEnableLowDensityTrimming( true )
    , m_zhuBridsonLowDensityTrimmingThreshold( 1.f )
    , m_zhuBridsonLowDensityTrimmingStrength( 15.f )
    , m_anisotropicRadiusScale( 4.f )
    , m_anisotropicIsosurfaceLevel( 0.5f )
    , m_anisotropicMaxAnisotropy( 4.f )
    , m_anisotropicMinNeighborCount( 25 )
    , m_anisotropicPositionSmoothingWeight( 0.9f )
    , m_vertRefinementIterations( 3 )
    , m_meshingResolutionMode( meshing_resolution_mode::subdivide_max_radius )
    , m_meshingResolution( 1.5f )
    , m_meshingVoxelLength( 0.167f ) {}

int frost_parameters::get_meshing_method() const { return m_meshingMethod; }
void frost_parameters::set_meshing_method( int meshingMethod ) { m_meshingMethod = meshingMethod; }

float frost_parameters::get_metaball_radius_scale() const { return m_metaballRadiusScale; }
void frost_parameters::set_metaball_radius_scale( float metaballRadiusScale ) {
    m_metaballRadiusScale = metaballRadiusScale;
}
float frost_parameters::get_metaball_isosurface_level() const { return m_metaballIsosurfaceLevel; }
void frost_parameters::set_metaball_isosurface_level( float metaballIsosurfaceLevel ) {
    m_metaballIsosurfaceLevel = metaballIsosurfaceLevel;
}

float frost_parameters::get_zhu_bridson_blend_radius_scale() const { return m_zhuBridsonBlendRadiusScale; }
void frost_parameters::set_zhu_bridson_blend_radius_scale( float zhuBridsonBlendRadiusScale ) {
    m_zhuBridsonBlendRadiusScale = zhuBridsonBlendRadiusScale;
}
bool frost_parameters::get_zhu_bridson_enable_low_density_trimming() const {
    return m_zhuBridsonEnableLowDensityTrimming;
}
void frost_parameters::set_zhu_bridson_enable_low_density_trimming( bool zhuBridsonEnableLowDensityTrimming ) {
    m_zhuBridsonEnableLowDensityTrimming = zhuBridsonEnableLowDensityTrimming;
}
float frost_parameters::get_zhu_bridson_low_density_trimming_threshold() const {
    return m_zhuBridsonLowDensityTrimmingThreshold;
}
void frost_parameters::set_zhu_bridson_low_density_trimming_threshold( float zhuBridsonLowDensityTrimmingThreshold ) {
    m_zhuBridsonLowDensityTrimmingThreshold = zhuBridsonLowDensityTrimmingThreshold;
}
float frost_parameters::get_zhu_bridson_low_density_trimming_strength() const {
    return m_zhuBridsonLowDensityTrimmingStrength;
}
void frost_parameters::set_zhu_bridson_low_density_trimming_strength( float zhuBridsonLowDensityTrimmingStrength ) {
    m_zhuBridsonLowDensityTrimmingStrength = zhuBridsonLowDensityTrimmingStrength;
}

float frost_parameters::get_anisotropic_radius_scale() const { return m_anisotropicRadiusScale; }
void frost_parameters::set_anisotropic_radius_scale( float anisotropicRadiusScale ) {
    m_anisotropicRadiusScale = anisotropicRadiusScale;
}
float frost_parameters::get_anisotropic_isosurface_level() const { return m_anisotropicIsosurfaceLevel; }
void frost_parameters::set_anisotropic_isosurface_level( float anisotropicIsosurfaceLevel ) {
    m_anisotropicIsosurfaceLevel = anisotropicIsosurfaceLevel;
}
float frost_parameters::get_anisotropic_max_anisotropy() const { return m_anisotropicMaxAnisotropy; }
void frost_parameters::set_anisotropic_max_anisotropy( float anisotropicMaxAnisotropy ) {
    m_anisotropicMaxAnisotropy = anisotropicMaxAnisotropy;
}
int frost_parameters::get_anisotropic_min_neighbor_count() const { return m_anisotropicMinNeighborCount; }
void frost_parameters::set_anisotropic_min_neighbor_count( int anisotropicMinNeighborCount ) {
    m_anisotropicMinNeighborCount = anisotropicMinNeighborCount;
}
float frost_parameters::get_anisotropic_position_smoothing_weight() const {
    return m_anisotropicPositionSmoothingWeight;
}
void frost_parameters::set_anisotropic_position_smoothing_weight( float anisotropicPositionSmoothingWeight ) {
    m_anisotropicPositionSmoothingWeight = anisotropicPositionSmoothingWeight;
}

int frost_parameters::get_vert_refinement_iterations() const { return m_vertRefinementIterations; }
void frost_parameters::set_vert_refinement_iterations( int vertRefinementIterations ) {
    m_vertRefinementIterations = vertRefinementIterations;
}

int frost_parameters::get_meshing_resolution_mode() const { return m_meshingResolutionMode; }
void frost_parameters::set_meshing_resolution_mode( int meshingResolutionMode ) {
    m_meshingResolutionMode = meshingResolutionMode;
}

float frost_parameters::get_meshing_resolution() const { return m_meshingResolution; }
void frost_parameters::set_meshing_resolution( float meshingResolution ) { m_meshingResolution = meshingResolution; }

float frost_parameters::get_meshing_voxel_length() const { return m_meshingVoxelLength; }
void frost_parameters::set_meshing_voxel_length( float meshingVoxelLength ) {
    m_meshingVoxelLength = meshingVoxelLength;
}
