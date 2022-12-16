// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <boost/shared_ptr.hpp>

#include <frantic/channels/channel_map.hpp>
#include <frantic/channels/channel_propagation_policy.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/logging/progress_logger.hpp>
#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/particle_grid_tree.hpp>

#include "sampler.hpp"

class meshing_policy {
  public:
    typedef boost::shared_ptr<meshing_policy> ptr_type;

    virtual frantic::channels::channel_map
    get_particle_buffer_channel_map( const frantic::channels::channel_map& channelMap );
    virtual void postprocess_particle_buffer( frantic::particles::particle_array& particles,
                                              frantic::logging::progress_logger& progressLogger );

    virtual float get_effect_radius_scale( const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                           float minimumParticleRadius ) = 0;
    virtual void build_mesh( frantic::geometry::trimesh3& outMesh, frantic::particles::particle_grid_tree& particles,
                             const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                             float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                             int vertRefinement, frantic::logging::progress_logger& progressLogger ) = 0;

    virtual void populate_mesh_channels( frantic::geometry::trimesh3& outMesh,
                                         frantic::particles::particle_grid_tree& particles,
                                         const frantic::channels::channel_propagation_policy& cpp,
                                         float maxParticleRadius, float effectRadiusScale,
                                         const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertRefinement,
                                         frantic::logging::progress_logger& progressLogger ) = 0;

    virtual void evolve_mesh( frantic::geometry::trimesh3& outMesh, frantic::particles::particle_grid_tree& particles,
                              const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                              float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                              int vertRefinement, int iterations, float spacing, float relaxWeight,
                              frantic::logging::progress_logger& progressLogger ) = 0;

    virtual frost::sampler::ptr_type create_sampler( frantic::particles::particle_grid_tree& particles,
                                                     float maxParticleRadius, float effectRadiusScale,
                                                     const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                                     int vertRefinement ) const = 0;
};

class union_of_spheres_meshing_policy : public meshing_policy {
  public:
    float get_effect_radius_scale( const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                   float minimumParticleRadius );
    void build_mesh( frantic::geometry::trimesh3& outMesh, frantic::particles::particle_grid_tree& particles,
                     const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                     float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                     int vertRefinement, frantic::logging::progress_logger& progressLogger );

    void populate_mesh_channels( frantic::geometry::trimesh3& outMesh,
                                 frantic::particles::particle_grid_tree& particles,
                                 const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                                 float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                 int vertRefinement, frantic::logging::progress_logger& progressLogger );

    void evolve_mesh( frantic::geometry::trimesh3& outMesh, frantic::particles::particle_grid_tree& particles,
                      const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                      float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                      int vertRefinement, int iterations, float spacing, float relaxWeight,
                      frantic::logging::progress_logger& progressLogger );

    frost::sampler::ptr_type create_sampler( frantic::particles::particle_grid_tree& particles, float maxParticleRadius,
                                             float effectRadiusScale,
                                             const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                             int vertRefinement ) const;
};

class metaball_meshing_policy : public meshing_policy {
  public:
    metaball_meshing_policy( float radiusScale, float implicitThreshold );

    float get_effect_radius_scale( const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                   float minimumParticleRadius );
    void build_mesh( frantic::geometry::trimesh3& outMesh, frantic::particles::particle_grid_tree& particles,
                     const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                     float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                     int vertRefinement, frantic::logging::progress_logger& progressLogger );

    void populate_mesh_channels( frantic::geometry::trimesh3& outMesh,
                                 frantic::particles::particle_grid_tree& particles,
                                 const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                                 float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                 int vertRefinement, frantic::logging::progress_logger& progressLogger );

    void evolve_mesh( frantic::geometry::trimesh3& outMesh, frantic::particles::particle_grid_tree& particles,
                      const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                      float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                      int vertRefinement, int iterations, float spacing, float relaxWeight,
                      frantic::logging::progress_logger& progressLogger );

    frost::sampler::ptr_type create_sampler( frantic::particles::particle_grid_tree& particles, float maxParticleRadius,
                                             float effectRadiusScale,
                                             const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                             int vertRefinement ) const;

  private:
    float m_radiusScale;
    float m_implicitThreshold;
};

class zhu_bridson_meshing_policy : public meshing_policy {
  public:
    zhu_bridson_meshing_policy( float blendRadiusScale, bool enableLowDensityTrimming,
                                float lowDensityTrimmingThreshold, float lowDensityTrimmingStrength );

    float get_effect_radius_scale( const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                   float minimumParticleRadius );
    void build_mesh( frantic::geometry::trimesh3& outMesh, frantic::particles::particle_grid_tree& particles,
                     const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                     float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                     int vertRefinement, frantic::logging::progress_logger& progressLogger );

    void populate_mesh_channels( frantic::geometry::trimesh3& outMesh,
                                 frantic::particles::particle_grid_tree& particles,
                                 const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                                 float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                 int vertRefinement, frantic::logging::progress_logger& progressLogger );

    void evolve_mesh( frantic::geometry::trimesh3& outMesh, frantic::particles::particle_grid_tree& particles,
                      const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                      float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                      int vertRefinement, int iterations, float spacing, float relaxWeight,
                      frantic::logging::progress_logger& progressLogger );

    frost::sampler::ptr_type create_sampler( frantic::particles::particle_grid_tree& particles, float maxParticleRadius,
                                             float effectRadiusScale,
                                             const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                             int vertRefinement ) const;

  private:
    float m_blendRadiusScale;
    bool m_enableLowDensityTrimming;
    float m_lowDensityTrimmingThreshold;
    float m_lowDensityTrimmingStrength;
};

class anisotropic_meshing_policy : public meshing_policy {
  public:
    anisotropic_meshing_policy( float radiusScale, float isosurfaceLevel, float maxAnisotropy,
                                std::size_t minNeighborCount, float positionSmoothingWeight );

    frantic::channels::channel_map get_particle_buffer_channel_map( const frantic::channels::channel_map& channelMap );
    void postprocess_particle_buffer( frantic::particles::particle_array& particles,
                                      frantic::logging::progress_logger& progressLogger );

    float get_effect_radius_scale( const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                   float minimumParticleRadius );
    void build_mesh( frantic::geometry::trimesh3& outMesh, frantic::particles::particle_grid_tree& particles,
                     const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                     float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                     int vertRefinement, frantic::logging::progress_logger& progressLogger );

    void populate_mesh_channels( frantic::geometry::trimesh3& outMesh,
                                 frantic::particles::particle_grid_tree& particles,
                                 const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                                 float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                 int vertRefinement, frantic::logging::progress_logger& progressLogger );

    void evolve_mesh( frantic::geometry::trimesh3& outMesh, frantic::particles::particle_grid_tree& particles,
                      const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                      float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                      int vertRefinement, int iterations, float spacing, float relaxWeight,
                      frantic::logging::progress_logger& progressLogger );

    frost::sampler::ptr_type create_sampler( frantic::particles::particle_grid_tree& particles, float maxParticleRadius,
                                             float effectRadiusScale,
                                             const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                             int vertRefinement ) const;

  private:
    float m_radiusScale;
    float m_isosurfaceLevel;
    float m_windowScale;
    float m_maxAnisotropy;
    std::size_t m_minNeighborCount;
    bool m_enablePositionSmoothing;
    float m_positionSmoothingWindowScale;
    float m_positionSmoothingWeight;
    frantic::tstring m_volumeChannelName;
};
