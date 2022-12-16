// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include "stdafx.h"

#include <frost/meshing_policy.hpp>

#include <boost/make_shared.hpp>

#include <frantic/geometry/relaxation.hpp>

#include <frantic/volumetrics/implicitsurface/calculate_particle_anisotropic_params.hpp>
#include <frantic/volumetrics/implicitsurface/implicit_surface_to_trimesh3.hpp>
#include <frantic/volumetrics/implicitsurface/particle_implicit_surface_policies.hpp>
#include <frantic/volumetrics/run_tree.hpp>

#ifdef min
#undef min
#endif

#define SPARSE_MESHING_VOXEL_COUNT_THRESHOLD 500000

using namespace frantic::volumetrics;

namespace {

double get_volume_double( const frantic::graphics::boundbox3& bounds ) {
    if( bounds.is_empty() ) {
        return 0;
    } else {
        return double( bounds.xsize() ) * double( bounds.ysize() ) * double( bounds.zsize() );
    }
}

bool use_sparse_meshing( const frantic::particles::particle_grid_tree& particles, float maximumEffectDistance,
                         const frantic::volumetrics::voxel_coord_system& meshingVCS ) {
    // frantic::graphics::boundbox3f worldParticleBounds = particles.compute_particle_bounds();
    // worldParticleBounds.expand( maximumEffectDistance );
    // const frantic::graphics::boundbox3 voxelParticleBounds = meshingVCS.get_voxel_bounds( worldParticleBounds );
    // return get_volume_double( voxelParticleBounds ) > SPARSE_MESHING_VOXEL_COUNT_THRESHOLD;
    return true;
}

template <class Derived>
class sampler_impl : public frost::sampler {
  public:
    const frantic::channels::channel_map& get_channel_map() const { return derived().m_isp.get_channel_map(); }

    void get_particles_in_range( const frantic::graphics::vector3f& position, std::vector<char*>& outParticles ) {
        return derived().m_isp.get_particles_in_range( position, outParticles );
    }

    void get_sample_weights( const frantic::graphics::vector3f& position, const std::vector<char*>& particles,
                             std::vector<float>& outWeights ) {
        return derived().m_isp.get_sample_weights( position, particles, outWeights );
    }

    Derived& derived() { return *static_cast<Derived*>( this ); }
    const Derived& derived() const { return *static_cast<const Derived*>( this ); }
};

class union_of_spheres_sampler : public sampler_impl<union_of_spheres_sampler> {
  public:
    union_of_spheres_sampler( frantic::particles::particle_grid_tree& particles, float maximumParticleRadius,
                              float particleRadiusToEffectRadiusScale, float implicitThreshold,
                              const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertexRefinement )
        : m_isp( particles, maximumParticleRadius, particleRadiusToEffectRadiusScale, implicitThreshold, meshingVCS,
                 vertexRefinement ) {}

    frantic::volumetrics::implicitsurface::particle_union_of_spheres_is_policy m_isp;
};

class metaball_sampler : public sampler_impl<metaball_sampler> {
  public:
    metaball_sampler( frantic::particles::particle_grid_tree& particles, float maximumParticleRadius,
                      float particleRadiusToEffectRadiusScale, float implicitThreshold,
                      const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertexRefinement )
        : m_isp( particles, maximumParticleRadius, particleRadiusToEffectRadiusScale, implicitThreshold, meshingVCS,
                 vertexRefinement ) {}

    frantic::volumetrics::implicitsurface::particle_metaball_is_policy m_isp;
};

class zhu_bridson_sampler : public sampler_impl<zhu_bridson_sampler> {
  public:
    zhu_bridson_sampler( frantic::particles::particle_grid_tree& particles, float maxParticleRadius, float effectRadius,
                         float lowDensityTrimmingDensity, float lowDensityTrimmingStrength,
                         const voxel_coord_system& meshingVCS, int vertexRefinement )
        : m_isp( particles, maxParticleRadius, effectRadius, lowDensityTrimmingDensity, lowDensityTrimmingStrength,
                 meshingVCS, vertexRefinement ) {}

    frantic::volumetrics::implicitsurface::particle_zhu_bridson_is_policy m_isp;
};

class anisotropic_sampler : public sampler_impl<anisotropic_sampler> {
  public:
    anisotropic_sampler( frantic::particles::particle_grid_tree& particles, float implicitThreshold,
                         const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertexRefinement )
        : m_isp( particles, implicitThreshold, meshingVCS, vertexRefinement ) {}

    frantic::volumetrics::implicitsurface::particle_anisotropic_is_policy m_isp;
};

} // anonymous namespace

frantic::channels::channel_map
meshing_policy::get_particle_buffer_channel_map( const frantic::channels::channel_map& channelMap ) {
    return channelMap;
}

void meshing_policy::postprocess_particle_buffer( frantic::particles::particle_array& particles,
                                                  frantic::logging::progress_logger& progressLogger ) {}

float union_of_spheres_meshing_policy::get_effect_radius_scale(
    const frantic::volumetrics::voxel_coord_system& meshingVCS, float minimumParticleRadius ) {
    return ( minimumParticleRadius == 0 )
               ? 2.f
               : std::min( 2.f, 1.f + 2.f * meshingVCS.voxel_length() / minimumParticleRadius );
}

void union_of_spheres_meshing_policy::build_mesh( frantic::geometry::trimesh3& outMesh,
                                                  frantic::particles::particle_grid_tree& particles,
                                                  const frantic::channels::channel_propagation_policy& cpp,
                                                  float maximumParticleRadius, float effectRadiusScale,
                                                  const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                                  int vertRefinement,
                                                  frantic::logging::progress_logger& progressLogger ) {
    float implicitThreshold = 2.f * meshingVCS.voxel_length();
    bool useSparseMeshing = use_sparse_meshing( particles, effectRadiusScale * maximumParticleRadius, meshingVCS );
    if( useSparseMeshing ) {
        frantic::volumetrics::implicitsurface::union_of_spheres_convert_sparse_particles_to_trimesh3(
            particles, cpp, maximumParticleRadius, effectRadiusScale, implicitThreshold, meshingVCS, vertRefinement,
            outMesh, progressLogger );
    } else {
        frantic::volumetrics::implicitsurface::union_of_spheres_convert_particles_to_trimesh3(
            particles, cpp, maximumParticleRadius, effectRadiusScale, implicitThreshold, meshingVCS, vertRefinement,
            outMesh, progressLogger );
    }
}

void union_of_spheres_meshing_policy::populate_mesh_channels(
    frantic::geometry::trimesh3& outMesh, frantic::particles::particle_grid_tree& particles,
    const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius, float effectRadiusScale,
    const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertRefinement,
    frantic::logging::progress_logger& progressLogger ) {
    float implicitThreshold = 2.f * meshingVCS.voxel_length();
    frantic::volumetrics::implicitsurface::particle_union_of_spheres_is_policy policy(
        particles, maxParticleRadius, effectRadiusScale, implicitThreshold, meshingVCS, vertRefinement );
    policy.populate_mesh_channels( outMesh, cpp );
}

void union_of_spheres_meshing_policy::evolve_mesh( frantic::geometry::trimesh3& outMesh,
                                                   frantic::particles::particle_grid_tree& particles,
                                                   const frantic::channels::channel_propagation_policy& cpp,
                                                   float maxParticleRadius, float effectRadiusScale,
                                                   const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                                   int vertRefinement, int iterations, float spacing, float relaxWeight,
                                                   frantic::logging::progress_logger& progressLogger ) {
    float implicitThreshold = 2.f * meshingVCS.voxel_length();
    frantic::volumetrics::implicitsurface::particle_union_of_spheres_is_policy policy(
        particles, maxParticleRadius, effectRadiusScale, implicitThreshold, meshingVCS, vertRefinement );
    frantic::geometry::relaxation::evolve_mesh_to_implicit_surface( outMesh, policy, iterations, spacing, relaxWeight,
                                                                    &progressLogger );
}

frost::sampler::ptr_type union_of_spheres_meshing_policy::create_sampler(
    frantic::particles::particle_grid_tree& particles, float maxParticleRadius, float effectRadiusScale,
    const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertRefinement ) const {
    float implicitThreshold = 2.f * meshingVCS.voxel_length();
    return boost::make_shared<union_of_spheres_sampler>( boost::ref( particles ), maxParticleRadius, effectRadiusScale,
                                                         implicitThreshold, meshingVCS, vertRefinement );
}

metaball_meshing_policy::metaball_meshing_policy( float radiusScale, float implicitThreshold )
    : m_radiusScale( radiusScale )
    , m_implicitThreshold( implicitThreshold ) {}

float metaball_meshing_policy::get_effect_radius_scale( const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                                        float minimumParticleRadius ) {
    return m_radiusScale;
}

void metaball_meshing_policy::build_mesh( frantic::geometry::trimesh3& outMesh,
                                          frantic::particles::particle_grid_tree& particles,
                                          const frantic::channels::channel_propagation_policy& cpp,
                                          float maximumParticleRadius, float effectRadiusScale,
                                          const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                          int vertRefinement, frantic::logging::progress_logger& progressLogger ) {
    float implicitThreshold = m_implicitThreshold;
    bool useSparseMeshing = use_sparse_meshing( particles, effectRadiusScale * maximumParticleRadius, meshingVCS );
    if( useSparseMeshing ) {
        frantic::volumetrics::implicitsurface::metaball_convert_sparse_particles_to_trimesh3(
            particles, cpp, maximumParticleRadius, effectRadiusScale, implicitThreshold, meshingVCS, vertRefinement,
            outMesh, progressLogger );
    } else {
        frantic::volumetrics::implicitsurface::metaball_convert_particles_to_trimesh3(
            particles, cpp, maximumParticleRadius, effectRadiusScale, implicitThreshold, meshingVCS, vertRefinement,
            outMesh, progressLogger );
    }
}

void metaball_meshing_policy::populate_mesh_channels( frantic::geometry::trimesh3& outMesh,
                                                      frantic::particles::particle_grid_tree& particles,
                                                      const frantic::channels::channel_propagation_policy& cpp,
                                                      float maxParticleRadius, float effectRadiusScale,
                                                      const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                                      int vertRefinement,
                                                      frantic::logging::progress_logger& progressLogger ) {
    float implicitThreshold = m_implicitThreshold;
    frantic::volumetrics::implicitsurface::particle_metaball_is_policy policy(
        particles, maxParticleRadius, effectRadiusScale, implicitThreshold, meshingVCS, vertRefinement );
    policy.populate_mesh_channels( outMesh, cpp );
}

void metaball_meshing_policy::evolve_mesh( frantic::geometry::trimesh3& outMesh,
                                           frantic::particles::particle_grid_tree& particles,
                                           const frantic::channels::channel_propagation_policy& cpp,
                                           float maxParticleRadius, float effectRadiusScale,
                                           const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                           int vertRefinement, int iterations, float spacing, float relaxWeight,
                                           frantic::logging::progress_logger& progressLogger ) {
    float implicitThreshold = m_implicitThreshold;
    frantic::volumetrics::implicitsurface::particle_metaball_is_policy policy(
        particles, maxParticleRadius, effectRadiusScale, implicitThreshold, meshingVCS, vertRefinement );
    frantic::geometry::relaxation::evolve_mesh_to_implicit_surface( outMesh, policy, iterations, spacing, relaxWeight,
                                                                    &progressLogger );
}

frost::sampler::ptr_type metaball_meshing_policy::create_sampler(
    frantic::particles::particle_grid_tree& particles, float maxParticleRadius, float effectRadiusScale,
    const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertRefinement ) const {
    float implicitThreshold = m_implicitThreshold;
    return boost::make_shared<metaball_sampler>( boost::ref( particles ), maxParticleRadius, effectRadiusScale,
                                                 implicitThreshold, meshingVCS, vertRefinement );
}

zhu_bridson_meshing_policy::zhu_bridson_meshing_policy( float blendRadiusScale, bool enableLowDensityTrimming,
                                                        float lowDensityTrimmingThreshold,
                                                        float lowDensityTrimmingStrength )
    : m_blendRadiusScale( blendRadiusScale )
    , m_enableLowDensityTrimming( enableLowDensityTrimming )
    , m_lowDensityTrimmingThreshold( lowDensityTrimmingThreshold )
    , m_lowDensityTrimmingStrength( lowDensityTrimmingStrength ) {}

float zhu_bridson_meshing_policy::get_effect_radius_scale( const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                                           float minimumParticleRadius ) {
    return m_blendRadiusScale;
}

void zhu_bridson_meshing_policy::build_mesh( frantic::geometry::trimesh3& outMesh,
                                             frantic::particles::particle_grid_tree& particles,
                                             const frantic::channels::channel_propagation_policy& cpp,
                                             float maximumParticleRadius, float effectRadiusScale,
                                             const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                             int vertRefinement, frantic::logging::progress_logger& progressLogger ) {
    const bool enableLowDensityTrimming =
        m_enableLowDensityTrimming; // m_zhuBridsonEnableLowDensityTrimming.at_time( t );
    float lowDensityTrimmingDensity = enableLowDensityTrimming ? m_lowDensityTrimmingThreshold : 0.f;
    float lowDensityTrimmingStrength = m_lowDensityTrimmingStrength;
    bool useSparseMeshing = use_sparse_meshing( particles, effectRadiusScale * maximumParticleRadius, meshingVCS );
    if( useSparseMeshing ) {
        frantic::volumetrics::implicitsurface::zhu_bridson_convert_sparse_particles_to_trimesh3(
            particles, cpp, maximumParticleRadius, effectRadiusScale, lowDensityTrimmingDensity,
            lowDensityTrimmingStrength, meshingVCS, vertRefinement, outMesh, progressLogger );
    } else {
        frantic::volumetrics::implicitsurface::zhu_bridson_convert_particles_to_trimesh3(
            particles, cpp, maximumParticleRadius, effectRadiusScale, lowDensityTrimmingDensity,
            lowDensityTrimmingStrength, meshingVCS, vertRefinement, outMesh, progressLogger );
    }
}

void zhu_bridson_meshing_policy::populate_mesh_channels( frantic::geometry::trimesh3& outMesh,
                                                         frantic::particles::particle_grid_tree& particles,
                                                         const frantic::channels::channel_propagation_policy& cpp,
                                                         float maxParticleRadius, float effectRadiusScale,
                                                         const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                                         int vertRefinement,
                                                         frantic::logging::progress_logger& progressLogger ) {
    const bool enableLowDensityTrimming =
        m_enableLowDensityTrimming; // m_zhuBridsonEnableLowDensityTrimming.at_time( t );
    float lowDensityTrimmingDensity = enableLowDensityTrimming ? m_lowDensityTrimmingThreshold : 0.f;
    float lowDensityTrimmingStrength = m_lowDensityTrimmingStrength;
    frantic::volumetrics::implicitsurface::particle_zhu_bridson_is_policy policy(
        particles, maxParticleRadius, effectRadiusScale, lowDensityTrimmingDensity, lowDensityTrimmingStrength,
        meshingVCS, vertRefinement );
    policy.populate_mesh_channels( outMesh, cpp );
}

void zhu_bridson_meshing_policy::evolve_mesh( frantic::geometry::trimesh3& outMesh,
                                              frantic::particles::particle_grid_tree& particles,
                                              const frantic::channels::channel_propagation_policy& cpp,
                                              float maxParticleRadius, float effectRadiusScale,
                                              const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                              int vertRefinement, int iterations, float spacing, float relaxWeight,
                                              frantic::logging::progress_logger& progressLogger ) {
    const bool enableLowDensityTrimming =
        m_enableLowDensityTrimming; // m_zhuBridsonEnableLowDensityTrimming.at_time( t );
    float lowDensityTrimmingDensity = enableLowDensityTrimming ? m_lowDensityTrimmingThreshold : 0.f;
    float lowDensityTrimmingStrength = m_lowDensityTrimmingStrength;
    frantic::volumetrics::implicitsurface::particle_zhu_bridson_is_policy policy(
        particles, maxParticleRadius, effectRadiusScale, lowDensityTrimmingDensity, lowDensityTrimmingStrength,
        meshingVCS, vertRefinement );
    frantic::geometry::relaxation::evolve_mesh_to_implicit_surface( outMesh, policy, iterations, spacing, relaxWeight,
                                                                    &progressLogger );
}

frost::sampler::ptr_type zhu_bridson_meshing_policy::create_sampler(
    frantic::particles::particle_grid_tree& particles, float maxParticleRadius, float effectRadiusScale,
    const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertRefinement ) const {
    const bool enableLowDensityTrimming =
        m_enableLowDensityTrimming; // m_zhuBridsonEnableLowDensityTrimming.at_time( t );
    float lowDensityTrimmingDensity = enableLowDensityTrimming ? m_lowDensityTrimmingThreshold : 0.f;
    float lowDensityTrimmingStrength = m_lowDensityTrimmingStrength;
    return boost::make_shared<zhu_bridson_sampler>( boost::ref( particles ), maxParticleRadius, effectRadiusScale,
                                                    lowDensityTrimmingDensity, lowDensityTrimmingStrength, meshingVCS,
                                                    vertRefinement );
}

anisotropic_meshing_policy::anisotropic_meshing_policy( float radiusScale, float isosurfaceLevel, float maxAnisotropy,
                                                        std::size_t minNeighborCount, float positionSmoothingWeight )
    : m_radiusScale( radiusScale )
    , m_isosurfaceLevel( isosurfaceLevel )
    , m_windowScale( 2.f )
    , m_maxAnisotropy( maxAnisotropy )
    , m_minNeighborCount( minNeighborCount )
    , m_enablePositionSmoothing( true )
    , m_positionSmoothingWindowScale( 2.f )
    , m_positionSmoothingWeight( positionSmoothingWeight )
    , m_volumeChannelName( _T("__Volume") ) {}

frantic::channels::channel_map anisotropic_meshing_policy::get_particle_buffer_channel_map(
    const frantic::channels::channel_map& originalChannelMap ) {
    frantic::channels::channel_map result( originalChannelMap );

    result.append_channel<float>( m_volumeChannelName );
    result.append_channel<float>( _T("__MaxDistance") );
    result.append_channel( _T("__Anisotropy"), 6, frantic::channels::data_type_float32 );
    result.append_channel<float>( _T("__invcsv") );

    return result;
}

void anisotropic_meshing_policy::postprocess_particle_buffer( frantic::particles::particle_array& particles,
                                                              frantic::logging::progress_logger& progressLogger ) {
    using namespace frantic::volumetrics::implicitsurface;

    const float compactSupportScale = m_radiusScale;
    const float anisotropyWindowScale = m_windowScale;

    const float kr = m_maxAnisotropy;
    const std::size_t ne = m_minNeighborCount;
    const float implicitThreshold = m_isosurfaceLevel;

    // Calculate anisotropy matrices
    frantic::logging::progress_logger_subinterval_tracker tracker( progressLogger, 0, 33 );
    calculate_anisotropy( particles, compactSupportScale, compactSupportScale * anisotropyWindowScale, kr,
                          /*kn, ks,*/ ne, progressLogger );
    progressLogger.update_progress( 100.f );

    // Smooth particle positions
    const float smoothingEffectRadiusScale = m_positionSmoothingWindowScale;
    const float lambda = m_positionSmoothingWeight;
    const bool enableSmoothing = m_enablePositionSmoothing && lambda > 0.01f;

    if( enableSmoothing ) {
        tracker.reset( 33, 67 );
        smooth_particle_positions( particles, smoothingEffectRadiusScale * compactSupportScale, lambda,
                                   progressLogger );
        progressLogger.update_progress( 100.f );
    }

    // Calculate particle density
    tracker.reset( 67, 100 );
    calculate_volume_with_anisotropic_kernel( particles, m_volumeChannelName, progressLogger );
    progressLogger.update_progress( 100.f );
}

float anisotropic_meshing_policy::get_effect_radius_scale( const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                                           float minimumParticleRadius ) {
    return m_radiusScale;
}

void anisotropic_meshing_policy::build_mesh( frantic::geometry::trimesh3& outMesh,
                                             frantic::particles::particle_grid_tree& particles,
                                             const frantic::channels::channel_propagation_policy& cpp,
                                             float maximumParticleRadius, float effectRadiusScale,
                                             const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                             int vertRefinement, frantic::logging::progress_logger& progressLogger ) {
    using namespace frantic::volumetrics::implicitsurface;
    const float implicitThreshold = m_isosurfaceLevel;
    bool useSparseMeshing = use_sparse_meshing( particles, effectRadiusScale * maximumParticleRadius, meshingVCS );
    if( useSparseMeshing ) {
        anisotropic_convert_sparse_particles_to_trimesh3( particles, cpp, implicitThreshold, meshingVCS, vertRefinement,
                                                          outMesh, progressLogger );
    } else {
        anisotropic_convert_particles_to_trimesh3( particles, cpp, implicitThreshold, meshingVCS, vertRefinement,
                                                   outMesh, progressLogger );
    }
}

void anisotropic_meshing_policy::populate_mesh_channels( frantic::geometry::trimesh3& outMesh,
                                                         frantic::particles::particle_grid_tree& particles,
                                                         const frantic::channels::channel_propagation_policy& cpp,
                                                         float maxParticleRadius, float effectRadiusScale,
                                                         const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                                         int vertRefinement,
                                                         frantic::logging::progress_logger& progressLogger ) {
    using namespace frantic::volumetrics::implicitsurface;
    const float implicitThreshold = m_isosurfaceLevel;
    frantic::volumetrics::implicitsurface::particle_anisotropic_is_policy policy( particles, implicitThreshold,
                                                                                  meshingVCS, vertRefinement );
    policy.populate_mesh_channels( outMesh, cpp );
}

void anisotropic_meshing_policy::evolve_mesh( frantic::geometry::trimesh3& outMesh,
                                              frantic::particles::particle_grid_tree& particles,
                                              const frantic::channels::channel_propagation_policy& cpp,
                                              float maxParticleRadius, float /*effectRadiusScale*/,
                                              const frantic::volumetrics::voxel_coord_system& meshingVCS,
                                              int vertRefinement, int iterations, float spacing, float relaxWeight,
                                              frantic::logging::progress_logger& progressLogger ) {
    using namespace frantic::volumetrics::implicitsurface;
    const float implicitThreshold = m_isosurfaceLevel;
    frantic::volumetrics::implicitsurface::particle_anisotropic_is_policy policy( particles, implicitThreshold,
                                                                                  meshingVCS, vertRefinement );
    frantic::geometry::relaxation::evolve_mesh_to_implicit_surface( outMesh, policy, iterations, spacing, relaxWeight,
                                                                    &progressLogger );
}

frost::sampler::ptr_type anisotropic_meshing_policy::create_sampler(
    frantic::particles::particle_grid_tree& particles, float maxParticleRadius, float effectRadiusScale,
    const frantic::volumetrics::voxel_coord_system& meshingVCS, int vertRefinement ) const {
    const float implicitThreshold = m_isosurfaceLevel;
    return boost::make_shared<anisotropic_sampler>( boost::ref( particles ), implicitThreshold, meshingVCS,
                                                    vertRefinement );
}
