// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include "stdafx.h"

#include <frost/frost.hpp>

#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/make_shared.hpp>

#include <frantic/particles/streams/particle_array_particle_istream.hpp>
#include <frantic/volumetrics/implicitsurface/implicit_surface_to_trimesh3.hpp>

#include <frost/channel_names.hpp>
#include <frost/create_meshing_policy.hpp>
#include <frost/meshing_policy.hpp>
#include <frost/positive_radius_culled_particle_istream.hpp>

namespace frost {

namespace {

std::pair<float, float>
get_particle_radius_extrema( boost::shared_ptr<frantic::particles::streams::particle_istream> pin ) {
    if( pin->particle_count() == 0 ) {
        return std::make_pair( 0.f, 0.f );
    }
    frantic::channels::channel_map pcm = pin->get_channel_map();
    if( !pcm.has_channel( Frost_RADIUSCHANNEL ) ) {
        return std::make_pair( 0.f, 0.f );
    }
    frantic::channels::channel_const_cvt_accessor<float> radius =
        pcm.get_const_cvt_accessor<float>( Frost_RADIUSCHANNEL );
    std::vector<char> p( pcm.structure_size() );

    if( !pin->get_particle( p ) ) {
        return std::make_pair( 0.f, 0.f );
    }
    float minRadius = radius( p );
    float maxRadius = radius( p );
    while( pin->get_particle( p ) ) {
        const float r = radius( p );
        if( r > maxRadius ) {
            maxRadius = r;
        }
        if( r < minRadius ) {
            minRadius = r;
        }
    }

    return std::make_pair( minRadius, maxRadius );
}

std::pair<float, float> get_particle_radius_extrema( const frantic::particles::particle_array& pa ) {
    boost::shared_ptr<frantic::particles::streams::particle_istream> pin(
        new frantic::particles::streams::particle_array_particle_istream( pa ) );
    return get_particle_radius_extrema( pin );
}

void init_meshing_policy_operation_parameters( frantic::particles::particle_grid_tree& outParticles,
                                               float& outMaximumParticleRadius, float& outEffectRadiusScale,
                                               frantic::volumetrics::voxel_coord_system& outMeshingVCS,
                                               const frost_parameter_interface& params,
                                               boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                                               const frantic::channels::channel_propagation_policy& cpp,
                                               frantic::logging::progress_logger& progressLogger ) {
    boost::shared_ptr<meshing_policy> meshingPolicy( create_meshing_policy( params ) );

    outParticles.clear();

    pin.reset( new positive_radius_culled_particle_istream( pin, boost::none ) );
    pin->set_channel_map( pin->get_native_channel_map() );

    // Force position channel to be of type float since Frost is not designed to work with doubles yet
    frantic::channels::channel_map useChannelMap;
    const frantic::channels::channel_map& sourceChannelMap = pin->get_channel_map();
    for( std::size_t i = 0, ie = sourceChannelMap.channel_count(); i != ie; ++i ) {
        const frantic::channels::channel& channel = sourceChannelMap[i];
        if( channel.name() == Frost_POSITIONCHANNEL || channel.name() == Frost_RADIUSCHANNEL ) {
            useChannelMap.define_channel( channel.name(), channel.arity(), frantic::channels::data_type_float32 );
        } else if( cpp.is_channel_included( channel.name() ) ) {
            useChannelMap.define_channel( channel.name(), channel.arity(), channel.data_type() );
        }
    }
    useChannelMap.end_channel_definition();

    useChannelMap = frantic::volumetrics::implicitsurface::create_optimized_channel_map( useChannelMap );

    frantic::particles::particle_array pa( meshingPolicy->get_particle_buffer_channel_map( useChannelMap ) );
    frantic::logging::progress_logger_subinterval_tracker tracker( progressLogger, 0, 50 );
    pa.insert_particles( pin, progressLogger );
    progressLogger.update_progress( 100.f );

    tracker.reset( 50, 60 );
    meshingPolicy->postprocess_particle_buffer( pa, progressLogger );
    progressLogger.update_progress( 100.f );

    std::pair<float, float> particleRadiusExtrema = get_particle_radius_extrema( pa );
    const float minimumParticleRadius = particleRadiusExtrema.first;
    const float maximumParticleRadius = particleRadiusExtrema.second;

    // bool useRenderMeshing = get_boolean_attribute( inRender );

    float meshingResolution = 1.f;
    const int vertRefinement = params.get_vert_refinement_iterations();

    float maximumParticleRadiusForVoxelLength = maximumParticleRadius;
    // if( useRenderMeshing ){
    //	vertRefinement = get_int_attribute( inRenderVertRefinementIterations );
    // }else{
    //	vertRefinement = get_int_attribute( inViewportVertRefinementIterations );
    // }

    float meshingVoxelLength = 0;
    const int meshingResolutionMode = params.get_meshing_resolution_mode();
    if( meshingResolutionMode == meshing_resolution_mode::subdivide_max_radius ) {
        // const float meshingResolution = useRenderMeshing ? get_float_attribute( inRenderMeshingResolution ) :
        // get_float_attribute( inViewportMeshingResolution );
        const float meshingResolution = params.get_meshing_resolution();
        if( meshingResolution <= 0 ) {
            throw std::runtime_error( "Meshing resolution must be greater than zero" );
        }
        if( pa.particle_count() > 0 ) {
            meshingVoxelLength = maximumParticleRadius / meshingResolution;
        } else {
            // In case there are no particles, use a default value different
            // from zero to avoid errors below.
            meshingVoxelLength = 1;
        }
    } else if( meshingResolutionMode == meshing_resolution_mode::voxel_length ) {
        meshingVoxelLength = params.get_meshing_voxel_length();
        // if( useRenderMeshing ) {
        //	meshingVoxelLength = get_float_attribute( inRenderVoxelLength );
        // } else {
        //	meshingVoxelLength = get_float_attribute( inViewportVoxelLength );
        // }
        if( meshingVoxelLength <= 0 ) {
            throw std::runtime_error( "Voxel Length must be greater than zero." );
        }
    } else {
        throw std::runtime_error( "Unrecognized inMeshingResolutionMode: " +
                                  boost::lexical_cast<std::string>( meshingResolutionMode ) );
    }
    if( meshingVoxelLength <= 0 ) {
        throw std::runtime_error( "Voxel Length must be greater than zero." );
    }

    frantic::volumetrics::voxel_coord_system meshingVCS( frantic::graphics::vector3f( 0 ), meshingVoxelLength );
    if( meshingVCS.voxel_length() == 0 ) {
        return;
    }

    const float effectRadiusScale = meshingPolicy->get_effect_radius_scale( meshingVCS, minimumParticleRadius );

    // build the pgt
    const float particleVoxelLength =
        std::max<float>( maximumParticleRadius * effectRadiusScale, 0.5f * meshingVoxelLength );
    frantic::volumetrics::voxel_coord_system particleVCS( frantic::graphics::vector3f( 0 ), particleVoxelLength );
    if( particleVCS.voxel_length() == 0 ) {
        return;
    }

    tracker.reset( 60, 100 );
    outParticles.reset( pa.get_channel_map(), particleVCS );
    { // scope for radiusAcc and positionAcc
        frantic::channels::channel_accessor<float> radiusAcc =
            pa.get_channel_map().get_accessor<float>( Frost_RADIUSCHANNEL );
        frantic::channels::channel_accessor<frantic::graphics::vector3f> positionAcc =
            pa.get_channel_map().get_accessor<frantic::graphics::vector3f>( Frost_POSITIONCHANNEL );
        std::size_t infinitePositionCount = 0;
        const std::size_t particleCount = pa.size();
        for( size_t i = 0; i < particleCount; ++i ) {
            const char* p = pa[i];
            if( positionAcc( p ).is_finite() ) {
                if( radiusAcc( p ) > 0 ) {
                    outParticles.insert( p );
                }
            } else {
                ++infinitePositionCount;
            }

            if( i % 65536 == 0 ) {
                progressLogger.update_progress( i, particleCount );
            }
        }
        if( infinitePositionCount > 0 ) {
            FF_LOG( warning ) << "Removed " << boost::lexical_cast<frantic::tstring>( infinitePositionCount )
                              << " particle" << ( infinitePositionCount > 1 ? "s" : "" ) << " with non-finite Position"
                              << std::endl;
        }
    }
    progressLogger.update_progress( 100.f );

    outMaximumParticleRadius = maximumParticleRadius;
    outEffectRadiusScale = effectRadiusScale;
    outMeshingVCS = meshingVCS;
}

class invoke_evolve_mesh {
  public:
    invoke_evolve_mesh( int iterations, float relativeSpacing, float relaxWeight )
        : m_iterations( iterations )
        , m_relativeSpacing( relativeSpacing )
        , m_relaxWeight( relaxWeight ) {}

    void operator()( boost::shared_ptr<meshing_policy> policy, frantic::geometry::trimesh3& mesh,
                     frantic::particles::particle_grid_tree& particles,
                     const frantic::channels::channel_propagation_policy& cpp, float maxParticleRadius,
                     float effectRadiusScale, const frantic::volumetrics::voxel_coord_system& meshingVCS,
                     int vertRefinement, frantic::logging::progress_logger& progressLogger ) {
        policy->evolve_mesh( mesh, particles, cpp, maxParticleRadius, effectRadiusScale, meshingVCS, vertRefinement,
                             m_iterations, m_relativeSpacing * meshingVCS.voxel_length(), m_relaxWeight,
                             progressLogger );
    }

  private:
    int m_iterations;
    float m_relativeSpacing;
    float m_relaxWeight;
};

class delegated_sampler : public frost::sampler {
  public:
    delegated_sampler( boost::shared_ptr<frantic::particles::particle_grid_tree> particles,
                       frost::sampler::ptr_type sampler )
        : m_particles( particles )
        , m_sampler( sampler ) {}

    const frantic::channels::channel_map& get_channel_map() const { return m_particles->get_channel_map(); }

    void get_particles_in_range( const frantic::graphics::vector3f& position, std::vector<char*>& outParticles ) {
        return m_sampler->get_particles_in_range( position, outParticles );
    }

    virtual void get_sample_weights( const frantic::graphics::vector3f& position, const std::vector<char*>& particles,
                                     std::vector<float>& outWeights ) {
        return m_sampler->get_sample_weights( position, particles, outWeights );
    }

  private:
    boost::shared_ptr<frantic::particles::particle_grid_tree> m_particles;
    frost::sampler::ptr_type m_sampler;
};

} // anonymous namespace

typedef boost::function<void(
    boost::shared_ptr<meshing_policy>, frantic::geometry::trimesh3&, frantic::particles::particle_grid_tree&,
    const frantic::channels::channel_propagation_policy&, float, float, const frantic::volumetrics::voxel_coord_system&,
    int, frantic::logging::progress_logger& )>
    meshing_policy_operation_t;

template <class MeshingOperation>
void invoke_meshing_policy_operation( MeshingOperation operation, frantic::geometry::trimesh3& outMesh,
                                      std::size_t& outParticleCount, frost_parameter_interface& params,
                                      boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                                      const frantic::channels::channel_propagation_policy& cpp,
                                      frantic::logging::progress_logger& progressLogger ) {
    frantic::logging::progress_logger_subinterval_tracker progressInterval( progressLogger, 0, 100 );

    progressInterval.reset( 0, 10 );

    frantic::particles::particle_grid_tree particles;
    float maximumParticleRadius;
    float effectRadiusScale;
    frantic::volumetrics::voxel_coord_system meshingVCS;

    init_meshing_policy_operation_parameters( particles, maximumParticleRadius, effectRadiusScale, meshingVCS, params,
                                              pin, cpp, progressLogger );

    progressInterval.reset( 10, 100 );

    if( particles.particle_count() > 0 ) {
        boost::shared_ptr<meshing_policy> meshingPolicy( create_meshing_policy( params ) );
        const int vertRefinement = params.get_vert_refinement_iterations();

        operation( meshingPolicy, outMesh, particles, cpp, maximumParticleRadius, effectRadiusScale, meshingVCS,
                   vertRefinement, progressLogger );
    } else {
        // Propagate channels to empty mesh.
        // I'm doing this because inconsistent channels between mesh files
        // in a sequence was causing problems for Sequoia's mesh loader.
        const frantic::channels::channel_map& channelMap = particles.get_channel_map();
        for( std::size_t i = 0, ie = channelMap.channel_count(); i < ie; ++i ) {
            const frantic::channels::channel& ch = channelMap[i];
            if( cpp.is_channel_included( ch.name() ) ) {
                outMesh.add_vertex_channel_raw( ch.name(), ch.arity(), ch.data_type() );
            }
        }
    }

    progressLogger.update_progress( 100 );
}

void build_trimesh3( frantic::geometry::trimesh3& outMesh, std::size_t& outParticleCount,
                     frost_parameter_interface& params,
                     boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                     frantic::logging::progress_logger& progressLogger ) {
    frantic::channels::channel_propagation_policy cpp( true );
    cpp.add_channel( _T( "Color" ) );
    cpp.add_channel( _T( "TextureCoord" ) );
    cpp.add_channel( _T( "Mapping2" ) );
    cpp.add_channel( _T( "Velocity" ) );

    outMesh.clear();

    meshing_policy_operation_t operation(
        boost::bind( &meshing_policy::build_mesh, _1, _2, _3, _4, _5, _6, _7, _8, _9 ) );

    invoke_meshing_policy_operation( operation, outMesh, outParticleCount, params, pin, cpp, progressLogger );
}

void populate_mesh_channels( frantic::geometry::trimesh3& outMesh, frost_parameter_interface& params,
                             boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                             const frantic::channels::channel_propagation_policy& cpp,
                             frantic::logging::progress_logger& progressLogger ) {
    meshing_policy_operation_t operation(
        boost::bind( &meshing_policy::populate_mesh_channels, _1, _2, _3, _4, _5, _6, _7, _8, _9 ) );

    std::size_t dummy;
    invoke_meshing_policy_operation( operation, outMesh, dummy, params, pin, cpp, progressLogger );
}

void evolve_mesh( frantic::geometry::trimesh3& outMesh, frost_parameter_interface& params, int iterations,
                  float relativeSpacing, float relaxWeight,
                  boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                  const frantic::channels::channel_propagation_policy& cpp,
                  frantic::logging::progress_logger& progressLogger ) {
    invoke_evolve_mesh f( iterations, relativeSpacing, relaxWeight );

    std::size_t dummy;
    invoke_meshing_policy_operation( f, outMesh, dummy, params, pin, cpp, progressLogger );
}

frost::sampler::ptr_type create_sampler( frost_parameter_interface& params,
                                         boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                                         frantic::logging::progress_logger& progressLogger ) {
    boost::shared_ptr<frantic::particles::particle_grid_tree> particles =
        boost::make_shared<frantic::particles::particle_grid_tree>();

    float maximumParticleRadius;
    float effectRadiusScale;
    frantic::volumetrics::voxel_coord_system meshingVCS;
    frantic::channels::channel_propagation_policy cpp;

    init_meshing_policy_operation_parameters( *particles, maximumParticleRadius, effectRadiusScale, meshingVCS, params,
                                              pin, cpp, progressLogger );

    boost::shared_ptr<meshing_policy> meshingPolicy = create_meshing_policy( params );
    const int vertRefinement = params.get_vert_refinement_iterations();

    frost::sampler::ptr_type sampler = meshingPolicy->create_sampler( *particles, maximumParticleRadius,
                                                                      effectRadiusScale, meshingVCS, vertRefinement );
    return boost::make_shared<delegated_sampler>( particles, sampler );
}

} // namespace frost
