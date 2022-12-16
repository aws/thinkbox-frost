// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include "geometry_meshing_parameters_interface.hpp"

#include <frantic/channels/channel_map.hpp>
#include <frantic/geometry/trimesh3.hpp>
#include <frantic/math/hash.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

frantic::particles::particle_istream_ptr prepare_radius_channel( frantic::particles::particle_istream_ptr pin,
                                                                 float radius, bool overwrite );

frantic::particles::particle_istream_ptr prepare_color_channel( frantic::particles::particle_istream_ptr pin );

void conform_mesh_channel_types( frantic::geometry::trimesh3& mesh,
                                 const frantic::channels::channel_map& vertexChannelMap, std::size_t skipChannelIndex,
                                 const frantic::channels::channel_map& faceChannelMap );

void get_channel_accessors(
    frantic::geometry::trimesh3& mesh, const frantic::channels::channel_map& vertexChannelMap,
    const frantic::channels::channel_map& faceChannelMap,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outVertexAccessors,
    std::vector<frantic::geometry::trimesh3_face_channel_general_accessor>& outFaceAccessors );

void build_plane_mesh( frantic::geometry::trimesh3& outMesh, const frantic::channels::channel_map& vertexChannelMap,
                       const frantic::channels::channel_map& faceChannelMap,
                       frost::hard_edge_type::option hardEdgeType );

void build_sprite_mesh( frantic::geometry::trimesh3& outMesh, const frantic::channels::channel_map& vertexChannelMap,
                        const frantic::channels::channel_map& faceChannelMap,
                        frost::hard_edge_type::option hardEdgeType );

void build_tetrahedron_mesh( frantic::geometry::trimesh3& outMesh,
                             const frantic::channels::channel_map& vertexChannelMap,
                             const frantic::channels::channel_map& faceChannelMap,
                             frost::hard_edge_type::option hardEdgeType );

// void build_pyramid_mesh( frantic::geometry::trimesh3 & outMesh, const std::vector<std::string> & channelNames, const
// std::vector<frantic::channels::channel_general_accessor> & inputChannels );

void build_box_mesh( frantic::geometry::trimesh3& outMesh, const frantic::channels::channel_map& vertexChannelMap,
                     const frantic::channels::channel_map& faceChannelMap, frost::hard_edge_type::option hardEdgeType );

void build_sphere_mesh( frantic::geometry::trimesh3& outMesh, const frantic::channels::channel_map& vertexChannelMap,
                        const frantic::channels::channel_map& faceChannelMap,
                        frost::hard_edge_type::option hardEdgeType );

boost::shared_ptr<frantic::particles::streams::particle_istream>
randomize_radius_by_id( const boost::shared_ptr<frantic::particles::streams::particle_istream> inParticles,
                        const frantic::tstring& idChannel, const frantic::tstring& radiusChannel, const float variation,
                        const int seed );

inline float get_unit_random_from_id( boost::int32_t id, boost::uint32_t randomSeed ) {
    boost::uint32_t hash = hashword( reinterpret_cast<boost::uint32_t*>( &id ), 1, randomSeed );
    return static_cast<float>( hash ) / static_cast<float>( std::numeric_limits<boost::uint32_t>::max() );
}

inline boost::int32_t id_randomized_geometry( boost::int32_t id, boost::int32_t elements, boost::int32_t randomSeed ) {
    boost::uint32_t hashedID = hashword( reinterpret_cast<boost::uint32_t*>( &id ), 1, randomSeed );
    float unitRandom =
        static_cast<float>( hashedID ) / static_cast<float>( std::numeric_limits<boost::uint32_t>::max() );
    return int( elements * ( 1.f - unitRandom ) );
}

inline float id_randomized_rotation( boost::int32_t id, float maxRotation, boost::int32_t randomSeed = 12345 ) {
    boost::uint32_t hashedID = hashword( reinterpret_cast<boost::uint32_t*>( &id ), 1, randomSeed );
    float unitRandom =
        static_cast<float>( hashedID ) / static_cast<float>( std::numeric_limits<boost::uint32_t>::max() );
    return maxRotation * ( 1.f - unitRandom );
}

inline float id_randomized_radius( boost::int32_t id, float radius, float variation, boost::uint32_t randomSeed ) {
    boost::uint32_t hashedID = hashword( reinterpret_cast<boost::uint32_t*>( &id ), 1, randomSeed );
    float unitRandom =
        static_cast<float>( hashedID ) / static_cast<float>( std::numeric_limits<boost::uint32_t>::max() );
    // Only shrink the particles, because we are using the original radius as the basis for acceleration grids
    return radius * ( 1.f - unitRandom * variation );
}

template <class T>
T mod_trunc_neginf( T a, T b ) {
    int result = a % b;
    if( result < 0 ) {
        result += b;
    }
    return result;
}
