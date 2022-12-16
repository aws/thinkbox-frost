// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include "stdafx.h"

#include <frost/channel_names.hpp>
#include <frost/utility.hpp>

#include <boost/array.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>

#include <frantic/graphics/transform4f.hpp>
#include <frantic/graphics/vector3f.hpp>

#include <frantic/particles/streams/apply_function_particle_istream.hpp>
#include <frantic/particles/streams/set_channel_particle_istream.hpp>

using namespace frantic::graphics;
using namespace frantic::geometry;

frantic::particles::particle_istream_ptr prepare_radius_channel( frantic::particles::particle_istream_ptr pin,
                                                                 float radius, bool overwrite ) {
    if( !pin ) {
        throw std::runtime_error( "prepare_radius_channel Error: particle stream is NULL" );
    }

    frantic::channels::channel_map pcm = pin->get_channel_map();
    if( overwrite || !pcm.has_channel( Frost_RADIUSCHANNEL ) ) {
        pin.reset(
            new frantic::particles::streams::set_channel_particle_istream<float>( pin, Frost_RADIUSCHANNEL, radius ) );
        pin->set_channel_map( pin->get_native_channel_map() );
        pcm = pin->get_channel_map();
    }

    return pin;
}

namespace {

frantic::graphics::vector3f create_vector3f_from_float( float value ) { return frantic::graphics::vector3f( value ); }

} // anonymous namespace

frantic::particles::particle_istream_ptr prepare_color_channel( frantic::particles::particle_istream_ptr pin ) {
    if( !pin ) {
        throw std::runtime_error( "prepare_color_channel Error: particle stream is NULL" );
    }

    frantic::channels::channel_map pcm = pin->get_channel_map();
    if( !pcm.has_channel( Frost_COLORCHANNEL ) && pcm.has_channel( Frost_INTENSITYCHANNEL ) ) {
        static const boost::array<frantic::tstring, 1> intensityChannel = { Frost_INTENSITYCHANNEL };
        pin.reset(
            new frantic::particles::streams::apply_function_particle_istream<frantic::graphics::vector3f( float )>(
                pin, create_vector3f_from_float, Frost_COLORCHANNEL, intensityChannel ) );
        pin->set_channel_map( pin->get_native_channel_map() );
    }

    return pin;
}

namespace {

void add_mesh_channels( trimesh3& mesh, const frantic::channels::channel_map& vertexChannelMap,
                        const frantic::channels::channel_map& faceChannelMap, bool textureChannelHasCustomFaces,
                        int& outTextureChannel ) {
    using std::size_t;

    for( size_t i = 0; i < vertexChannelMap.channel_count(); ++i ) {
        const frantic::channels::channel& ch = vertexChannelMap[i];
        bool hasCustomFaces = false;

        if( ch.name() == Frost_TEXTURECOORDCHANNEL ) {
            outTextureChannel = static_cast<int>( i );
            hasCustomFaces = textureChannelHasCustomFaces;
        }
        if( hasCustomFaces ) {
            mesh.add_vertex_channel_raw( ch.name(), ch.arity(), ch.data_type(), mesh.vertex_count(), hasCustomFaces );
        } else {
            mesh.add_vertex_channel_raw( ch.name(), ch.arity(), ch.data_type() );
        }
    }

    for( size_t i = 0; i < faceChannelMap.channel_count(); ++i ) {
        const frantic::channels::channel& ch = faceChannelMap[i];

        mesh.add_face_channel_raw( ch.name(), ch.arity(), ch.data_type() );
    }
}

void add_faceted_vertex_normal_channel( frantic::geometry::trimesh3& mesh ) {
    mesh.add_vertex_channel<vector3f>( _T("Normal"), mesh.face_count(), true );

    trimesh3_vertex_channel_accessor<vector3f> normalAcc = mesh.get_vertex_channel_accessor<vector3f>( _T("Normal") );

    for( std::size_t i = 0, ie = mesh.face_count(); i < ie; ++i ) {
        const vector3& face = mesh.get_face( i );
        const vector3f a = mesh.get_vertex( face[0] );
        const vector3f b = mesh.get_vertex( face[1] );
        const vector3f c = mesh.get_vertex( face[2] );
        normalAcc[i] = triangle_normal( a, b, c );

        normalAcc.face( i ).set( static_cast<boost::int32_t>( i ) );
    }
}

} // anonymous namespace

void build_plane_mesh( frantic::geometry::trimesh3& outMesh, const frantic::channels::channel_map& vertexChannelMap,
                       const frantic::channels::channel_map& faceChannelMap,
                       frost::hard_edge_type::option hardEdgeType ) {
    frantic::geometry::trimesh3 primitiveMesh;

    const float f = 1.f; // 0.5f;

    std::vector<vector3f> verts;
    verts.push_back( vector3f( -f, -f, 0 ) );
    verts.push_back( vector3f( f, -f, 0 ) );
    verts.push_back( vector3f( -f, f, 0 ) );
    verts.push_back( vector3f( f, f, 0 ) );

    for( unsigned i = 0; i < verts.size(); ++i )
        primitiveMesh.add_vertex( verts[i] );
    primitiveMesh.add_face( 0, 1, 2 );
    primitiveMesh.add_face( 1, 3, 2 );

    // add the channels, noting if there exists a texture channel
    int textureChannel = -1;
    add_mesh_channels( primitiveMesh, vertexChannelMap, faceChannelMap, false, textureChannel );

    // if there wasnt a texture channel, add it
    if( textureChannel < 0 ) {
        textureChannel = (int)vertexChannelMap.channel_count(); // channelNames.size();
        primitiveMesh.add_vertex_channel<vector3f>( Frost_TEXTURECOORDCHANNEL );
    }

    // populate the texture channel with a standard set of coords
    trimesh3_vertex_channel_cvt_accessor<vector3f> texAcc =
        primitiveMesh.get_vertex_channel_cvt_accessor<vector3f>( Frost_TEXTURECOORDCHANNEL );
    texAcc.set( 0, vector3f( 0, 0, 0 ) );
    texAcc.set( 1, vector3f( 1, 0, 0 ) );
    texAcc.set( 2, vector3f( 0, 1, 0 ) );
    texAcc.set( 3, vector3f( 1, 1, 0 ) );

    // TODO: is this really necessary?
    if( hardEdgeType == frost::hard_edge_type::vertex_normal ) {
        add_faceted_vertex_normal_channel( primitiveMesh );
    }

    outMesh.swap( primitiveMesh );
}

void build_sprite_mesh( frantic::geometry::trimesh3& outMesh, const frantic::channels::channel_map& vertexChannelMap,
                        const frantic::channels::channel_map& faceChannelMap,
                        frost::hard_edge_type::option hardEdgeType ) {
    trimesh3 primitiveMesh;

    trimesh3 planeMesh;
    build_plane_mesh( planeMesh, vertexChannelMap, faceChannelMap, hardEdgeType );

    // Create vector of channels that should be transformed. Velocity should not be transformed.
    std::vector<std::pair<frantic::tstring, frantic::geometry::trimesh3::vector_type>> transformChannels;
    transformChannels.push_back( std::pair<frantic::tstring, frantic::geometry::trimesh3::vector_type>(
        _T("Normal"), frantic::geometry::trimesh3::NORMAL ) );

    primitiveMesh.combine( transform4f::identity(), planeMesh, transformChannels );
    primitiveMesh.combine( transform4f::from_cubeface( cube_face::CF_X_POS ), planeMesh, transformChannels );
    primitiveMesh.combine( transform4f::from_cubeface( cube_face::CF_Y_POS ), planeMesh, transformChannels );

    primitiveMesh.swap( outMesh );
}

void build_tetrahedron_mesh( frantic::geometry::trimesh3& outMesh,
                             const frantic::channels::channel_map& vertexChannelMap,
                             const frantic::channels::channel_map& faceChannelMap,
                             frost::hard_edge_type::option hardEdgeType ) {
    trimesh3 primitiveMesh;

    const float radius = 1.f;

    const float f = radius / sqrtf( 3 );

    primitiveMesh.add_vertex( f, f, f );
    primitiveMesh.add_vertex( -f, -f, f );
    primitiveMesh.add_vertex( -f, f, -f );
    primitiveMesh.add_vertex( f, -f, -f );

    primitiveMesh.add_face( 0, 2, 1 );
    primitiveMesh.add_face( 0, 1, 3 );
    primitiveMesh.add_face( 0, 3, 2 );
    primitiveMesh.add_face( 1, 2, 3 );

    // add the channels, noting if there exists a texture channel
    int textureChannel = -1;
    add_mesh_channels( primitiveMesh, vertexChannelMap, faceChannelMap, true, textureChannel );

    // if there wasnt a texture channel, add it
    if( textureChannel < 0 ) {
        textureChannel = (int)vertexChannelMap.channel_count(); // channelNames.size();
        primitiveMesh.add_vertex_channel<vector3f>( Frost_TEXTURECOORDCHANNEL, 4, true );
    }

    // populate the texture channel with a standard set of coords
    trimesh3_vertex_channel_cvt_accessor<vector3f> texAcc =
        primitiveMesh.get_vertex_channel_cvt_accessor<vector3f>( Frost_TEXTURECOORDCHANNEL );
    texAcc.set_vertex_count( 4 );
    texAcc.set( 0, vector3f( 0, 0, 0 ) );
    texAcc.set( 1, vector3f( 1.f, 0, 0 ) );
    texAcc.set( 2, vector3f( 0.5f, 0.866f, 0 ) );
    texAcc.set( 3, vector3f( 0, 0, 0 ) );

    texAcc.face( 0 ) = vector3( 0, 1, 2 );
    texAcc.face( 1 ) = vector3( 0, 1, 2 );
    texAcc.face( 2 ) = vector3( 0, 1, 2 );
    texAcc.face( 3 ) = vector3( 0, 1, 2 );

    switch( hardEdgeType ) {
    case frost::hard_edge_type::vertex_normal:
        add_faceted_vertex_normal_channel( primitiveMesh );
        break;
    case frost::hard_edge_type::smoothing_group: {
        primitiveMesh.add_face_channel<int>( Frost_SMOOTHINGGROUPCHANNEL );
        trimesh3_face_channel_accessor<int> sgAcc =
            primitiveMesh.get_face_channel_accessor<int>( Frost_SMOOTHINGGROUPCHANNEL );
        sgAcc[0] = 1;
        sgAcc[1] = 2;
        sgAcc[2] = 4;
        sgAcc[3] = 8;
    } break;
    }

    outMesh.swap( primitiveMesh );
}

void build_box_mesh( frantic::geometry::trimesh3& outMesh, const frantic::channels::channel_map& vertexChannelMap,
                     const frantic::channels::channel_map& faceChannelMap,
                     frost::hard_edge_type::option hardEdgeType ) {
    typedef boost::uint16_t material_id_t;
    trimesh3 primitiveMesh;

    const float r = 1.f;

    primitiveMesh.set_to_box( frantic::graphics::boundbox3f( -r, r, -r, r, -r, r ) );

    // populate the texture channel with a standard set of coords
    const boost::array<vector3f, 12> textureCoords = { vector3f( 0, 0, 0 ), vector3f( 1, 0, 0 ), vector3f( 0, 1, 0 ),
                                                       vector3f( 1, 1, 0 ), vector3f( 0, 0, 0 ), vector3f( 1, 0, 0 ),
                                                       vector3f( 0, 1, 0 ), vector3f( 1, 1, 0 ), vector3f( 0, 0, 0 ),
                                                       vector3f( 1, 0, 0 ), vector3f( 0, 1, 0 ), vector3f( 1, 1, 0 ) };
    const boost::array<vector3, 12> textureCoordFaces = {
        vector3( 8, 9, 11 ), vector3( 11, 10, 8 ), vector3( 10, 8, 9 ), vector3( 9, 11, 10 ),
        vector3( 7, 6, 4 ),  vector3( 4, 5, 7 ),   vector3( 5, 7, 6 ),  vector3( 6, 4, 5 ),
        vector3( 1, 3, 2 ),  vector3( 2, 0, 1 ),   vector3( 3, 2, 0 ),  vector3( 0, 1, 3 ) };

    // add the channels, noting if there exists a texture channel
    int textureChannel = -1;
    add_mesh_channels( primitiveMesh, vertexChannelMap, faceChannelMap, true, textureChannel );

    // if there wasnt a texture channel, add it
    if( textureChannel < 0 ) {
        textureChannel = (int)vertexChannelMap.channel_count();
        primitiveMesh.add_vertex_channel<vector3f>( Frost_TEXTURECOORDCHANNEL, textureCoords.size(), true );
    }

    trimesh3_vertex_channel_cvt_accessor<frantic::graphics::vector3f> texAcc =
        primitiveMesh.get_vertex_channel_cvt_accessor<frantic::graphics::vector3f>( Frost_TEXTURECOORDCHANNEL );
    if( textureCoordFaces.size() != texAcc.face_count() ) {
        throw std::runtime_error( "build_box_mesh: Mismatch between face data array and face count" );
    }
    for( std::size_t i = 0; i < textureCoordFaces.size(); ++i ) {
        texAcc.face( i ) = textureCoordFaces[i];
    }
    texAcc.set_vertex_count( textureCoords.size() );
    for( std::size_t i = 0; i < textureCoords.size(); ++i ) {
        texAcc.set( i, textureCoords[i] );
    }

    switch( hardEdgeType ) {
    case frost::hard_edge_type::vertex_normal:
        add_faceted_vertex_normal_channel( primitiveMesh );
        break;
    case frost::hard_edge_type::smoothing_group: {
        primitiveMesh.add_face_channel<int>( Frost_SMOOTHINGGROUPCHANNEL );
        trimesh3_face_channel_accessor<int> sgAcc =
            primitiveMesh.get_face_channel_accessor<int>( Frost_SMOOTHINGGROUPCHANNEL );
        sgAcc[0] = 1;
        sgAcc[1] = 1;
        sgAcc[2] = 2;
        sgAcc[3] = 2;
        sgAcc[4] = 4;
        sgAcc[5] = 4;
        sgAcc[6] = 8;
        sgAcc[7] = 8;
        sgAcc[8] = 16;
        sgAcc[9] = 16;
        sgAcc[10] = 32;
        sgAcc[11] = 32;
    } break;
    }

    if( !primitiveMesh.has_face_channel( Frost_MESH_MTLINDEX_CHANNEL_NAME ) ) {
        primitiveMesh.add_face_channel<material_id_t>( Frost_MESH_MTLINDEX_CHANNEL_NAME );

        const boost::array<boost::uint8_t, 12> faceMaterials = { 1, 1, 0, 0, 4, 4, 5, 5, 2, 2, 3, 3 };

        trimesh3_face_channel_accessor<material_id_t> mtlAcc =
            primitiveMesh.get_face_channel_accessor<material_id_t>( Frost_MESH_MTLINDEX_CHANNEL_NAME );
        for( std::size_t i = 0; i < faceMaterials.size(); ++i ) {
            mtlAcc[i] = faceMaterials[i];
        }
    }

    outMesh.swap( primitiveMesh );
}

void build_sphere_mesh( frantic::geometry::trimesh3& outMesh, const frantic::channels::channel_map& vertexChannelMap,
                        const frantic::channels::channel_map& faceChannelMap,
                        frost::hard_edge_type::option hardEdgeType ) {
    trimesh3 primitiveMesh;

    const float radius = 1.f;
    primitiveMesh.set_to_icosahedron( vector3f( 0 ), radius );

    // add the channels, noting if there exists a texture channel
    int textureChannel = -1;
    add_mesh_channels( primitiveMesh, vertexChannelMap, faceChannelMap, false, textureChannel );

    // if there wasnt a texture channel, add it
    if( textureChannel < 0 ) {
        textureChannel = (int)vertexChannelMap.channel_count(); // channelNames.size();
        primitiveMesh.add_vertex_channel<vector3f>( Frost_TEXTURECOORDCHANNEL );
    }

    // populate the texture channel with a standard set of coords
    trimesh3_vertex_channel_cvt_accessor<vector3f> texAcc =
        primitiveMesh.get_vertex_channel_cvt_accessor<vector3f>( Frost_TEXTURECOORDCHANNEL );

    for( size_t vertexIndex = 0; vertexIndex < primitiveMesh.vertex_count(); ++vertexIndex ) {
        vector3f v;
        for( int i = 0; i < 3; ++i ) {
            v[i] = 0.5f * ( primitiveMesh.get_vertex( vertexIndex )[i] + 1.f ) / radius;
        }
        texAcc.set( vertexIndex, v );
    }

    if( hardEdgeType == frost::hard_edge_type::vertex_normal ) {
        primitiveMesh.build_vertex_normals();
    }

    outMesh.swap( primitiveMesh );
}

boost::shared_ptr<frantic::particles::streams::particle_istream>
randomize_radius_by_id( const boost::shared_ptr<frantic::particles::streams::particle_istream> inParticles,
                        const frantic::tstring& idChannel, const frantic::tstring& radiusChannel, const float variation,
                        const int seed ) {
    boost::array<frantic::tstring, 2> inputParamNames = { idChannel, radiusChannel };
    float randomVariation = 0.01f * variation;
    inParticles->set_channel_map( inParticles->get_native_channel_map() );
    frantic::channels::channel_map channels = inParticles->get_channel_map();
    if( channels.has_channel( idChannel ) ) {
        frantic::channels::data_type_t idType;
        std::size_t idArity;
        channels.get_channel_definition( idChannel, idType, idArity );

        if( idArity != 1 ) {
            throw std::runtime_error( "Error: randomize_radius_by_id: ID channel must have arity of 1." );
        }

        return boost::make_shared<
            frantic::particles::streams::apply_function_particle_istream<float( const boost::int32_t&, const float& )>>(
            inParticles, boost::bind( &id_randomized_radius, _1, _2, randomVariation, seed ), radiusChannel,
            inputParamNames );
    }

    throw std::runtime_error( "Error: randomize_radius_by_id: input particle stream has no Radius channel" );
}

void conform_mesh_channel_types( frantic::geometry::trimesh3& mesh,
                                 const frantic::channels::channel_map& vertexChannelMap, std::size_t skipChannelIndex,
                                 const frantic::channels::channel_map& faceChannelMap ) {
    for( std::size_t i = 0, ie = faceChannelMap.channel_count(); i < ie; ++i ) {
        const frantic::channels::channel& ch = faceChannelMap[i];
        const frantic::tstring& channelName = ch.name();
        const std::size_t arity = ch.arity();
        const frantic::channels::data_type_t dataType = ch.data_type();

        if( mesh.has_face_channel( channelName ) ) {
            frantic::geometry::trimesh3_face_channel_general_accessor acc =
                mesh.get_face_channel_general_accessor( channelName );
            if( acc.arity() != arity || acc.data_type() != dataType ) {
                mesh.erase_face_channel( channelName );
                mesh.add_face_channel_raw( channelName, arity, dataType );
            }
        } else {
            mesh.add_face_channel_raw( channelName, arity, dataType );
        }
    }

    for( std::size_t i = 0, ie = vertexChannelMap.channel_count(); i < ie; ++i ) {
        const frantic::channels::channel& ch = vertexChannelMap[i];
        const frantic::tstring& channelName = ch.name();
        const std::size_t arity = ch.arity();
        const frantic::channels::data_type_t dataType = ch.data_type();

        if( i == skipChannelIndex ) {
            if( mesh.has_vertex_channel( channelName ) ) {
                // do nothing
            } else {
                mesh.add_vertex_channel_raw( channelName, arity, dataType );
                trimesh3_vertex_channel_general_accessor acc = mesh.get_vertex_channel_general_accessor( channelName );
                const std::size_t primitiveSize = acc.primitive_size();
                for( std::size_t vertexIndex = 0, vertexIndexEnd = acc.size(); vertexIndex < vertexIndexEnd;
                     ++vertexIndex ) {
                    memset( acc.data( vertexIndex ), 0, primitiveSize );
                }
            }
        } else {
            if( mesh.has_vertex_channel( channelName ) ) {
                trimesh3_vertex_channel_general_accessor acc = mesh.get_vertex_channel_general_accessor( channelName );

                if( acc.arity() != arity || acc.data_type() != dataType ) {
                    const std::size_t dataCount = acc.size();
                    const bool hasCustomFaces = acc.has_custom_faces();
                    std::vector<vector3> customFaceBuffer;

                    if( hasCustomFaces ) {
                        const std::size_t faceCount = acc.face_count();
                        if( faceCount != mesh.face_count() ) {
                            throw std::runtime_error(
                                "conform_vertex_channel_types Error: mismatch between number of faces in mesh (" +
                                boost::lexical_cast<std::string>( acc.arity() ) +
                                ") and number of faces in vertex channel \'" +
                                frantic::strings::to_string( channelName ) + "\' (" +
                                boost::lexical_cast<std::string>( arity ) + ")." );
                        }
                        customFaceBuffer.resize( faceCount );
                        for( std::size_t faceIndex = 0; faceIndex < faceCount; ++faceIndex ) {
                            customFaceBuffer[faceIndex] = acc.face( faceIndex );
                        }
                    }

                    mesh.erase_vertex_channel( channelName );

                    mesh.add_vertex_channel_raw( channelName, arity, dataType, dataCount, hasCustomFaces );

                    if( hasCustomFaces ) {
                        acc = mesh.get_vertex_channel_general_accessor( channelName );
                        for( std::size_t faceIndex = 0, faceIndexEnd = mesh.face_count(); faceIndex < faceIndexEnd;
                             ++faceIndex ) {
                            acc.face( faceIndex ) = customFaceBuffer[faceIndex];
                        }
                    }
                }
            } else {
                mesh.add_vertex_channel_raw( channelName, arity, dataType );
            }
        }
    }
}

void get_channel_accessors(
    frantic::geometry::trimesh3& mesh, const frantic::channels::channel_map& vertexChannelMap,
    const frantic::channels::channel_map& faceChannelMap,
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor>& outVertexAccessors,
    std::vector<frantic::geometry::trimesh3_face_channel_general_accessor>& outFaceAccessors ) {
    outVertexAccessors.clear();
    outFaceAccessors.clear();

    for( std::size_t i = 0; i < vertexChannelMap.channel_count(); ++i ) {
        const frantic::channels::channel& ch = vertexChannelMap[i];
        const frantic::tstring& channelName = ch.name();
        outVertexAccessors.push_back( mesh.get_vertex_channel_general_accessor( channelName ) );
        trimesh3_vertex_channel_general_accessor acc = outVertexAccessors.back();
        if( acc.data_type() != ch.data_type() ) {
            throw std::runtime_error(
                "get_channel_accessors Error: data type mismatch between mesh (" +
                frantic::strings::to_string(
                    frantic::channels::channel_data_type_str( acc.arity(), acc.data_type() ) ) +
                ") and particles (" +
                frantic::strings::to_string( frantic::channels::channel_data_type_str( ch.arity(), ch.data_type() ) ) +
                ") in vertex channel \'" + frantic::strings::to_string( channelName ) + "\'." );
        }
        if( acc.arity() != ch.arity() ) {
            throw std::runtime_error( "get_channel_accessors Error: arity mismatch between mesh (" +
                                      boost::lexical_cast<std::string>( acc.arity() ) + ") and particles (" +
                                      boost::lexical_cast<std::string>( ch.arity() ) + ") in vertex channel \'" +
                                      frantic::strings::to_string( channelName ) + "\'." );
        }
    }

    for( std::size_t i = 0; i < faceChannelMap.channel_count(); ++i ) {
        const frantic::channels::channel& ch = faceChannelMap[i];
        const frantic::tstring& channelName = ch.name();
        outFaceAccessors.push_back( mesh.get_face_channel_general_accessor( channelName ) );
        trimesh3_face_channel_general_accessor acc = outFaceAccessors.back();
        if( acc.data_type() != ch.data_type() ) {
            throw std::runtime_error(
                "get_channel_accessors Error: data type mismatch between mesh (" +
                frantic::strings::to_string(
                    frantic::channels::channel_data_type_str( acc.arity(), acc.data_type() ) ) +
                ") and particles (" +
                frantic::strings::to_string( frantic::channels::channel_data_type_str( ch.arity(), ch.data_type() ) ) +
                ") in face channel \'" + frantic::strings::to_string( channelName ) + "\'." );
        }
        if( acc.arity() != ch.arity() ) {
            throw std::runtime_error( "get_channel_accessors Error: arity mismatch between mesh (" +
                                      boost::lexical_cast<std::string>( acc.arity() ) + ") and particles (" +
                                      boost::lexical_cast<std::string>( ch.arity() ) + ") in face channel \'" +
                                      frantic::strings::to_string( channelName ) + "\'." );
        }
    }
}
