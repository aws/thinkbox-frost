// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include <stdafx.h>

#include <boost/array.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/scoped_ptr.hpp>

#include <ImathVec.h>

#include <frantic/channels/channel_propagation_policy.hpp>
#include <frantic/graphics/units.hpp>

#include <frost/channel_names.hpp>
#include <frost/geometry_meshing.hpp>
#include <frost/utility.hpp>

using namespace frantic::geometry;
using namespace frantic::graphics;

namespace {

// TODO: I used Imath's normalize() here because, at the time I wrote this
// (December 2013), our own normalize function had strange behaviour for small
// values due to a hard-coded epsilon.  I think we could adapt the code to use
// our own normalize function now.
inline frantic::graphics::vector3f normalized( const frantic::graphics::vector3f& v ) {
    Imath::Vec3<float> imfVec( v.x, v.y, v.z );
    imfVec.normalize();
    return frantic::graphics::vector3f( imfVec.x, imfVec.y, imfVec.z );
}

} // anonymous namespace

get_shape_number_strategy::~get_shape_number_strategy() {}

boost::int32_t use_first_shape_strategy::get_shape_number( boost::int64_t /*particleNumber*/,
                                                           const std::vector<char>& /*particle*/ ) {
    return 0;
}

cycle_shapes_strategy::cycle_shapes_strategy()
    : shapeCount( 0 ) {}
cycle_shapes_strategy::cycle_shapes_strategy( boost::int32_t shapeCount )
    : shapeCount( shapeCount ) {}
boost::int32_t cycle_shapes_strategy::get_shape_number( boost::int64_t particleNumber,
                                                        const std::vector<char>& /*particle*/ ) {
    return static_cast<boost::int32_t>( particleNumber % shapeCount );
}

get_random_shape_by_id_strategy::get_random_shape_by_id_strategy()
    : m_seed( 0 )
    , m_count( 0 )
    , m_idAcc( 0 ) {}
get_random_shape_by_id_strategy::get_random_shape_by_id_strategy(
    frantic::channels::channel_const_cvt_accessor<boost::int32_t>& idAcc, boost::int32_t shapeCount,
    boost::int32_t seed )
    : m_seed( seed )
    , m_count( shapeCount )
    , m_idAcc( idAcc ) {}
boost::int32_t get_random_shape_by_id_strategy::get_shape_number( boost::int64_t /*particleNumber*/,
                                                                  const std::vector<char>& particle ) {
    return id_randomized_geometry( m_idAcc( particle ), m_count, m_seed );
}

use_shapeindex_channel_strategy::use_shapeindex_channel_strategy()
    : m_shapeIndexAcc( 0 )
    , m_count( 0 ) {}
use_shapeindex_channel_strategy::use_shapeindex_channel_strategy(
    frantic::channels::channel_const_cvt_accessor<boost::int32_t> shapeIndexAcc, boost::int32_t shapeCount )
    : m_shapeIndexAcc( shapeIndexAcc )
    , m_count( shapeCount ) {}
boost::int32_t use_shapeindex_channel_strategy::get_shape_number( boost::int64_t /*particleNumber*/,
                                                                  const std::vector<char>& particle ) {
    return mod_trunc_neginf( m_shapeIndexAcc( particle ) - 1, m_count );
}

instance_info_accessor::instance_info_accessor( const frantic::channels::channel_map& channelMap,
                                                boost::shared_ptr<frost::geometry_meshing_parameters_interface> params )
    : m_orientationMethod( params->get_geometry_orientation_mode() )
    , m_geometryLookAtPosition( params->get_geometry_look_at_position() )
    , m_coordinateSystem( params->get_coordinate_system() )
    , m_geometryLookAtOrientation( params->get_geometry_look_at_orientation() )
    , m_geometryOrientation( params->get_geometry_orientation() )
    , m_orientationDivergenceAxis( params->get_orientation_divergence_axis() )
    , m_orientationRestrictDivergenceAxis( params->get_orientation_restrict_divergence_axis() )
    , m_geometryOrientationDivergenceAxisSpace( params->get_geometry_orientation_divergence_axis_space() )
    , m_divergence( params->get_geometry_orientation_divergence() ) // maximum rotation
    , m_divergenceX2( m_divergence * 2 ) // Must generate random in 0 to 360 (max) then shift it to -180 to 180
    , m_divergenceRad( boost::math::constants::pi<float>() *
                       frantic::math::clamp<float>( fabs( m_divergence ), 0, 180 ) / 180 )
    , m_idAcc( 0 )
    , m_vecAcc( vector3f( 0, 0, 1 ) )
    , m_orientationAcc( vector4f( 0, 0, 0, 1 ) )
    , m_positionAcc( vector3f( 0 ) )
    , m_radiusAcc( 0 )
    , m_radiusXYZAcc()
    , m_shapeIndexAcc( 0 ) {

    using namespace frost;

    if( channelMap.has_channel( Frost_IDCHANNEL ) )
        m_idAcc = channelMap.get_const_cvt_accessor<boost::int32_t>( Frost_IDCHANNEL );
    if( channelMap.has_channel( Frost_POSITIONCHANNEL ) )
        m_positionAcc = channelMap.get_const_cvt_accessor<vector3f>( Frost_POSITIONCHANNEL );
    if( channelMap.has_channel( Frost_RADIUSCHANNEL ) )
        m_radiusAcc = channelMap.get_const_cvt_accessor<float>( Frost_RADIUSCHANNEL );
    // Temporarily disable RadiusXYZ until we have proper support in the Frost UI
    /*
    if( channelMap.has_channel( Frost_RADIUSXYZCHANNEL ) )
     m_radiusXYZAcc = channelMap.get_const_cvt_accessor<vector3f>( Frost_RADIUSXYZCHANNEL );
    */

    const frantic::tstring orientationVector = params->get_geometry_orientation_vector_channel();

    switch( m_orientationMethod ) {
    case geometry_orientation_mode::look_at:
        break;
    case geometry_orientation_mode::match_object:
        break;
    case geometry_orientation_mode::specify:
        break;
    case geometry_orientation_mode::orientation_channel:
        if( channelMap.has_channel( _T("Orientation") ) )
            m_orientationAcc = channelMap.get_const_cvt_accessor<vector4f>( _T("Orientation") );
        else
            throw std::runtime_error( "The particles don\'t have an \'Orientation\' channel, which was requested for "
                                      "instance geometry orientation." );
        break;
    case geometry_orientation_mode::vector_channel:
        if( channelMap.has_channel( orientationVector ) )
            m_vecAcc = channelMap.get_const_cvt_accessor<vector3f>( orientationVector );
        else
            throw std::runtime_error( "The particles don\'t have a \'" +
                                      frantic::strings::to_string( orientationVector ) +
                                      "\' channel, which was requested for instance geometry orientation." );
        break;
    default:
        throw std::runtime_error( "Unrecognized Instance Geometry Orientation Type: " +
                                  boost::lexical_cast<std::string>( m_orientationMethod ) );
    }
}

transform4f instance_info_accessor::get_transform( const std::vector<char>& particle ) const {
    using namespace frost;

    const boost::int32_t id = m_idAcc( particle );

    // Orientation
    transform4f rotationTransform = transform4f::identity(); // set to identity in case invalid value is found
    switch( m_orientationMethod ) {
    case geometry_orientation_mode::look_at: {
        const frantic::graphics::vector3f direction( m_geometryLookAtPosition - m_positionAcc( particle ) );
        if( !direction.is_zero() ) {
            switch( m_coordinateSystem ) {
            case frantic::graphics::coordinate_system::right_handed_yup:
                rotationTransform = transform4f::from_normal_y( normalized( direction ) );
                break;
            case frantic::graphics::coordinate_system::right_handed_zup:
            default:
                rotationTransform = transform4f::from_normal_z( normalized( direction ) );
                break;
            }
        }
    } break;
    case geometry_orientation_mode::match_object: {
        rotationTransform = transform4f::from_euler_xyz_rotation( m_geometryLookAtOrientation );
    } break;
    case geometry_orientation_mode::orientation_channel: {
        const vector4f q = m_orientationAcc.get( particle );
        if( !q.is_zero() )
            rotationTransform = q.quaternion_to_matrix();
        // vector4f::quaternion_to_matrix() currently returns a zero matrix if the quaternion's magnitude is too small.
        // In this case, we reset to identity.
        if( rotationTransform.is_zero() )
            rotationTransform.set_to_identity();
    } break;
    case geometry_orientation_mode::vector_channel: {
        const vector3f v = m_vecAcc.get( particle );
        if( v != vector3f( 0 ) ) {
            vector3f vn = v.to_normalized();
            for( int axis = 0; axis < 3; ++axis ) {
                if( fabs( vn[axis] ) <= std::numeric_limits<float>::epsilon() ) {
                    vn[axis] = 0;
                }
            }
            if( vn.get_magnitude_squared() < 0.99f ) {
                switch( m_coordinateSystem ) {
                case frantic::graphics::coordinate_system::right_handed_yup:
                    vn = vector3f::from_yaxis();
                    break;
                case frantic::graphics::coordinate_system::right_handed_zup:
                default:
                    vn = vector3f::from_zaxis();
                    break;
                }
            }
            switch( m_coordinateSystem ) {
            case frantic::graphics::coordinate_system::right_handed_yup:
                rotationTransform = transform4f::from_normal_y( vn );
                break;
            case frantic::graphics::coordinate_system::right_handed_zup:
            default:
                rotationTransform = transform4f::from_normal_z( vn );
                break;
            }
        }
    } break;
    case geometry_orientation_mode::specify:
        rotationTransform = transform4f::from_euler_xyz_rotation( m_geometryOrientation );
        break;
    }

    transform4f divergenceTransform;
    const vector3f divergenceAxis( m_orientationDivergenceAxis );

    if( m_orientationRestrictDivergenceAxis ) {
        // TODO: give option of object space or world space axis
        const float angleRad = id_randomized_rotation( id, 2 * m_divergenceRad, 177720997 ) - m_divergenceRad;
        divergenceTransform = transform4f::from_angle_axis( angleRad, divergenceAxis );
    } else {
        // Set the random divergence (-180 to 180 Max)
        // amounts to rotate in x, y, z direction
        float orientationX = id_randomized_rotation( id, m_divergenceX2, 1423264495 ) - m_divergence;
        float orientationY = id_randomized_rotation( id, m_divergenceX2, 778162823 ) - m_divergence;
        float orientationZ = id_randomized_rotation( id, m_divergenceX2, 171472941 ) - m_divergence;
        divergenceTransform =
            transform4f::from_euler_xyz_rotation_degrees( vector3f( orientationX, orientationY, orientationZ ) );
    }

    const geometry_orientation_divergence_axis_space::option divergenceAxisSpace =
        m_geometryOrientationDivergenceAxisSpace;
    if( divergenceAxisSpace != geometry_orientation_divergence_axis_space::world &&
        divergenceAxisSpace != geometry_orientation_divergence_axis_space::local ) {
        throw std::runtime_error( "Unrecognized Instance Geometry Divergence Axis Space: " +
                                  boost::lexical_cast<std::string>( divergenceAxisSpace ) );
    }
    const transform4f rotationAndDivergenceTransform =
        ( divergenceAxisSpace == geometry_orientation_divergence_axis_space::local )
            ? rotationTransform * divergenceTransform
            : divergenceTransform * rotationTransform;

    frantic::graphics::transform4f xform;
    if( m_radiusXYZAcc.is_valid() ) {
        const vector3f scaleXYZ = m_radiusXYZAcc( particle );
        xform = transform4f::from_translation( m_positionAcc( particle ) ) * transform4f::from_scale( scaleXYZ ) *
                rotationAndDivergenceTransform;
    } else {
        const float scale = m_radiusAcc( particle );
        xform = transform4f::from_translation( m_positionAcc( particle ) ) * transform4f::from_scale( scale ) *
                rotationAndDivergenceTransform;
    }

    return xform;
}

void prepare_mesh_material_channel( frantic::geometry::trimesh3& mesh, int materialMode, int shapeNumber,
                                    const std::vector<int>& meshEndMaterialId,
                                    const std::vector<geometry_material_id_map_t>& geometryMaterialIdMaps,
                                    const boost::uint16_t undefinedMaterialID ) {
    typedef boost::uint16_t material_id_t;

    switch( materialMode ) {
    case frost::material_mode::single:
        if( mesh.has_face_channel( Frost_MESH_MTLINDEX_CHANNEL_NAME ) ) {
            mesh.erase_face_channel( Frost_MESH_MTLINDEX_CHANNEL_NAME );
        }
        break;
    case frost::material_mode::mtlindex_channel: {
        if( !mesh.has_face_channel( Frost_MESH_MTLINDEX_CHANNEL_NAME ) ) {
            mesh.add_face_channel<material_id_t>( Frost_MESH_MTLINDEX_CHANNEL_NAME );
        }
        frantic::geometry::trimesh3_face_channel_general_accessor acc =
            mesh.get_face_channel_general_accessor( Frost_MESH_MTLINDEX_CHANNEL_NAME );
        for( std::size_t faceNumber = 0; faceNumber < acc.size(); ++faceNumber ) {
            memset( acc.data( faceNumber ), 0, acc.primitive_size() );
        }
    } break;
    case frost::material_mode::shape_number: {
        if( mesh.has_face_channel( Frost_MESH_MTLINDEX_CHANNEL_NAME ) ) {
            frantic::geometry::trimesh3_face_channel_general_accessor acc =
                mesh.get_face_channel_general_accessor( Frost_MESH_MTLINDEX_CHANNEL_NAME );
            if( acc.arity() != frantic::channels::channel_data_type_traits<material_id_t>::arity() ||
                acc.data_type() != frantic::channels::channel_data_type_traits<material_id_t>::data_type() ) {
                mesh.erase_face_channel( Frost_MESH_MTLINDEX_CHANNEL_NAME );
            }
        }
        if( !mesh.has_face_channel( Frost_MESH_MTLINDEX_CHANNEL_NAME ) ) {
            mesh.add_face_channel<material_id_t>( Frost_MESH_MTLINDEX_CHANNEL_NAME );
        }
        frantic::geometry::trimesh3_face_channel_accessor<material_id_t> acc =
            mesh.get_face_channel_accessor<material_id_t>( Frost_MESH_MTLINDEX_CHANNEL_NAME );
        const material_id_t mtlId = static_cast<material_id_t>(
            shapeNumber % ( static_cast<int>( std::numeric_limits<material_id_t>::max() ) + 1 ) );
        for( std::size_t faceNumber = 0; faceNumber < acc.size(); ++faceNumber ) {
            acc[faceNumber] = mtlId;
        }
    } break;
    case frost::material_mode::material_id_from_geometry:
        break;
    case frost::material_mode::material_from_geometry: {
        if( !mesh.has_face_channel( Frost_MESH_MTLINDEX_CHANNEL_NAME ) ) {
            mesh.add_face_channel<material_id_t>( Frost_MESH_MTLINDEX_CHANNEL_NAME );
            frantic::geometry::trimesh3_face_channel_accessor<material_id_t> acc =
                mesh.get_face_channel_accessor<material_id_t>( Frost_MESH_MTLINDEX_CHANNEL_NAME );
            for( std::size_t faceNumber = 0; faceNumber < acc.size(); ++faceNumber ) {
                acc[faceNumber] = undefinedMaterialID;
            }
            return;
        }
        frantic::geometry::trimesh3_face_channel_cvt_accessor<material_id_t> acc =
            mesh.get_face_channel_cvt_accessor<material_id_t>( Frost_MESH_MTLINDEX_CHANNEL_NAME );
        const int endMtlId = meshEndMaterialId[shapeNumber];
        if( shapeNumber < static_cast<int>( geometryMaterialIdMaps.size() ) ) {
            const geometry_material_id_map_t& materialIdMap = geometryMaterialIdMaps[shapeNumber];
            for( std::size_t faceNumber = 0; faceNumber < acc.size(); ++faceNumber ) {
                material_id_t mtlId = acc.get( faceNumber );
                if( mtlId >= endMtlId ) {
                    mtlId = static_cast<material_id_t>( mtlId % endMtlId );
                }
                geometry_material_id_map_t::const_iterator i = materialIdMap.find( mtlId );
                if( i == materialIdMap.end() ) {
                    acc.set( faceNumber, undefinedMaterialID );
                } else {
                    acc.set( faceNumber, i->second );
                }
            }
        } else {
            for( std::size_t faceNumber = 0; faceNumber < acc.size(); ++faceNumber ) {
                acc.set( faceNumber, undefinedMaterialID );
            }
        }
    } break;
    default:
        throw std::runtime_error( "Unrecognized Material Mode: " + boost::lexical_cast<std::string>( materialMode ) );
    }
}

template <class Derived>
class build_mesh_geometry_factory : public instance_geometry_factory {
    mesh_and_accessors m_mesh;
    std::size_t m_denylistChannel;
    build_mesh_geometry_factory();
    const frantic::graphics::transform4f m_coordinateSystem;

  public:
    build_mesh_geometry_factory( const frantic::channels::channel_map& vertexChannelMap,
                                 const frantic::channels::channel_map& faceChannelMap,
                                 const material_id_mode_info& materialModeInfo,
                                 frost::hard_edge_type::option hardEdgeType,
                                 const frantic::graphics::transform4f& coordinateSystem )
        : m_coordinateSystem( coordinateSystem ) {
        Derived::build_mesh( m_mesh.mesh, vertexChannelMap, faceChannelMap, hardEdgeType );
        get_channel_accessors( m_mesh.mesh, vertexChannelMap, faceChannelMap, m_mesh.vertexChannels,
                               m_mesh.faceChannels );
        prepare_mesh_material_channel( m_mesh.mesh, materialModeInfo.materialMode, 0,
                                       materialModeInfo.meshEndMaterialId, materialModeInfo.geometryMaterialIdMaps,
                                       materialModeInfo.undefinedMaterialID );
        m_denylistChannel = std::numeric_limits<std::size_t>::max();
        for( size_t i = 0; i < vertexChannelMap.channel_count(); ++i ) {
            if( vertexChannelMap[i].name() == Frost_TEXTURECOORDCHANNEL ) {
                m_denylistChannel = i;
            }
        }

        m_mesh.mesh.transform( m_coordinateSystem );
    }
    std::size_t get_keep_vertex_channel_number() { return m_denylistChannel; }
    mesh_and_accessors* get_instance_geometry( double /*timeInSeconds*/, std::size_t /*shapeNumber*/ ) {
        return &m_mesh;
    }
    std::size_t get_shape_count() const { return 1; }
    bool is_animated() const { return false; }
};

class plane_geometry_factory : public build_mesh_geometry_factory<plane_geometry_factory> {
  public:
    plane_geometry_factory( const frantic::channels::channel_map& vertexChannelMap,
                            const frantic::channels::channel_map& faceChannelMap,
                            const material_id_mode_info& materialModeInfo, frost::hard_edge_type::option hardEdgeType,
                            const frantic::graphics::transform4f& coordinateSystem )
        : build_mesh_geometry_factory<plane_geometry_factory>( vertexChannelMap, faceChannelMap, materialModeInfo,
                                                               hardEdgeType, coordinateSystem ) {}

    static void build_mesh( frantic::geometry::trimesh3& mesh, const frantic::channels::channel_map& vertexChannelMap,
                            const frantic::channels::channel_map& faceChannelMap,
                            frost::hard_edge_type::option hardEdgeType ) {
        build_plane_mesh( mesh, vertexChannelMap, faceChannelMap, hardEdgeType );
    }
};

class sprite_geometry_factory : public build_mesh_geometry_factory<sprite_geometry_factory> {
  public:
    sprite_geometry_factory( const frantic::channels::channel_map& vertexChannelMap,
                             const frantic::channels::channel_map& faceChannelMap,
                             const material_id_mode_info& materialModeInfo, frost::hard_edge_type::option hardEdgeType,
                             const frantic::graphics::transform4f& coordinateSystem )
        : build_mesh_geometry_factory<sprite_geometry_factory>( vertexChannelMap, faceChannelMap, materialModeInfo,
                                                                hardEdgeType, coordinateSystem ) {}

    static void build_mesh( frantic::geometry::trimesh3& mesh, const frantic::channels::channel_map& vertexChannelMap,
                            const frantic::channels::channel_map& faceChannelMap,
                            frost::hard_edge_type::option hardEdgeType ) {
        build_sprite_mesh( mesh, vertexChannelMap, faceChannelMap, hardEdgeType );
    }
};

class tetrahedron_geometry_factory : public build_mesh_geometry_factory<tetrahedron_geometry_factory> {
  public:
    tetrahedron_geometry_factory( const frantic::channels::channel_map& vertexChannelMap,
                                  const frantic::channels::channel_map& faceChannelMap,
                                  const material_id_mode_info& materialModeInfo,
                                  frost::hard_edge_type::option hardEdgeType,
                                  const frantic::graphics::transform4f& coordinateSystem )
        : build_mesh_geometry_factory<tetrahedron_geometry_factory>( vertexChannelMap, faceChannelMap, materialModeInfo,
                                                                     hardEdgeType, coordinateSystem ) {}

    static void build_mesh( frantic::geometry::trimesh3& mesh, const frantic::channels::channel_map& vertexChannelMap,
                            const frantic::channels::channel_map& faceChannelMap,
                            frost::hard_edge_type::option hardEdgeType ) {
        build_tetrahedron_mesh( mesh, vertexChannelMap, faceChannelMap, hardEdgeType );
    }
};

class box_geometry_factory : public build_mesh_geometry_factory<box_geometry_factory> {
  public:
    box_geometry_factory( const frantic::channels::channel_map& vertexChannelMap,
                          const frantic::channels::channel_map& faceChannelMap,
                          const material_id_mode_info& materialModeInfo, frost::hard_edge_type::option hardEdgeType,
                          const frantic::graphics::transform4f& coordinateSystem )
        : build_mesh_geometry_factory<box_geometry_factory>( vertexChannelMap, faceChannelMap, materialModeInfo,
                                                             hardEdgeType, coordinateSystem ) {}

    static void build_mesh( frantic::geometry::trimesh3& mesh, const frantic::channels::channel_map& vertexChannelMap,
                            const frantic::channels::channel_map& faceChannelMap,
                            frost::hard_edge_type::option hardEdgeType ) {
        build_box_mesh( mesh, vertexChannelMap, faceChannelMap, hardEdgeType );
    }
};

class sphere20_geometry_factory : public build_mesh_geometry_factory<sphere20_geometry_factory> {
  public:
    sphere20_geometry_factory( const frantic::channels::channel_map& vertexChannelMap,
                               const frantic::channels::channel_map& faceChannelMap,
                               const material_id_mode_info& materialModeInfo,
                               frost::hard_edge_type::option hardEdgeType,
                               const frantic::graphics::transform4f& coordinateSystem )
        : build_mesh_geometry_factory<sphere20_geometry_factory>( vertexChannelMap, faceChannelMap, materialModeInfo,
                                                                  hardEdgeType, coordinateSystem ) {}

    static void build_mesh( frantic::geometry::trimesh3& mesh, const frantic::channels::channel_map& vertexChannelMap,
                            const frantic::channels::channel_map& faceChannelMap,
                            frost::hard_edge_type::option hardEdgeType ) {
        build_sphere_mesh( mesh, vertexChannelMap, faceChannelMap, hardEdgeType );
    }
};

class custom_geometry_factory : public instance_geometry_factory {
    typedef std::map<double, mesh_and_accessors*> mesh_sample_map_t;

    std::vector<boost::shared_ptr<mesh_and_accessors>> m_storage;
    std::vector<mesh_sample_map_t> m_samples;

    boost::shared_ptr<geometry_source_interface> m_geometrySource;

    frantic::channels::channel_map m_vertexChannelMap;
    frantic::channels::channel_map m_faceChannelMap;

    material_id_mode_info m_materialModeInfo;

    std::size_t m_denylistChannel;

    custom_geometry_factory& operator=( const custom_geometry_factory& ); // not implemented

  public:
    custom_geometry_factory( boost::shared_ptr<geometry_source_interface> geometrySource,
                             const frantic::channels::channel_map& vertexChannelMap,
                             const frantic::channels::channel_map& faceChannelMap,
                             const material_id_mode_info& materialModeInfo )
        : m_vertexChannelMap( vertexChannelMap )
        , m_faceChannelMap( faceChannelMap )
        , m_materialModeInfo( materialModeInfo )
        , m_geometrySource( geometrySource )
        , m_denylistChannel( std::numeric_limits<std::size_t>::max() ) {
        m_samples.resize( m_geometrySource->get_geometry_count() );
    }
    std::size_t get_keep_vertex_channel_number() { return m_denylistChannel; }
    mesh_and_accessors* get_instance_geometry( double timeInSeconds, std::size_t shapeNumber ) {
        if( shapeNumber >= m_samples.size() ) {
            return 0;
        }

        mesh_sample_map_t::iterator i = m_samples[shapeNumber].find( timeInSeconds );
        if( i == m_samples[shapeNumber].end() ) {
            boost::shared_ptr<mesh_and_accessors> meshAndAccessors( new mesh_and_accessors() );
            bool gotMesh = m_geometrySource->get_geometry( meshAndAccessors->mesh, timeInSeconds, shapeNumber );
            if( gotMesh ) {
                conform_mesh_channel_types( meshAndAccessors->mesh, m_vertexChannelMap, m_denylistChannel,
                                            m_faceChannelMap );
                prepare_mesh_material_channel(
                    meshAndAccessors->mesh, m_materialModeInfo.materialMode, static_cast<boost::int32_t>( shapeNumber ),
                    m_materialModeInfo.meshEndMaterialId, m_materialModeInfo.geometryMaterialIdMaps,
                    m_materialModeInfo.undefinedMaterialID );
                get_channel_accessors( meshAndAccessors->mesh, m_vertexChannelMap, m_faceChannelMap,
                                       meshAndAccessors->vertexChannels, meshAndAccessors->faceChannels );
                m_samples[shapeNumber][timeInSeconds] = meshAndAccessors.get();
                m_storage.push_back( meshAndAccessors );
                return meshAndAccessors.get();
            } else {
                return 0;
            }
        } else {
            return i->second;
        }
    }
    std::size_t get_shape_count() const { return m_samples.size(); }
    bool is_animated() const { return true; }
};

std::size_t get_instance_geometry_shape_count( frost::geometry_type::option geometryType,
                                               boost::shared_ptr<geometry_source_interface> customGeometrySource ) {
    switch( geometryType ) {
    case frost::geometry_type::plane:
        return 1;
    case frost::geometry_type::sprite:
        return 1;
    case frost::geometry_type::tetrahedron:
        return 1;
    case frost::geometry_type::box:
        return 1;
    case frost::geometry_type::sphere20:
        return 1;
    case frost::geometry_type::custom_geometry:
        return customGeometrySource->get_geometry_count();
    default:
        throw std::runtime_error( "get_instance_geometry_shape_count Error: unrecognized instance geometry source: " +
                                  boost::lexical_cast<std::string>( geometryType ) );
    }
}

boost::shared_ptr<instance_geometry_factory> get_instance_geometry_factory(
    frost::geometry_type::option geometryType, const frantic::channels::channel_map& vertexChannelMap,
    const frantic::channels::channel_map& faceChannelMap, const material_id_mode_info& materialModeInfo,
    boost::shared_ptr<geometry_source_interface> customGeometrySource, frost::hard_edge_type::option hardEdgeType,
    const frantic::graphics::transform4f& coordinateSystem ) {
    switch( geometryType ) {
    case frost::geometry_type::plane:
        return boost::shared_ptr<instance_geometry_factory>( new plane_geometry_factory(
            vertexChannelMap, faceChannelMap, materialModeInfo, hardEdgeType, coordinateSystem ) );
    case frost::geometry_type::sprite:
        return boost::shared_ptr<instance_geometry_factory>( new sprite_geometry_factory(
            vertexChannelMap, faceChannelMap, materialModeInfo, hardEdgeType, coordinateSystem ) );
    case frost::geometry_type::tetrahedron:
        return boost::shared_ptr<instance_geometry_factory>( new tetrahedron_geometry_factory(
            vertexChannelMap, faceChannelMap, materialModeInfo, hardEdgeType, coordinateSystem ) );
    case frost::geometry_type::box:
        return boost::shared_ptr<instance_geometry_factory>( new box_geometry_factory(
            vertexChannelMap, faceChannelMap, materialModeInfo, hardEdgeType, coordinateSystem ) );
    case frost::geometry_type::sphere20:
        return boost::shared_ptr<instance_geometry_factory>( new sphere20_geometry_factory(
            vertexChannelMap, faceChannelMap, materialModeInfo, hardEdgeType, coordinateSystem ) );
    case frost::geometry_type::custom_geometry:
        return boost::shared_ptr<instance_geometry_factory>(
            new custom_geometry_factory( customGeometrySource, vertexChannelMap, faceChannelMap, materialModeInfo ) );
    default:
        throw std::runtime_error( "get_instance_geometry_factory Error: unrecognized instance geometry source: " +
                                  boost::lexical_cast<std::string>( geometryType ) );
    }
}

void geometry_convert_particles_to_trimesh3( boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                                             const frantic::channels::channel_propagation_policy& cpp,
                                             boost::shared_ptr<frost::geometry_meshing_parameters_interface> params,
                                             boost::shared_ptr<geometry_source_interface> sceneGeometrySource,
                                             material_id_mode_info& materialModeInfo,
                                             frantic::geometry::trimesh3& outMesh, boost::int64_t& outParticleCount,
                                             frantic::logging::progress_logger& progressLogger ) {
    using namespace frost;
    using namespace frantic::graphics;
    using namespace frantic::geometry;
    using namespace frantic::channels;
    using std::vector;

    outMesh.clear();
    outParticleCount = 0;

    std::vector<char> particle( pin->particle_size() );

    std::vector<frantic::tstring> particleChannelNames;
    for( std::size_t i = 0; i < pin->get_channel_map().channel_count(); ++i ) {
        const frantic::channels::channel ch = pin->get_channel_map()[i];
        const frantic::tstring& channelName = ch.name();

        if( cpp.is_channel_included( channelName ) ) {
            particleChannelNames.push_back( channelName );
        }
    }

    // Create accessors and channel converters for the particle input and the mesh output.
    vector<channel_general_accessor> inputChannels;
    for( unsigned i = 0; i < particleChannelNames.size(); ++i ) {
        inputChannels.push_back( pin->get_channel_map().get_general_accessor( particleChannelNames[i] ) );
    }

    typedef std::map<frantic::tstring, frantic::tstring> channel_name_map_t;
    channel_name_map_t particleToMeshChannelNames;
    particleToMeshChannelNames[_T("MtlIndex")] = _T("MaterialID");

    frantic::channels::channel_map vertexChannelMap;
    frantic::channels::channel_map faceChannelMap;

    vector<channel_general_accessor> particleVertexChannels;
    vector<channel_general_accessor> particleFaceChannels;

    for( std::size_t particleChannelNumber = 0; particleChannelNumber < particleChannelNames.size();
         ++particleChannelNumber ) {
        const frantic::tstring& particleChannelName = particleChannelNames[particleChannelNumber];
        const channel_general_accessor& ch = inputChannels[particleChannelNumber];

        channel_name_map_t::const_iterator particleToMeshChannelName =
            particleToMeshChannelNames.find( particleChannelName );

        if( particleToMeshChannelName != particleToMeshChannelNames.end() ) {
            faceChannelMap.define_channel( particleToMeshChannelName->second, ch.arity(), ch.data_type() );
            particleFaceChannels.push_back( ch );
        } else {
            vertexChannelMap.define_channel( particleChannelName, ch.arity(), ch.data_type() );
            particleVertexChannels.push_back( ch );
        }
    }

    vertexChannelMap.end_channel_definition( 4, true );
    faceChannelMap.end_channel_definition( 4, true );

    const geometry_type::option instanceGeometryType = params->get_geometry_type();
    const hard_edge_type::option hardEdgeType = params->get_hard_edge_type();
    frantic::graphics::transform4f coordinateSystemXForm;
    // the transform is created with an assumed right_handed_zup coordinate system, so we need to convert it to use the
    // native coordinate system
    frantic::graphics::coordinate_system::create_transform( coordinateSystemXForm,
                                                            frantic::graphics::coordinate_system::right_handed_zup,
                                                            params->get_coordinate_system() );

    boost::shared_ptr<instance_geometry_factory> instanceGeometryFactory =
        get_instance_geometry_factory( instanceGeometryType, vertexChannelMap, faceChannelMap, materialModeInfo,
                                       sceneGeometrySource, hardEdgeType, coordinateSystemXForm );

    if( !instanceGeometryFactory ) {
        throw std::runtime_error( "Internal Error: instance geometry factory is NULL" );
    }

    const std::size_t denylistChannel =
        instanceGeometryFactory->get_keep_vertex_channel_number(); // do not copy this channel number from the particles

    const std::size_t shapeCount = instanceGeometryFactory->get_shape_count();

    if( shapeCount > 0 ) {

        // TODO: Reserve the average number of vertices and faces.
        // This should reserve the vertex and face channels too.

        // Generate the mesh
        const bool isOrientationRestricted = params->get_orientation_restrict_divergence_axis();
        const boost::uint32_t randomSelectionSeed = params->get_geometry_selection_seed();

        frantic::channels::channel_map pinChannelMap = pin->get_channel_map();
        frantic::channels::channel_const_cvt_accessor<boost::int32_t> idAcc( 0 );
        if( pinChannelMap.has_channel( Frost_IDCHANNEL ) )
            idAcc = pinChannelMap.get_const_cvt_accessor<boost::int32_t>( Frost_IDCHANNEL );
        frantic::channels::channel_const_cvt_accessor<float> radiusAcc( 0 );
        if( pinChannelMap.has_channel( Frost_RADIUSCHANNEL ) )
            radiusAcc = pinChannelMap.get_const_cvt_accessor<float>( Frost_RADIUSCHANNEL );
        frantic::channels::channel_const_cvt_accessor<vector3f> radiusXYZAcc;
        // Temporarily disable RadiusXYZ until we have proper support in the Frost UI
        /*
        if( pinChannelMap.has_channel( Frost_RADIUSXYZCHANNEL ) )
         radiusXYZAcc = pinChannelMap.get_const_cvt_accessor<vector3f>( Frost_RADIUSXYZCHANNEL );
        */
        frantic::channels::channel_const_cvt_accessor<boost::int32_t> shapeIndexAcc( 0 );

        boost::scoped_ptr<get_shape_number_strategy> randomSelectionStrategy;
        if( instanceGeometryFactory->get_shape_count() > 1 ) {
            const geometry_selection_mode::option randomSelectionModeVal = params->get_geometry_selection_mode();
            switch( randomSelectionModeVal ) {
            case geometry_selection_mode::cycle:
                randomSelectionStrategy.reset( new cycle_shapes_strategy( static_cast<boost::int32_t>( shapeCount ) ) );
                break;
            case geometry_selection_mode::random_by_id:
                randomSelectionStrategy.reset( new get_random_shape_by_id_strategy(
                    idAcc, static_cast<boost::int32_t>( shapeCount ), randomSelectionSeed ) );
                break;
            case geometry_selection_mode::shapeindex_channel:
                if( !pinChannelMap.has_channel( Frost_SHAPEINDEXCHANNEL ) ) {
                    throw std::runtime_error( "The particles don\'t have a \"ShapeIndex\" channel, which is required "
                                              "when selecting instance geometry by ShapeIndex." );
                }
                shapeIndexAcc = pinChannelMap.get_const_cvt_accessor<boost::int32_t>( Frost_SHAPEINDEXCHANNEL );
                randomSelectionStrategy.reset(
                    new use_shapeindex_channel_strategy( shapeIndexAcc, static_cast<boost::int32_t>( shapeCount ) ) );
                break;
            default:
                throw std::runtime_error( "Unrecognized Instance Geometry Random Selection Mode: " +
                                          boost::lexical_cast<std::string>( randomSelectionModeVal ) );
            }
        } else {
            randomSelectionStrategy.reset( new use_first_shape_strategy() );
        }
        if( randomSelectionStrategy == 0 ) {
            throw std::runtime_error( "Internal Error: No Instance Geometry Selection Mode." );
        }

        bool useAbsoluteTimeChannel = false;
        frantic::channels::channel_const_cvt_accessor<float> absoluteTimeAcc;

        bool useTimeOffsetChannel = false;
        frantic::channels::channel_const_cvt_accessor<float> timeOffsetAcc;

        bool useGeomTimeChannel = false;
        frantic::channels::channel_const_cvt_accessor<float> geomTimeAcc;

        bool randomizeTimeOffsetByID = false;
        double maxRandomizedTimeOffsetInSeconds = 0;
        boost::uint32_t geometrySampleTimeSeed = 0;

        double baseTimeInSeconds = params->get_time();
        const geometry_sample_time_base_mode::option baseTimeMode = params->get_geometry_sample_time_base_mode();
        switch( baseTimeMode ) {
        case geometry_sample_time_base_mode::time_0:
            baseTimeInSeconds = 0;
            break;
        case geometry_sample_time_base_mode::current_time:
            baseTimeInSeconds = params->get_time();
            break;
        default:
            throw std::runtime_error( "unrecognized geometrySampleTimeBaseMode: " +
                                      boost::lexical_cast<std::string>( baseTimeMode ) );
        }

        if( instanceGeometryFactory->is_animated() ) {
            const int instanceGeometryAnimationMode = params->get_geometry_sample_time_offset_mode();
            switch( instanceGeometryAnimationMode ) {
            case geometry_sample_time_offset_mode::no_offset:
                break;
            case geometry_sample_time_offset_mode::random_by_id:
                maxRandomizedTimeOffsetInSeconds = params->get_geometry_sample_time_max_random_offset();
                randomizeTimeOffsetByID = true;
                geometrySampleTimeSeed = params->get_geometry_sample_time_seed();
                break;
            case geometry_sample_time_offset_mode::abstime_channel:
                if( !pinChannelMap.has_channel( _T("AbsTime") ) ) {
                    throw std::runtime_error( "The particles don\'t have an \"AbsTime\" channel, which is required "
                                              "when animating geometry by AbsTime channel." );
                }
                absoluteTimeAcc = pinChannelMap.get_const_cvt_accessor<float>( _T("AbsTime") );
                useAbsoluteTimeChannel = true;
                break;
            case geometry_sample_time_offset_mode::timeoffset_channel:
                if( !pinChannelMap.has_channel( _T("TimeOffset") ) ) {
                    throw std::runtime_error( "The particles don\'t have a \"TimeOffset\" channel, which is required "
                                              "when animating geometry by TimeOffset channel." );
                }
                timeOffsetAcc = pinChannelMap.get_const_cvt_accessor<float>( _T("TimeOffset") );
                useTimeOffsetChannel = true;
                break;
            case geometry_sample_time_offset_mode::geomtime_channel:
                if( !pinChannelMap.has_channel( _T("GeomTime") ) ) {
                    throw std::runtime_error( "The particles don\'t have a \"GeomTime\" channel, which is required "
                                              "when animating geometry by GeomTime channel." );
                }
                geomTimeAcc = pinChannelMap.get_const_cvt_accessor<float>( _T("GeomTime") );
                useGeomTimeChannel = true;
                break;
            default:
                throw std::runtime_error( "Unrecognized Particle Geometry Animation Timing Mode: " +
                                          boost::lexical_cast<std::string>( instanceGeometryAnimationMode ) );
            }
        }

        // Create vector of channels that should be transformed. Velocity should not be transformed.
        std::vector<std::pair<frantic::tstring, frantic::geometry::trimesh3::vector_type>> transformNoChannels;
        std::vector<std::pair<frantic::tstring, frantic::geometry::trimesh3::vector_type>> transformNormalChannel;
        transformNormalChannel.push_back( std::pair<frantic::tstring, frantic::geometry::trimesh3::vector_type>(
            _T("Normal"), frantic::geometry::trimesh3::NORMAL ) );

        // update progress on the main particle loop
        long long total = pin->particle_count();
        long long count = 0;

        while( pin->get_particle( particle ) ) {

            const boost::int32_t id = idAcc( particle );

            // Random Shape Selection
            boost::int32_t currentElement =
                randomSelectionStrategy->get_shape_number( pin->particle_index(), particle );

            double requestTimeInSeconds = baseTimeInSeconds;
            if( useAbsoluteTimeChannel ) {
                requestTimeInSeconds = params->frames_to_seconds( absoluteTimeAcc( particle ) );
            }
            if( useTimeOffsetChannel ) {
                requestTimeInSeconds += params->frames_to_seconds( timeOffsetAcc( particle ) );
            }
            if( useGeomTimeChannel ) {
                requestTimeInSeconds += params->frames_to_seconds( geomTimeAcc( particle ) );
            }
            if( randomizeTimeOffsetByID ) {
                requestTimeInSeconds +=
                    get_unit_random_from_id( id, geometrySampleTimeSeed ) * maxRandomizedTimeOffsetInSeconds;
            }

            if( currentElement < 0 || static_cast<std::size_t>( currentElement ) >= shapeCount ) {
                throw std::runtime_error( "Internal Error: Random shape selection is out of range." );
            }

            mesh_and_accessors* pMeshAndAccessors =
                instanceGeometryFactory->get_instance_geometry( requestTimeInSeconds, currentElement );
            if( pMeshAndAccessors ) {
                { // scope for vertex accessors
                    std::vector<trimesh3_vertex_channel_general_accessor>& particleMeshAccessors =
                        pMeshAndAccessors->vertexChannels;
                    // Copy all the channels to the primitive mesh
                    for( unsigned i = 0; i < particleVertexChannels.size(); ++i ) {
                        if( i == denylistChannel ) {
                            continue;
                        }
                        const void* inputChannelValue = particleVertexChannels[i].get_channel_data_pointer( particle );
                        size_t channelSize = particleVertexChannels[i].primitive_size();
                        trimesh3_vertex_channel_general_accessor& meshChannel = particleMeshAccessors[i];
                        for( unsigned j = 0; j < meshChannel.size(); ++j ) {
                            memcpy( meshChannel.data( j ), inputChannelValue, channelSize );
                        }
                    }
                }
                { // scope for face accessors
                    std::vector<trimesh3_face_channel_general_accessor>& particleMeshAccessors =
                        pMeshAndAccessors->faceChannels;
                    // Copy all the channels to the primitive mesh
                    for( unsigned i = 0; i < particleFaceChannels.size(); ++i ) {
                        const void* inputChannelValue = particleFaceChannels[i].get_channel_data_pointer( particle );
                        size_t channelSize = particleFaceChannels[i].primitive_size();
                        trimesh3_face_channel_general_accessor& meshChannel = particleMeshAccessors[i];
                        for( unsigned j = 0; j < meshChannel.size(); ++j ) {
                            memcpy( meshChannel.data( j ), inputChannelValue, channelSize );
                        }
                    }
                }

                const instance_info_accessor instanceInfoAccessor( pin->get_channel_map(), params );
                const frantic::graphics::transform4f xform = instanceInfoAccessor.get_transform( particle );
                // Transforming the Normal channel requires the inverse transform matrix.
                // If scale is 0 for any component, the transform matrix is not invertible, so we want to avoid
                // transforming the Normal channel.
                bool hasValidNormalTransform = false;
                if( radiusXYZAcc.is_valid() ) {
                    const frantic::graphics::vector3f radiusXYZ = radiusXYZAcc( particle );
                    hasValidNormalTransform = ( radiusXYZ.x != 0.f && radiusXYZ.y != 0.f && radiusXYZ.z != 0.f );
                } else {
                    hasValidNormalTransform = ( radiusAcc( particle ) != 0 );
                }

                outMesh.combine( xform, pMeshAndAccessors->mesh,
                                 hasValidNormalTransform ? transformNormalChannel : transformNoChannels );
            }

            ++count;
            ++outParticleCount;
            if( total > 0 ) {
                progressLogger.update_progress( count, total );
            } else {
                progressLogger.update_progress( 1, 2 ); // TODO: what ?
            }
        }
    } else {
        // No geometry available for particles.
        outMesh.clear();
    }
}
