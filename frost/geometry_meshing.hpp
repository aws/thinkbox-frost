// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <frantic/channels/channel_accessor.hpp>
#include <frantic/channels/channel_map.hpp>
#include <frantic/logging/progress_logger.hpp>
#include <frantic/particles/streams/particle_istream.hpp>

#include "geometry_meshing_parameters_interface.hpp"
#include "geometry_source_interface.hpp"

#include <frost/channel_names.hpp>
#include <frost/material_mode.hpp>
#include <frost/utility.hpp>

#include <boost/math/constants/constants.hpp>

class get_shape_number_strategy {
  public:
    virtual ~get_shape_number_strategy();
    virtual boost::int32_t get_shape_number( boost::int64_t particleNumber, const std::vector<char>& particle ) = 0;
};

class use_first_shape_strategy : public get_shape_number_strategy {
  public:
    virtual boost::int32_t get_shape_number( boost::int64_t particleNumber, const std::vector<char>& particle );
};

class cycle_shapes_strategy : public get_shape_number_strategy {
    boost::int32_t shapeCount;
    cycle_shapes_strategy();

  public:
    cycle_shapes_strategy( boost::int32_t shapeCount );
    virtual boost::int32_t get_shape_number( boost::int64_t particleNumber, const std::vector<char>& particle );
};

class get_random_shape_by_id_strategy : public get_shape_number_strategy {
    boost::int32_t m_seed;
    boost::int32_t m_count;
    frantic::channels::channel_const_cvt_accessor<boost::int32_t> m_idAcc;
    get_random_shape_by_id_strategy();

  public:
    get_random_shape_by_id_strategy( frantic::channels::channel_const_cvt_accessor<boost::int32_t>& idAcc,
                                     boost::int32_t shapeCount, boost::int32_t seed );
    virtual boost::int32_t get_shape_number( boost::int64_t particleNumber, const std::vector<char>& particle );
};

class use_shapeindex_channel_strategy : public get_shape_number_strategy {
    frantic::channels::channel_const_cvt_accessor<boost::int32_t> m_shapeIndexAcc;
    boost::int32_t m_count;
    use_shapeindex_channel_strategy();

  public:
    use_shapeindex_channel_strategy( frantic::channels::channel_const_cvt_accessor<boost::int32_t> shapeIndexAcc,
                                     boost::int32_t shapeCount );
    virtual boost::int32_t get_shape_number( boost::int64_t particleNumber, const std::vector<char>& particle );
};

typedef std::map<boost::uint16_t, boost::uint16_t> geometry_material_id_map_t;

struct material_id_mode_info {
    typedef boost::uint16_t material_id_t;
    int materialMode;
    std::vector<int> meshEndMaterialId;
    std::vector<geometry_material_id_map_t> geometryMaterialIdMaps;
    material_id_t undefinedMaterialID;

    material_id_mode_info( int materialMode, material_id_t undefinedMaterialID )
        : materialMode( materialMode )
        , undefinedMaterialID( undefinedMaterialID ) {}
};

struct mesh_and_accessors {
    frantic::geometry::trimesh3 mesh;
    std::vector<frantic::geometry::trimesh3_vertex_channel_general_accessor> vertexChannels;
    std::vector<frantic::geometry::trimesh3_face_channel_general_accessor> faceChannels;
};

class instance_geometry_factory {
  public:
    virtual ~instance_geometry_factory() {}
    virtual std::size_t get_keep_vertex_channel_number() { return std::numeric_limits<std::size_t>::max(); }
    virtual mesh_and_accessors* get_instance_geometry( double /*timeInSeconds*/, std::size_t /*shapeNumber*/ ) {
        return 0;
    }
    virtual std::size_t get_shape_count() const { return 0; }
    virtual bool is_animated() const { return false; }
};

/**
 * Function object for getting the transform of a mesh instance of a particle.
 */
class instance_info_accessor {
  public:
    instance_info_accessor( const frantic::channels::channel_map& channelMap,
                            boost::shared_ptr<frost::geometry_meshing_parameters_interface> params );

    frantic::graphics::transform4f get_transform( const std::vector<char>& particle ) const;

  private:
    const frost::geometry_orientation_mode::option m_orientationMethod;
    const frantic::graphics::vector3f m_geometryLookAtPosition;
    const frantic::graphics::coordinate_system::option m_coordinateSystem;
    const frantic::graphics::vector3f m_geometryLookAtOrientation;
    const frantic::graphics::vector3f m_geometryOrientation;
    const frantic::graphics::vector3f m_orientationDivergenceAxis;
    const bool m_orientationRestrictDivergenceAxis;
    const frost::geometry_orientation_divergence_axis_space::option m_geometryOrientationDivergenceAxisSpace;

    const float m_divergence;
    const float m_divergenceX2;
    const float m_divergenceRad;

    frantic::channels::channel_const_cvt_accessor<boost::int32_t> m_idAcc;
    frantic::channels::channel_const_cvt_accessor<frantic::graphics::vector3f> m_vecAcc;
    frantic::channels::channel_const_cvt_accessor<frantic::graphics::vector4f> m_orientationAcc;
    frantic::channels::channel_const_cvt_accessor<frantic::graphics::vector3f> m_positionAcc;
    frantic::channels::channel_const_cvt_accessor<frantic::graphics::vector3f> m_radiusXYZAcc;
    frantic::channels::channel_const_cvt_accessor<float> m_radiusAcc;
    frantic::channels::channel_const_cvt_accessor<boost::int32_t> m_shapeIndexAcc;
};

void prepare_mesh_material_channel( frantic::geometry::trimesh3& mesh, frost::material_mode::option materialMode,
                                    int shapeNumber, const std::vector<int>& meshEndMaterialId,
                                    const std::vector<geometry_material_id_map_t>& geometryMaterialIdMaps,
                                    const boost::uint16_t undefinedMaterialID );

std::size_t get_instance_geometry_shape_count( frost::geometry_type::option geometryType,
                                               boost::shared_ptr<geometry_source_interface> customGeometrySource );

boost::shared_ptr<instance_geometry_factory> get_instance_geometry_factory(
    frost::geometry_type::option geometryType, const frantic::channels::channel_map& vertexChannelMap,
    const frantic::channels::channel_map& faceChannelMap, const material_id_mode_info& materialModeInfo,
    boost::shared_ptr<geometry_source_interface> customGeometrySource, frost::hard_edge_type::option hardEdgeType,
    const frantic::graphics::transform4f& coordinateSystem = frantic::graphics::transform4f::identity() );

void geometry_convert_particles_to_trimesh3( boost::shared_ptr<frantic::particles::streams::particle_istream> pin,
                                             const frantic::channels::channel_propagation_policy& cpp,
                                             boost::shared_ptr<frost::geometry_meshing_parameters_interface> params,
                                             boost::shared_ptr<geometry_source_interface> sceneGeometrySource,
                                             material_id_mode_info& materialModeInfo,
                                             frantic::geometry::trimesh3& outMesh, boost::int64_t& outParticleCount,
                                             frantic::logging::progress_logger& progressLogger );
