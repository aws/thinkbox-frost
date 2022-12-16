// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
#include "stdafx.h"

#include <frost/frost.hpp>
#include <frost/frost_parameters.hpp>
#include <frost/geometry_meshing.hpp>
#include <frost/geometry_meshing_parameters_interface.hpp>
#include <frost/geometry_source_interface.hpp>

#include "UnitTests/geometry_meshing_parameters.hpp"

#include <frantic/particles/particle_array.hpp>
#include <frantic/particles/streams/empty_particle_istream.hpp>
#include <frantic/particles/streams/particle_array_particle_istream.hpp>

#include <boost/make_shared.hpp>
#include <boost/shared_array.hpp>

TEST( frost, no_particles ) {
    using namespace frantic::graphics;
    using namespace frantic::particles::streams;

    frantic::geometry::trimesh3 mesh;

    std::size_t particleCount;

    frost_parameters params;

    frantic::channels::channel_map channelMap;
    channelMap.define_channel<vector3f>( _T("Position") );
    channelMap.define_channel<float>( _T("Radius") );
    channelMap.end_channel_definition();

    particle_istream_ptr pin( new empty_particle_istream( channelMap, channelMap ) );

    frantic::logging::null_progress_logger progress;

    ASSERT_NO_THROW( frost::build_trimesh3( mesh, particleCount, params, pin, progress ) );
}

using namespace frost;
// The following classes are used to mock a geometry_meshing_parameters object, which are defined
// in the FrostMY project as maya_geometry_meshing_parameters, and in FrostMax as
// max3d_geometry_meshing_parameters

class geometry_source : public geometry_source_interface {
    std::size_t get_geometry_count() { return 0; }
};

class mesh_orientation_test : public testing::TestWithParam<frantic::graphics::coordinate_system::option> {
  public:
    boost::shared_ptr<geometry_meshing_parameters> params;
    frantic::geometry::trimesh3 mesh;
    boost::int64_t particleCount;
    frantic::logging::null_progress_logger progress;
    frantic::channels::channel_propagation_policy cpp;
    boost::shared_ptr<geometry_source> geometrySource;
    material_id_mode_info materialModeInfo;
    boost::shared_ptr<frantic::particles::streams::particle_array_particle_istream> pin;

    mesh_orientation_test()
        : mesh()
        , particleCount()
        , progress()
        , cpp()
        , materialModeInfo( 0, 0 ) {
        using namespace frantic::graphics;
        using namespace frantic::particles;
        using namespace frantic::particles::streams;

        params = boost::make_shared<geometry_meshing_parameters>();
        geometrySource = boost::shared_ptr<geometry_source>();

        params->set_geometry_type( geometry_type::plane );
        params->set_geometry_orientation_mode( geometry_orientation_mode::look_at );
        params->set_coordinate_system( TestWithParam<frantic::graphics::coordinate_system::option>::GetParam() );

        vector3f position0( 1.f, 0.f, 0.f );
        vector3f position1( 0.f, 0.f, -1.f );
        vector3f position2( -1.f, 0.f, 0.f );
        vector3f position3( 0.f, 0.f, 1.f );
        float radius( 1.f );

        frantic::channels::channel_map channelMap;
        channelMap.define_channel<vector3f>( _T( "Position" ) );
        channelMap.define_channel<float>( _T( "Radius" ) );
        channelMap.end_channel_definition();

        frantic::channels::channel_accessor<vector3f> posAcc = channelMap.get_accessor<vector3f>( _T( "Position" ) );
        frantic::channels::channel_accessor<float> radAcc = channelMap.get_accessor<float>( _T( "Radius" ) );

        particle_array particleArray( channelMap );

        boost::shared_array<char> particle( new char[channelMap.structure_size()] );
        posAcc.get( particle.get() ) = position0;
        radAcc.get( particle.get() ) = radius;
        particleArray.push_back( particle.get() );
        posAcc.get( particle.get() ) = position1;
        particleArray.push_back( particle.get() );
        posAcc.get( particle.get() ) = position2;
        particleArray.push_back( particle.get() );
        posAcc.get( particle.get() ) = position3;
        particleArray.push_back( particle.get() );

        pin = boost::make_shared<particle_array_particle_istream>( particleArray );

        geometry_convert_particles_to_trimesh3( pin, cpp, params, geometrySource, materialModeInfo, mesh, particleCount,
                                                progress );
    }
};

TEST_P( mesh_orientation_test, Check_Normals ) {
    using namespace frantic::graphics;
    using namespace frantic::particles;
    using namespace frantic::particles::streams;

    // params->set_coordinate_system( GetParam() );
    // geometry_convert_particles_to_trimesh3( pin, cpp, params, geometrySource, materialModeInfo, mesh, particleCount,
    //                                         progress );

    /**The output mesh in maya will consist of 4 planes centered around the origin.
     * The expected normals of these planes will be (The negative of each particles' position):
     * [-1, 0, 0]
     * [0, 0, 1]
     * [1, 0, 0]
     * [0, 0, -1]
     **/

    // Four particles will produce 4 planes, each with 2 faces.
    // We grab the first face from each plane to compare normals.
    vector3f mesh0_normal = mesh.calculate_face_normal( 0 );
    vector3f mesh1_normal = mesh.calculate_face_normal( 2 );
    vector3f mesh2_normal = mesh.calculate_face_normal( 4 );
    vector3f mesh3_normal = mesh.calculate_face_normal( 6 );

    vector3f expected_mesh0_normal( -1.f, 0.f, 0.f );
    vector3f expected_mesh1_normal( 0.f, 0.f, 1.f );
    vector3f expected_mesh2_normal( 1.f, 0.f, 0.f );
    vector3f expected_mesh3_normal( 0.f, 0.f, -1.f );

    ASSERT_TRUE( mesh0_normal == expected_mesh0_normal );
    ASSERT_TRUE( mesh1_normal == expected_mesh1_normal );
    ASSERT_TRUE( mesh2_normal == expected_mesh2_normal );
    ASSERT_TRUE( mesh3_normal == expected_mesh3_normal );
}

INSTANTIATE_TEST_CASE_P( Max_and_Maya_Mesh_Tests, mesh_orientation_test,
                         testing::Values( frantic::graphics::coordinate_system::option::right_handed_yup,
                                          frantic::graphics::coordinate_system::option::right_handed_zup ) );
