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

class geometry_source : public geometry_source_interface {
    std::size_t get_geometry_count() { return 1; }
};

TEST( frost, radiusXYZTest ) {
    using namespace frantic::graphics;
    using namespace frantic::particles::streams;

    boost::shared_ptr<geometry_meshing_parameters> params;
    frantic::geometry::trimesh3 mesh;
    boost::int64_t particleCount;
    frantic::logging::null_progress_logger progress;
    frantic::channels::channel_propagation_policy cpp;
    boost::shared_ptr<geometry_source> geometrySource;
    material_id_mode_info materialModeInfo( 0, 0 );
    boost::shared_ptr<frantic::particles::streams::particle_array_particle_istream> pin;

    params = boost::make_shared<geometry_meshing_parameters>();
    geometrySource = boost::shared_ptr<geometry_source>();

    params->set_geometry_type( geometry_type::box );
    params->set_geometry_orientation_mode( geometry_orientation_mode::specify );

    float inX = 1.f;
    float inY = 2.f;
    float inZ = 3.f;
    vector3f position0( 0.f, 0.f, 0.f );
    float radius( 1.f );
    vector3f radiusXYZ( inX, inY, inZ );

    frantic::channels::channel_map channelMap;
    channelMap.define_channel<vector3f>( _T("Position") );
    channelMap.define_channel<float>( _T("Radius") );
    channelMap.define_channel<vector3f>( _T("RadiusXYZ") );
    channelMap.end_channel_definition();

    frantic::channels::channel_accessor<vector3f> posAcc = channelMap.get_accessor<vector3f>( _T("Position") );
    frantic::channels::channel_accessor<float> radAcc = channelMap.get_accessor<float>( _T("Radius") );
    frantic::channels::channel_accessor<vector3f> radXYZAcc = channelMap.get_accessor<vector3f>( _T("RadiusXYZ") );

    frantic::particles::particle_array particleArray( channelMap );

    boost::shared_array<char> particle( new char[channelMap.structure_size()] );
    posAcc.get( particle.get() ) = position0;
    radAcc.get( particle.get() ) = radius;
    radXYZAcc.get( particle.get() ) = radiusXYZ;
    particleArray.push_back( particle.get() );

    pin = boost::make_shared<particle_array_particle_istream>( particleArray );

    frantic::geometry::trimesh3 outMesh;

    boost::int64_t outParticleCount;

    geometry_convert_particles_to_trimesh3( pin, cpp, params, geometrySource, materialModeInfo, outMesh,
                                            outParticleCount, progress );
    boundbox3f outBox = outMesh.compute_bound_box();
    boundbox3f expectedOutput( inX, inX, inY, inY, inZ, inZ );

    ASSERT_FLOAT_EQ( outBox.xminimum(), expectedOutput.xminimum() * -1 );
    ASSERT_FLOAT_EQ( outBox.yminimum(), expectedOutput.yminimum() * -1 );
    ASSERT_FLOAT_EQ( outBox.zminimum(), expectedOutput.zminimum() * -1 );
    ASSERT_FLOAT_EQ( outBox.xmaximum(), expectedOutput.xmaximum() );
    ASSERT_FLOAT_EQ( outBox.ymaximum(), expectedOutput.ymaximum() );
    ASSERT_FLOAT_EQ( outBox.zmaximum(), expectedOutput.zmaximum() );
}

TEST( frost, radiusTest ) {
    using namespace frantic::graphics;
    using namespace frantic::particles::streams;

    boost::shared_ptr<geometry_meshing_parameters> params;
    frantic::geometry::trimesh3 mesh;
    boost::int64_t particleCount;
    frantic::logging::null_progress_logger progress;
    frantic::channels::channel_propagation_policy cpp;
    boost::shared_ptr<geometry_source> geometrySource;
    material_id_mode_info materialModeInfo( 0, 0 );
    boost::shared_ptr<frantic::particles::streams::particle_array_particle_istream> pin;

    params = boost::make_shared<geometry_meshing_parameters>();
    geometrySource = boost::shared_ptr<geometry_source>();

    params->set_geometry_type( geometry_type::box );
    params->set_geometry_orientation_mode( geometry_orientation_mode::specify );

    vector3f position0( 0.f, 0.f, 0.f );
    float radius( 12.f );

    frantic::channels::channel_map channelMap;
    channelMap.define_channel<vector3f>( _T("Position") );
    channelMap.define_channel<float>( _T("Radius") );
    channelMap.end_channel_definition();

    frantic::channels::channel_accessor<vector3f> posAcc = channelMap.get_accessor<vector3f>( _T("Position") );
    frantic::channels::channel_accessor<float> radAcc = channelMap.get_accessor<float>( _T("Radius") );

    frantic::particles::particle_array particleArray( channelMap );

    boost::shared_array<char> particle( new char[channelMap.structure_size()] );
    posAcc.get( particle.get() ) = position0;
    radAcc.get( particle.get() ) = radius;
    particleArray.push_back( particle.get() );

    pin = boost::make_shared<particle_array_particle_istream>( particleArray );

    frantic::geometry::trimesh3 outMesh;

    boost::int64_t outParticleCount;

    geometry_convert_particles_to_trimesh3( pin, cpp, params, geometrySource, materialModeInfo, outMesh,
                                            outParticleCount, progress );
    boundbox3f outBox = outMesh.compute_bound_box();
    boundbox3f expectedOutput( radius, radius, radius, radius, radius, radius );

    ASSERT_FLOAT_EQ( outBox.xminimum(), expectedOutput.xminimum() * -1 );
    ASSERT_FLOAT_EQ( outBox.yminimum(), expectedOutput.yminimum() * -1 );
    ASSERT_FLOAT_EQ( outBox.zminimum(), expectedOutput.zminimum() * -1 );
    ASSERT_FLOAT_EQ( outBox.xmaximum(), expectedOutput.xmaximum() );
    ASSERT_FLOAT_EQ( outBox.ymaximum(), expectedOutput.ymaximum() );
    ASSERT_FLOAT_EQ( outBox.zmaximum(), expectedOutput.zmaximum() );
}

TEST( frost, invalidRadiusXYZAccessor ) {
    using namespace frantic::graphics;

    frantic::channels::channel_const_cvt_accessor<vector3f> radXYZAcc;

    ASSERT_FALSE( radXYZAcc.is_valid() );
}