# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
cmake_minimum_required( VERSION 3.20 FATAL_ERROR )

project( Frost )

find_package( thinkboxcmlibrary REQUIRED )
include( ThinkboxCMLibrary )
include( PrecompiledHeader )

option( BUILD_UNIT_TESTS "Build unit tests" OFF )

add_library( frost STATIC )

target_include_directories( frost PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
	$<INSTALL_INTERFACE:include>
)

file( GLOB_RECURSE H_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    "frost/*.h"
    "frost/*.hpp"
)

file( GLOB_RECURSE CXX_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    "src/*.cpp"
)

target_sources( frost PRIVATE
    stdafx.cpp
    stdafx.h
    ${H_FILES}
    ${CXX_FILES}
)

# The Conan version of Boost was built with this, and it changes the library names.
# As a result, we need to set this to tell Boost to look for the right libraries to
# link against.
target_compile_definitions( frost PUBLIC BOOST_AUTO_LINK_SYSTEM )

find_package( thinkboxlibrary REQUIRED )
find_package( Boost REQUIRED )
find_package( OpenEXR REQUIRED )
find_package( TBB REQUIRED )
find_package( tinyxml2 REQUIRED )
find_package( ZLIB REQUIRED )

target_include_directories( frost PUBLIC ${thinkboxlibrary_INCLUDE_DIRS} )
target_include_directories( frost PUBLIC ${Boost_INCLUDE_DIRS} )
target_include_directories( frost PUBLIC ${OpenEXR_INCLUDE_DIRS} )
target_include_directories( frost PUBLIC ${TBB_INCLUDE_DIRS} )
target_include_directories( frost PUBLIC ${tinyxml2_INCLUDE_DIRS} )
target_include_directories( frost PUBLIC ${ZLIB_INCLUDE_DIRS} )

target_link_libraries( frost INTERFACE thinkboxlibrary::thinkboxlibrary )
target_link_libraries( frost INTERFACE Boost::Boost )
target_link_libraries( frost INTERFACE OpenEXR::OpenEXR )
target_link_libraries( frost INTERFACE TBB::tbb )
target_link_libraries( frost INTERFACE tinyxml2::tinyxml2 )
target_link_libraries( frost INTERFACE ZLIB::ZLIB )

frantic_default_source_groups( frost HEADERDIR frost SOURCEDIR src )
frantic_common_platform_setup( frost )

add_precompiled_header( frost stdafx.h SOURCE_CXX stdafx.cpp )

if( BUILD_UNIT_TESTS )
	add_subdirectory( UnitTests )
endif()

install( DIRECTORY frost
		 DESTINATION include
		 FILES_MATCHING PATTERN "*.hpp"
)
install( TARGETS frost
		 RUNTIME DESTINATION bin
		 LIBRARY DESTINATION lib
		 ARCHIVE DESTINATION lib
)
